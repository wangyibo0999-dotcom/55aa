from pathlib import Path
import ast,json,sys
import numpy as np
from collections import defaultdict
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('candidate_scan')
# Load definitions without re-running the exhaustive scan.
tree=ast.parse((SCRIPTS/'12_scan_candidate_pairs.py').read_text())
ns={'__file__':str(SCRIPTS/'12_scan_candidate_pairs.py')}
keep=[]
for n in tree.body:
 if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef)):keep.append(n)
 elif isinstance(n,ast.Assign) and all(isinstance(t,ast.Name) and t.id in {'ROOT','SRC','seq52','seq55','SIDE','angles'} for t in n.targets):keep.append(n)
 elif isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='sys.path.insert':keep.append(n)
exec(compile(ast.Module(body=keep,type_ignores=[]),'screen_definitions','exec'),ns)
ns['angles']=np.arange(-180,180,5)
geom_node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='geom')
exec(ast.unparse(geom_node).replace('md >= 2', 'md >= 2.5'),ns)
parse=ns['parse_cif'];amap=ns['atom_map'];geom=ns['geom'];place=ns['place_sg'];angle=ns['angle'];dih=ns['dihedral']
scan=json.loads((ROOT/'scan.json').read_text());rank=json.loads((ROOT/'rank.json').read_text())
# Re-evaluate every initially three-copy-passing pair, including substituted sites, at all 15 expected interfaces.
candidates=[r for r in rank if r['samples_all3'] and r['side']!='other']
records=[];checks=[]
for inv in scan['reference']:
 f=Path(inv['path']);s=int(f.stem.rsplit('_',1)[1]);res=parse(f)
 allatoms=[(k,a) for k,atoms in res.items() for a in atoms];xyz=np.array([a.xyz for _,a in allatoms])
 for c in candidates:
  pcs='ABD' if c['side']=='primary' else 'BCE'
  for ac,pc in zip('FGH',pcs):
   ak=(ac,c['adapter']);pk=(pc,c['penton']);a=res[ak];p=res[pk];aM=amap(a);pM=amap(p)
   minimum=float(np.linalg.norm(np.array([v.xyz for v in a])[:,None]-np.array([v.xyz for v in p])[None,:],axis=-1).min())
   ca=float(np.linalg.norm(aM['CA'].xyz-pM['CA'].xyz));cb=float(np.linalg.norm(aM['CB'].xyz-pM['CB'].xyz))
   near=np.minimum(np.linalg.norm(xyz-aM['CB'].xyz,axis=1),np.linalg.norm(xyz-pM['CB'].xyz,axis=1))<6
   keepmask=np.array([k not in (ak,pk) or atom.name in ('N','CA','C','O','OXT') for k,atom in allatoms])
   g=geom(p,a,xyz[near&keepmask]) if ca<=9 and cb<=7.5 and minimum<=6 else None
   if g:
    # Explicit independent recomputation of construction invariants.
    for m,key in [(aM,'a'),(pM,'p')]:
     sg=np.array(g['sg_'+key]);length=float(np.linalg.norm(sg-m['CB'].xyz));ang=angle(m['CA'].xyz,m['CB'].xyz,sg)
     assert abs(length-1.81)<1e-8 and abs(ang-114)<1e-8
    g['chi2_p']=dih(pM['CA'].xyz,pM['CB'].xyz,np.array(g['sg_p']),np.array(g['sg_a']))
    g['chi2_a']=dih(aM['CA'].xyz,aM['CB'].xyz,np.array(g['sg_a']),np.array(g['sg_p']))
   records.append(dict(sample=s,adapter=c['adapter'],penton=c['penton'],side=c['side'],adapter_chain=ac,penton_chain=pc,min_heavy=minimum,ca_distance=ca,cb_distance=cb,adapter_plddt=aM['CA'].plddt,penton_plddt=pM['CA'].plddt,geometry=g))
 print('refined sample',s,flush=True)
summary=[]
for c in candidates:
 rr=[r for r in records if (r['adapter'],r['penton'],r['side'])==(c['adapter'],c['penton'],c['side'])]
 pp=[r for r in rr if r['geometry'] and r['geometry']['pass_environment']]
 summary.append(dict(**c,refined_passes=len(pp),refined_all3=[s for s in range(5) if sum(r['sample']==s for r in pp)==3],refined_mean_score=float(np.mean([r['geometry']['score'] for r in pp])) if pp else None,best_pose_environment_ge_2_5=sum(r['geometry']['environment_min']>=2.5 for r in pp),all_observations=15))
summary.sort(key=lambda c:(len(c['refined_all3']),c['refined_passes'],c['conserved'],-(c['refined_mean_score'] or 999)),reverse=True)
(ROOT/'refined_stricter.json').write_text(json.dumps({'grid_degrees':5,'ranking':summary,'observations':records,'construction_checks':'All evaluated SG atoms: CB-SG=1.81 A and CA-CB-SG=114 degrees, numerical tolerance 1e-8.'},indent=2))
for c in summary:print(c['adapter'],c['penton'],c['side'],c['refined_passes'],c['refined_all3'],round(c['refined_mean_score'],3),c['best_pose_environment_ge_2_5'])
