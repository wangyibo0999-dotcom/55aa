from pathlib import Path
import ast,sys,json,hashlib
from collections import defaultdict
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('candidate_geometry')
BASE=output_dir('reconstruction')
from _geometry import Atom,atom_map,place_sg,angle,dihedral
tree=ast.parse((SCRIPTS/'12_scan_candidate_pairs.py').read_text())
ns=dict(np=np,atom_map=atom_map,place_sg=place_sg,dihedral=dihedral,angles=np.arange(-180,180,2))
for name in ['vecangle','rot','geom']:
 node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
 exec(ast.unparse(node).replace('md >= 2', 'md >= 2.5'),ns)
def read(path):
 res=defaultdict(list)
 for l in path.read_text().splitlines():
  if not l.startswith('ATOM'):continue
  if l[76:78].strip() in ('H','D') or l[12:16].strip().startswith('H'):continue
  a=Atom(l[12:16].strip(),l[17:20],l[21],int(l[22:26]),np.array([float(l[i:i+8]) for i in (30,38,46)]),0)
  res[a.chain,a.seq].append(a)
 return dict(res)
models=[('minimized',BASE/'template55_source0_FLEX_NTERM_MINIMIZED.pdb')]+[(f'template{i}_unrelaxed',BASE/f'template55_source{i}_UNRELAXED.pdb') for i in range(5)]
candidates=[('V31C_Q457C',31,457,'VAL','GLN'),('E30C_L223C',30,223,'GLU','LEU')]
rows=[]
def write_pdb(res,changes,path):
 lines=['REMARK EXPLORATORY FIXED-BACKBONE CYS GEOMETRY; NOT RELAXED OR VALIDATED']
 serial=0
 for key,atoms in res.items():
  for a in atoms:
   if key in changes and a.name not in ('N','CA','C','O','OXT','CB'):continue
   serial+=1;rn='CYS' if key in changes else a.resname;x,y,z=a.xyz
   lines.append(f'ATOM  {serial:5d} {a.name:>4s} {rn:3s} {a.chain}{a.seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {a.name[0]:>2s}')
  if key in changes:
   serial+=1;x,y,z=changes[key];ch,ri=key
   lines.append(f'ATOM  {serial:5d}   SG CYS {ch}{ri:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}           S')
 lines.append('END');path.write_text('\n'.join(lines)+'\n')
for label,path in models:
 res=read(path);allatoms=[(k,a) for k,aa in res.items() for a in aa];xyz=np.array([a.xyz for k,a in allatoms])
 for cid,ap,pp,ar,pr in candidates:
  for ac in 'FGH':
   for pc in 'ABCDE':
    ak=(ac,ap);pk=(pc,pp);aa=res[ak];pa=res[pk];am=atom_map(aa);pm=atom_map(pa)
    assert aa[0].resname==ar and pa[0].resname==pr
    cb=float(np.linalg.norm(am['CB'].xyz-pm['CB'].xyz));ca=float(np.linalg.norm(am['CA'].xyz-pm['CA'].xyz))
    g=None
    # Triangle inequality: CB-SG + SG-SG + SG-CB cannot exceed 5.92 A.
    if cb<=5.92:
     near=np.minimum(np.linalg.norm(xyz-am['CB'].xyz,axis=1),np.linalg.norm(xyz-pm['CB'].xyz,axis=1))<6
     keep=np.array([k not in (ak,pk) or a.name in ('N','CA','C','O','OXT') for k,a in allatoms])
     g=ns['geom'](pa,aa,xyz[near&keep])
    row=dict(model=label,candidate=cid,adapter_chain=ac,penton_chain=pc,ca_distance=ca,cb_distance=cb,sg_distance_lower_bound=max(0,cb-3.62),geometry=g,passed=bool(g and g['pass_environment']))
    rows.append(row)
    if row['passed']:
     out=ROOT/f'{cid}_{label}_{ac}{pc}_UNRELAXED.pdb'
     write_pdb(res,{ak:g['sg_a'],pk:g['sg_p']},out);row['structure']=str(out)
 print(label,[(cid,sum(r['passed'] for r in rows if r['model']==label and r['candidate']==cid)) for cid,*_ in candidates],flush=True)
(ROOT/'audit.json').write_text(json.dumps(dict(grid_degrees=2,criteria='SG distance 1.8-2.3 A; CB-SG-SG angles 85-125 degrees; chi3 within30 of +/-90; SG-other retained heavy atoms >=2.5 A. Screening cutoffs, not experimental validation.',sources=[dict(label=l,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for l,p in models],observations=rows),indent=2))

assert all(r['geometry']['environment_min'] >= 2.5 for r in rows if r['passed'])
