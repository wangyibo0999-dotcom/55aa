from pathlib import Path
import ast,json,shlex
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('interface')
ns={'np':np,'shlex':shlex,'aa':dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))}
tree=ast.parse((SCRIPTS/'02_map_reference.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'defs','exec'),ns)
parse=ns['parse'];ca=ns['ca'];fit=ns['fit']
r=json.loads((ROOT/'spatial_audit.json').read_text());ref=parse(Path(r['reference_1X9T']));rca=ca(ref['A']);pep=[a for a in ref['B'] if 10<=a['auth']<=19];px=np.array([a['xyz'] for a in pep]);pca=np.array([next(a['xyz'] for a in pep if a['auth']==i and a['name']=='CA') for i in range(10,20)])
map0={};i=j=0
for a,b in zip(r['alignment_reference_gapped'],r['alignment_target_gapped']):
 i+=a!='-';j+=b!='-'
 if a!='-' and b!='-':map0[i]=(j,a==b)
# Local comparison uses reference CAs within 15 A of any resolved reference peptide heavy atom.
loc=[ri for ri,a in rca.items() if ri in map0 and map0[ri][1] and np.linalg.norm(px-a['xyz'],axis=1).min()<=15]
# Numerical control of proper rigid rotation, translation and reflections.
x=np.array([[0.,0,0],[1,2,0],[0,1,3],[2,0,1],[3,1,2]])
rot=np.array([[0.,-1,0],[1,0,0],[0,0,1]]);y=x@rot+np.array([9.,3.,-4.]);rr,tt,err=fit(x,y);assert err<1e-10 and np.linalg.det(rr)>0.999999
f0r,f0t,f0e=fit(pca,pca);assert f0e<1e-10
out={'numerical_rigid_fit_RMSD':err,'reference_self_fit_RMSD':f0e,'local_reference_CA_count':len(loc),'models':[]}
files=[('new55',i,RESULTS/'predictions'/f'sample_{i}.cif') for i in range(5)]
oldbase=(DATA/'reference52_original')
files += [('reference52',i,next(oldbase.rglob(f'*seed_64438_sample_{i}.cif'))) for i in range(5)]
for kind,sample,f in files:
 ch=parse(f);cc={c:ca(a) for c,a in ch.items()};placements={};fits=[]
 for pc in 'ABCDE':
  x=np.array([rca[ri]['xyz'] for ri in loc]);y=np.array([cc[pc][map0[ri][0]]['xyz'] for ri in loc]);R,t,e=fit(x,y);placements[pc]=pca@R+t;fits.append({'penton':pc,'local_site_CA_RMSD_A':e,'n':len(x)})
 ads=[]
 for ac in 'FGH':
  target=np.array([cc[ac][i]['xyz'] for i in range(43,53)]);vals=[{'penton':pc,'motif_CA_RMSD_A':float(np.sqrt(np.mean(np.sum((p-target)**2,axis=1))))} for pc,p in placements.items()]
  ads.append({'adapter':ac,'best':min(vals,key=lambda a:a['motif_CA_RMSD_A']),'all':vals})
 out['models'].append({'kind':kind,'sample':sample,'site_fits':fits,'adapters':ads})
 print(kind,sample,'site fit',round(max(v['local_site_CA_RMSD_A'] for v in fits),2),'motif',[round(v['best']['motif_CA_RMSD_A'],2) for v in ads])
(ROOT/'controls_and_local_fit.json').write_text(json.dumps(out,indent=2))
