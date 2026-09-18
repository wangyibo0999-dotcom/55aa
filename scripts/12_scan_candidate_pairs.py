from pathlib import Path
import sys,json,math,hashlib
from collections import defaultdict
import numpy as np
from _geometry import parse_cif,atom_map,place_sg,angle,dihedral
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('candidate_scan')
SRC=(DATA/'reference52_original')
seq52=json.loads((SRC/'inputs.json').read_text())[0]['sequences'][1]['proteinChain']['sequence']
seq55='GYIPEAPRDGQAYVRKDGEWVLLSTFLAKEVAKTKTKKVSRGTFDPVYPYDADNA'
SIDE={('F','A'):'primary',('G','B'):'primary',('H','D'):'primary',('F','B'):'adjacent',('G','C'):'adjacent',('H','E'):'adjacent'}
angles=np.arange(-180,180,10)
def vecangle(u,v):
 return np.degrees(np.arccos(np.clip(np.sum(u*v,axis=-1)/(np.linalg.norm(u,axis=-1)*np.linalg.norm(v,axis=-1)),-1,1)))
def rot(x):return min(abs((x-c+180)%360-180) for c in [-60,60,180])
def geom(pa,aa,env):
 pm=atom_map(pa);am=atom_map(aa)
 ps=np.array([place_sg(pa,float(k)) for k in angles]);ass=np.array([place_sg(aa,float(k)) for k in angles])
 d=np.linalg.norm(ps[:,None]-ass[None,:],axis=-1);ii,jj=np.where((d>=1.8)&(d<=2.3))
 if not len(ii):return None
 p=ps[ii];a=ass[jj];angp=vecangle(pm['CB'].xyz-p,a-p);anga=vecangle(am['CB'].xyz-a,p-a)
 good=(angp>=85)&(angp<=125)&(anga>=85)&(anga<=125);ii=ii[good];jj=jj[good];angp=angp[good];anga=anga[good]
 opts=[]
 for z,(i,j) in enumerate(zip(ii,jj)):
  chi3=dihedral(pm['CB'].xyz,ps[i],ass[j],am['CB'].xyz)
  if abs(abs(chi3)-90)>30:continue
  # Use actual N-CA-CB-SG dihedral rather than assuming placement parameter convention.
  c1=dihedral(pm['N'].xyz,pm['CA'].xyz,pm['CB'].xyz,ps[i]);c2=dihedral(am['N'].xyz,am['CA'].xyz,am['CB'].xyz,ass[j])
  md=float(min(np.linalg.norm(env-ps[i],axis=1).min(),np.linalg.norm(env-ass[j],axis=1).min()))
  score=((d[i,j]-2.03)/.25)**2+((angp[z]-104)/18)**2+((anga[z]-104)/18)**2+((abs(chi3)-90)/30)**2+((rot(c1)+rot(c2))/40)**2
  opts.append(dict(sg_distance=float(d[i,j]),angle_p=float(angp[z]),angle_a=float(anga[z]),chi3=chi3,chi1_p=c1,chi1_a=c2,parameter_p=int(angles[i]),parameter_a=int(angles[j]),environment_min=md,pass_environment=md>=2,score=float(score),sg_p=ps[i].tolist(),sg_a=ass[j].tolist()))
 return min(opts,key=lambda v:(not v['pass_environment'],v['score'])) if opts else None
rows=[];inventory=[]
for file in sorted(SRC.rglob('*sample_*.cif')):
 sample=int(file.stem.rsplit('_',1)[1]);res=parse_cif(file);inventory.append({'path':str(file),'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
 all_atoms=[(key,a) for key,atoms in res.items() for a in atoms]
 xyz=np.array([a.xyz for _,a in all_atoms]); keys=[key for key,_ in all_atoms]
 caches={}
 for ac in 'FGH':
  for pc in 'ABCDE':
   aa=[a for key,atoms in res.items() if key[0]==ac for a in atoms];pa=[a for key,atoms in res.items() if key[0]==pc for a in atoms]
   x=np.array([a.xyz for a in aa]);y=np.array([a.xyz for a in pa]);ds=np.linalg.norm(x[:,None]-y[None,:],axis=-1)
   ai,pi=np.where(ds<=6);pairs=defaultdict(lambda:float('inf'))
   for i,j in zip(ai,pi):
    k=(aa[i].seq,pa[j].seq);pairs[k]=min(pairs[k],float(ds[i,j]))
   for (ap,pp),nearest in sorted(pairs.items()):
    ak=(ac,ap);pk=(pc,pp);a=res[ak];p=res[pk];am=atom_map(a);pm=atom_map(p)
    row=dict(sample=sample,adapter_chain=ac,penton_chain=pc,adapter_position=ap,penton_position=pp,adapter_reference_resname=a[0].resname,penton_resname=p[0].resname,adapter_target_letter=seq55[ap-1],adapter_identity_conserved=seq52[ap-1]==seq55[ap-1],side=SIDE.get((ac,pc),'other'),min_heavy=nearest,ca_distance=float(np.linalg.norm(am['CA'].xyz-pm['CA'].xyz)),adapter_plddt=am['CA'].plddt,penton_plddt=pm['CA'].plddt,geometry=None)
    if 'CB' in am and 'CB' in pm:row['cb_distance']=float(np.linalg.norm(am['CB'].xyz-pm['CB'].xyz))
    row['broad_pass']=a[0].resname not in ('GLY','PRO','CYS') and p[0].resname not in ('GLY','PRO','CYS') and row['ca_distance']<=9 and row.get('cb_distance',99)<=7.5
    if row['broad_pass']:
     # All non-mutated atoms plus retained backbone N/CA/C/O of the two mutated residues.
     near=np.minimum(np.linalg.norm(xyz-am['CB'].xyz,axis=1),np.linalg.norm(xyz-pm['CB'].xyz,axis=1))<6
     keep=np.array([k not in (ak,pk) or atom.name in ('N','CA','C','O','OXT') for k,atom in all_atoms])
     env=xyz[near&keep]
     row['geometry']=geom(p,a,env)
    rows.append(row)
 print('sample',sample,'contact observations',len(rows),'geometry environment passes',sum(bool(r['geometry'] and r['geometry']['pass_environment']) for r in rows),flush=True)
(ROOT/'scan.json').write_text(json.dumps({'reference':inventory,'seq52':seq52,'seq55':seq55,'changes':[(i+1,a,b) for i,(a,b) in enumerate(zip(seq52,seq55)) if a!=b],'observations':rows},indent=2))
groups=defaultdict(list)
for r in rows:groups[(r['adapter_position'],r['penton_position'],r['side'])].append(r)
rank=[]
for k,v in groups.items():
 passed=[r for r in v if r['geometry'] and r['geometry']['pass_environment']]
 if not passed:continue
 samples=[s for s in range(5) if len({r['adapter_chain'] for r in passed if r['sample']==s})==3]
 rank.append(dict(adapter=k[0],penton=k[1],side=k[2],reference_resname=v[0]['adapter_reference_resname'],target_letter=v[0]['adapter_target_letter'],penton_resname=v[0]['penton_resname'],conserved=v[0]['adapter_identity_conserved'],contacts=len(v),passes=len(passed),samples_all3=samples,mean_adapter_plddt=float(np.mean([r['adapter_plddt'] for r in passed])),mean_penton_plddt=float(np.mean([r['penton_plddt'] for r in passed])),mean_score=float(np.mean([r['geometry']['score'] for r in passed]))))
rank.sort(key=lambda r:(len(r['samples_all3']),r['passes'],r['conserved'],-r['mean_score']),reverse=True)
(ROOT/'rank.json').write_text(json.dumps(rank,indent=2));print(json.dumps(rank[:20],indent=2))
