from pathlib import Path
import json,sys
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('reconstruction')
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
inp=json.loads((RESULTS/'predictions/inputs_exported.json').read_text())[0];seqs=[r['proteinChain']['sequence'] for r in inp['sequences']]
def parse(f):
 out={}
 for l in f.read_text().splitlines():
  if not l.startswith('ATOM'):continue
  name=l[12:16].strip();element=l[76:78].strip()
  if element in ('H','D') or name.startswith('H'):continue
  out.setdefault(l[21],[]).append({'pos':int(l[22:26]),'name':name,'res':l[17:20].strip(),'xyz':np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])})
 return out
results=[]
for f in sorted(ROOT.glob('template55_source*.pdb')):
 if '_prepared' in f.name:continue
 ch=parse(f);r={'file':f.name,'chains':{},'interfaces':[],'interprotein_atom_pairs_lt2A':0,'peptide_bond_warnings':[]}
 assert set(ch)==set('ABCDEFGH'),(f,list(ch))
 for c,atoms in ch.items():
  ca={a['pos']:a for a in atoms if a['name']=='CA'};seq=''.join(AA[ca[i]['res']] for i in sorted(ca));assert seq==seqs[0 if c in 'ABCDE' else 1],(f,c,seq)
  r['chains'][c]={'length':len(ca),'sequence_verified':True}
  m={(a['pos'],a['name']):a['xyz'] for a in atoms}
  for i in range(1,len(ca)):
   d=float(np.linalg.norm(m[i,'C']-m[i+1,'N']))
   if not 1.1<=d<=1.6:r['peptide_bond_warnings'].append({'chain':c,'residue':i,'C_N_A':d})
 for ac in 'FGH':
  a=ch[ac];x=np.array([v['xyz'] for v in a]);ar=np.array([v['pos'] for v in a])
  for pc in 'ABCDE':
   p=ch[pc];y=np.array([v['xyz'] for v in p]);pr=np.array([v['pos'] for v in p]);ds=np.linalg.norm(x[:,None]-y[None,:],axis=-1)
   r['interfaces'].append({'adapter':ac,'penton':pc,'atom_pairs_lt2A':int((ds<2).sum()),'minimum_A':float(ds.min()),'Y48_H466_min_heavy_A':float(ds[np.ix_(ar==48,pr==466)].min()),'Y50_M233_min_heavy_A':float(ds[np.ix_(ar==50,pr==233)].min())});r['interprotein_atom_pairs_lt2A']+=int((ds<2).sum())
 results.append(r);print(f.name,'clashes',r['interprotein_atom_pairs_lt2A'],'peptide warnings',len(r['peptide_bond_warnings']),flush=True)
(ROOT/'rebuilt_model_checks.json').write_text(json.dumps(results,indent=2))
