from pathlib import Path
import tarfile,json,hashlib,shlex
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('predictions')
ARCH=(DATA/'55aa_predictions.tar.gz')
aa=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def parse(s):
 headers=[]; out={}
 for l in s.splitlines():
  if l.startswith('_atom_site.'):headers.append(l.strip().split()[0].split('.',1)[1])
  elif l.startswith('ATOM '):
   d=dict(zip(headers,shlex.split(l)))
   c=d['label_asym_id'];r=int(d['label_seq_id']);name=d['label_atom_id'];elem=d['type_symbol']
   if elem in ('H','D'):continue
   out.setdefault(c,[]).append((r,name,aa[d['label_comp_id']],float(d['B_iso_or_equiv']),[float(d['Cartn_'+a]) for a in 'xyz']))
 return out
result={'archive_sha256':hashlib.sha256(ARCH.read_bytes()).hexdigest(),'samples':[]}
with tarfile.open(ARCH) as t:
 inp=json.load(t.extractfile('inputs.json'));result['input']=inp
 (ROOT/'inputs_exported.json').write_text(json.dumps(inp,indent=2))
 for i in range(5):
  base='55aa/seed_79470/predictions/'
  cif=t.extractfile(base+f'55aa_sample_{i}.cif').read();(ROOT/f'sample_{i}.cif').write_bytes(cif)
  chains=parse(cif.decode());conf=json.load(t.extractfile(base+f'55aa_summary_confidence_sample_{i}.json'))
  (ROOT/f'confidence_{i}.json').write_text(json.dumps(conf,indent=2))
  row={'sample':i,'global':{k:v for k,v in conf.items() if not isinstance(v,list)},'chains':{},'interfaces':[],'reference_site_distances':[]}
  for c,atoms in chains.items():
   ca=[a for a in atoms if a[1]=='CA'];seq=''.join(a[2] for a in ca);expected=inp[0]['sequences'][0 if len(ca)>100 else 1]['proteinChain']['sequence']
   row['chains'][c]={'length':len(ca),'sequence_matches_input':seq==expected,'ca_plddt':float(np.mean([a[3] for a in ca])),'tail43_52_plddt':float(np.mean([a[3] for a in ca if 43<=a[0]<=52]))}
  for ac in 'FGH':
   a=chains[ac];x=np.array([v[4] for v in a]);ar=np.array([v[0] for v in a])
   for pc in 'ABCDE':
    p=chains[pc];y=np.array([v[4] for v in p]);pr=np.array([v[0] for v in p])
    ds=np.sqrt(((x[:,None,:]-y[None,:,:])**2).sum(axis=2));ii,jj=np.where(ds<=6)
    row['interfaces'].append({'adapter':ac,'penton':pc,'min_heavy_A':float(ds.min()),'residue_pairs_6A':len(set(zip(ar[ii].tolist(),pr[jj].tolist()))),'atom_pairs_lt2A':int((ds<2).sum()),'pair_iptm':conf['chain_pair_iptm'][ord(ac)-65][ord(pc)-65]})
    for ap,pp in [(45,203),(48,466),(50,233)]:
     xx=x[ar==ap];yy=y[pr==pp];dd=np.sqrt(((xx[:,None,:]-yy[None,:,:])**2).sum(axis=2)).min()
     row['reference_site_distances'].append({'adapter':ac,'penton':pc,'adapter_res':ap,'penton_res':pp,'min_heavy_A':float(dd)})
  result['samples'].append(row)
(ROOT/'audit.json').write_text(json.dumps(result,indent=2))
for r in result['samples']:
 print('SAMPLE',r['sample'],'global',r['global'])
 print('chains',r['chains'])
 for c in 'FGH':
  rows=[v for v in r['interfaces'] if v['adapter']==c]
  print(c,'nearest',min(rows,key=lambda v:v['min_heavy_A']),'total_contact_pairs',sum(v['residue_pairs_6A'] for v in rows))
 print('best ANY-chain reference distances',[(ap,pp,round(min(v['min_heavy_A'] for v in r['reference_site_distances'] if v['adapter_res']==ap and v['penton_res']==pp),3)) for ap,pp in [(45,203),(48,466),(50,233)]])
