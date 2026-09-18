from pathlib import Path
import ast,json,shlex,collections
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('interface')
ns={'np':np,'shlex':shlex,'aa':dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))}
tree=ast.parse((SCRIPTS/'02_map_reference.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'defs','exec'),ns);parse=ns['parse'];ca=ns['ca']
BASE=RESULTS/'predictions'
rows=[];samples=[];sets={4:[],6:[]}
for sample in range(5):
 chains=parse(BASE/f'sample_{sample}.cif');conf=json.loads((BASE/f'confidence_{sample}.json').read_text());s={'sample':sample,'clashes_lt2_atom_pairs':0,'adapters':[]};sample_sets={4:set(),6:set()}
 for ac in 'FGH':
  atoms=chains[ac];x=np.array([a['xyz'] for a in atoms]);apos=np.array([a['pos'] for a in atoms]);am=ca(atoms);adrows=[]
  for pc in 'ABCDE':
   pat=chains[pc];y=np.array([a['xyz'] for a in pat]);ppos=np.array([a['pos'] for a in pat]);pm=ca(pat);d=np.linalg.norm(x[:,None]-y[None,:],axis=-1);ii,jj=np.where(d<=6)
   pairs={}
   for i,j in zip(ii,jj):
    key=(int(apos[i]),int(ppos[j]));pairs[key]=min(pairs.get(key,999),float(d[i,j]))
   for (ap,pp),dist in pairs.items():
    rr={'sample':sample,'adapter_chain':ac,'penton_chain':pc,'adapter_pos':ap,'penton_pos':pp,'adapter_aa':am[ap]['aa'],'penton_aa':pm[pp]['aa'],'min_heavy_A':dist,'adapter_CA_plddt':am[ap]['b'],'penton_CA_plddt':pm[pp]['b'],'pair_iptm':conf['chain_pair_iptm'][ord(ac)-65][ord(pc)-65]}
    rows.append(rr);adrows.append(rr);sample_sets[6].add((ap,pp))
    if dist<=4:sample_sets[4].add((ap,pp))
   s['clashes_lt2_atom_pairs']+=int((d<2).sum())
  contacted=sorted(set(v['adapter_pos'] for v in adrows));near4=[v for v in adrows if v['min_heavy_A']<=4]
  s['adapters'].append({'chain':ac,'contact_observations_6A':len(adrows),'contact_observations_4A':len(near4),'contact_adapter_residues':contacted,'min_penton_CA_plddt':min([v['penton_CA_plddt'] for v in adrows],default=None),'max_pair_iptm':max([v['pair_iptm'] for v in adrows],default=None)})
 for cutoff in (4,6):sets[cutoff].append(sample_sets[cutoff])
 samples.append(s)
 print('sample',sample,'contacts',len(sample_sets[4]),len(sample_sets[6]),'clashes',s['clashes_lt2_atom_pairs'],flush=True)
consensus=[]
for key in sorted(set.union(*sets[6])):
 vv=[v for v in rows if (v['adapter_pos'],v['penton_pos'])==key]
 c={'adapter_pos':key[0],'penton_pos':key[1],'adapter_aa':vv[0]['adapter_aa'],'penton_aa':vv[0]['penton_aa'],'samples_6A':sorted(set(v['sample'] for v in vv)),'samples_4A':sorted(set(v['sample'] for v in vv if v['min_heavy_A']<=4)),'per_sample':[]}
 for sample in range(5):
  v=[r for r in vv if r['sample']==sample];q=[r for r in v if 2<=r['min_heavy_A']<=4]
  c['per_sample'].append({'sample':sample,'adapters_4A':sorted(set(r['adapter_chain'] for r in v if r['min_heavy_A']<=4)),'adapters_2_to_4A':sorted(set(r['adapter_chain'] for r in q)),'min_distance':min([r['min_heavy_A'] for r in v],default=None),'penton_CA_plddt_mean':float(np.mean([r['penton_CA_plddt'] for r in v])) if v else None})
 c['samples_with_at_least2_adapters_2_to_4A']=[v['sample'] for v in c['per_sample'] if len(v['adapters_2_to_4A'])>=2]
 consensus.append(c)
consensus.sort(key=lambda c:(len(c['samples_with_at_least2_adapters_2_to_4A']),len(c['samples_4A']),len(c['samples_6A'])),reverse=True)
jaccard={}
for cutoff in (4,6):
 jaccard[cutoff]=[[len(a&b)/len(a|b) if a|b else 0 for b in sets[cutoff]] for a in sets[cutoff]]
out={'description':'Actual-contact consistency independent of chain labels and reference contact positions; model samples are same-seed diffusion samples, not independent experimental replicates. 2 A is a severe-short-contact flag, not a complete physical validation.','samples':samples,'contact_rows':rows,'consensus':consensus,'pair_jaccard':jaccard}
(ROOT/'actual_interface.json').write_text(json.dumps(out,indent=2))
print('pairs present all5 at4A',sum(len(c['samples_4A'])==5 for c in consensus))
print('pairs present >=3samples >=2adapters with minimum 2-4A',sum(len(c['samples_with_at_least2_adapters_2_to_4A'])>=3 for c in consensus))
print('top consensus',json.dumps(consensus[:8],indent=2))
