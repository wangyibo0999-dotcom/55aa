from pathlib import Path
import json,ast,shlex
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('interface')
OUT=output_dir('reconstruction');OUT.mkdir(exist_ok=True)
ns={'np':np,'shlex':shlex,'aa':dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))}
tree=ast.parse((SCRIPTS/'02_map_reference.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'defs','exec'),ns)
parse,ca,fit=ns['parse'],ns['ca'],ns['fit'];s=json.loads((ROOT/'spatial_audit.json').read_text());c=json.loads((ROOT/'controls_and_local_fit.json').read_text());a=json.loads((ROOT/'actual_interface.json').read_text())
base=RESULTS/'predictions';ref=parse(Path(s['reference_1X9T']));rca=ca(ref['A']);pep=[v for v in ref['B'] if 10<=v['auth']<=19];px=np.array([v['xyz'] for v in pep]);pcas=np.array([next(v['xyz'] for v in pep if v['auth']==i and v['name']=='CA') for i in range(10,20)])
mp={};i=j=0
for x,y in zip(s['alignment_reference_gapped'],s['alignment_target_gapped']):
 i+=x!='-';j+=y!='-'
 if x!='-' and y!='-':mp[i]=(j,x==y)
loc=[i for i,v in rca.items() if i in mp and mp[i][1] and np.linalg.norm(px-v['xyz'],axis=1).min()<=15]
# Select the existing target penton with the smallest worst-chain local reference-site RMSD.
cs=[m for m in c['models'] if m['kind']=='new55'];chosen=min(cs,key=lambda m:max(v['local_site_CA_RMSD_A'] for v in m['site_fits']))['sample'];ch=parse(base/f'sample_{chosen}.cif')
seq=json.loads((base/'inputs_exported.json').read_text())[0]['sequences'][1]['proteinChain']['sequence'];inv={v:k for k,v in ns['aa'].items()}
lines=['REMARK 900 REFERENCE POSITION GUIDES ONLY - NOT A COMPLETE ADAPTER MODEL',f'REMARK 900 PENTON A-E FROM ORIGINAL 55AA SAMPLE {chosen}; UNMODIFIED COORDINATES', 'REMARK 900 GUIDE CHAINS V-Z: CA-ONLY 1X9T PEPTIDE POSITIONS, MAPPED TO 43-52', 'REMARK 900 FIVE POSSIBLE SITES SHOWN; DOES NOT SPECIFY THREE-ADAPTER OCCUPANCY', 'REMARK 900 GUIDE B FACTORS ZERO ARE PLACEHOLDERS, NOT CONFIDENCE SCORES']
serial=1
for chain in 'ABCDE':
 for v in ch[chain]:
  x,y,z=v['xyz'];name=v['name'];lines.append(f'ATOM  {serial:5d} {name:>4s} {inv[v["aa"]]:>3s} {chain}{v["pos"]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{v["b"]:6.2f}          {name[0]:>2s}');serial+=1
 lines.append('TER')
guides=[]
for chain,gc in zip('ABCDE','VWXYZ'):
 cc=ca(ch[chain]);x=np.array([rca[i]['xyz'] for i in loc]);y=np.array([cc[mp[i][0]]['xyz'] for i in loc]);R,t,e=fit(x,y);positions=pcas@R+t
 for n,pos in enumerate(positions,43):
  x,y,z=pos;lines.append(f'ATOM  {serial:5d}   CA {inv[seq[n-1]]:>3s} {gc}{n:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}           C');serial+=1
 lines.append('TER');guides.append({'primary_penton':chain,'guide_chain':gc,'reference_local_CA_RMSD_A':e,'sequence_positions':list(range(43,53)),'CA_coordinates':positions.tolist()})
lines.append('END');(OUT/'penton_with_reference_CA_guides_NOT_full55_model.pdb').write_text('\n'.join(lines)+'\n')
(OUT/'guide_manifest.json').write_text(json.dumps({'source_penton_sample':chosen,'reference':s['reference_1X9T'],'fit_positions':loc,'guides':guides,'status':'local_reference_CA_guides_only; no full adapter rebuilding, loop closure, optimization or validation'},indent=2))
print('Reference guide sample:',chosen)
