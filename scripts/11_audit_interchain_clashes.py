from pathlib import Path
import ast,json,itertools,shlex
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('reconstruction')
# Read PDB without invoking earlier report scripts.
ns={'np':np};tree=ast.parse((SCRIPTS/'10_check_rebuilt_models.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'parser','exec'),ns);parse=ns['parse']
files=[ROOT/'template55_source0_prepared.pdb',ROOT/'template55_source0_LOCAL_MINIMIZED.pdb',ROOT/'template55_source0_FLEX_NTERM_MINIMIZED.pdb']
results=[]
for f in files:
 if not f.exists():continue
 ch=parse(f);grid={};rows=[]
 for chain,aa in ch.items():
  for a in aa:
   a=dict(a,chain=chain);x=a['xyz'];key=np.floor(x/2).astype(int)
   for off in itertools.product((-1,0,1),repeat=3):
    for b in grid.get(tuple(key+off),[]):
     if chain==b['chain']:continue
     dist=float(np.linalg.norm(x-b['xyz']))
     if dist<2:
      group='penton-penton' if chain in 'ABCDE' and b['chain'] in 'ABCDE' else ('adapter-adapter' if chain in 'FGH' and b['chain'] in 'FGH' else 'adapter-penton')
      rows.append({'type':group,'atom1':f'{chain}:{a["pos"]}:{a["name"]}','atom2':f'{b["chain"]}:{b["pos"]}:{b["name"]}','distance_A':dist})
   grid.setdefault(tuple(key),[]).append(a)
 results.append({'file':f.name,'counts_below2A':{k:sum(r['type']==k for r in rows) for k in ['penton-penton','adapter-adapter','adapter-penton']},'counts_below1A':{k:sum(r['type']==k and r['distance_A']<1 for r in rows) for k in ['penton-penton','adapter-adapter','adapter-penton']},'short_contacts':rows})
 print(f.name,results[-1]['counts_below2A'],flush=True)
(ROOT/'all_interchain_short_contacts.json').write_text(json.dumps(results,indent=2))
