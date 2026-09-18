"""Reference-guided intact-trimer rigid placements. NOT flexible refinement or validation."""
from pathlib import Path
import ast,json,shlex,itertools,hashlib
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('reconstruction')
BASE=RESULTS/'predictions'
ns={'np':np,'shlex':shlex,'aa':dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))}
tree=ast.parse((SCRIPTS/'02_map_reference.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'defs','exec'),ns)
parse,ca,fit=ns['parse'],ns['ca'],ns['fit'];inv={v:k for k,v in ns['aa'].items()}
guides=json.loads((ROOT/'guide_manifest.json').read_text());ps=guides['source_penton_sample'];penton=parse(BASE/f'sample_{ps}.cif');pat=[a for c in 'ABCDE' for a in penton[c]];px=np.array([a['xyz'] for a in pat]);ref={g['primary_penton']:np.array(g['CA_coordinates']) for g in guides['guides']}
inputs=json.loads((BASE/'inputs_exported.json').read_text())[0];seqs=[v['proteinChain']['sequence'] for v in inputs['sequences']]
# Spatial hash avoids constructing a dense whole-complex atom-distance matrix for every pose.
cells={};cellwidth=4.
for i,x in enumerate(px):cells.setdefault(tuple(np.floor(x/cellwidth).astype(int)),[]).append(i)
def overlaps(ax):
 count2=0;count25=0;minimum=999.;ar=set();pr=set()
 for i,x in enumerate(ax):
  cc=np.floor(x/cellwidth).astype(int);ids=[]
  for off in itertools.product((-1,0,1),repeat=3):ids.extend(cells.get(tuple(cc+off),[]))
  if not ids:continue
  ids=np.array(ids);d=np.linalg.norm(px[ids]-x,axis=1);minimum=min(minimum,float(d.min()));count2+=int((d<2).sum());count25+=int((d<2.5).sum())
  if (d<2).any():ar.add(i);pr.update(ids[d<2].tolist())
 return {'atom_pairs_lt2_A':count2,'atom_pairs_lt2_5_A':count25,'minimum_interprotein_heavy_A':minimum,'adapter_atoms_in_lt2_contacts':len(ar),'penton_atoms_in_lt2_contacts':len(pr)}
allrows=[];structures={}
for s in range(5):
 ch=parse(BASE/f'sample_{s}.cif');am={c:ca(ch[c]) for c in 'FGH'}
 for c in 'FGH':assert ''.join(a['aa'] for a in am[c].values())==seqs[1]
 allat=[(c,a) for c in 'FGH' for a in ch[c]];ax=np.array([a['xyz'] for c,a in allat]);motif=np.concatenate([np.array([am[c][i]['xyz'] for i in range(43,53)]) for c in 'FGH'])
 fits=[]
 for perm in itertools.permutations('ABCDE',3):
  target=np.concatenate([ref[c] for c in perm]);R,t,e=fit(motif,target)
  fits.append({'source_adapter_sample':s,'guide_assignment_FGH':''.join(perm),'motif30_CA_RMSD_A':e,'rotation':R.tolist(),'translation':t.tolist(),'individual_motif_RMSD_A':[float(np.sqrt(np.mean(np.sum((motif[k*10:(k+1)*10]@R+t-target[k*10:(k+1)*10])**2,axis=1)))) for k in range(3)]})
 # All 60 scored on location; only best three per source receive detailed clash checks.
 fits.sort(key=lambda r:r['motif30_CA_RMSD_A'])
 for index,r in enumerate(fits):
  if index<3:r['clash_audit']=overlaps(ax@np.array(r['rotation'])+np.array(r['translation']))
  allrows.append(r)
 structures[s]=(allat,ax,motif)
 print('sample',s,'best RMSD',round(fits[0]['motif30_CA_RMSD_A'],3),'assignment',fits[0]['guide_assignment_FGH'],'clashes',fits[0]['clash_audit']['atom_pairs_lt2_A'],flush=True)
# Save one best-location intact-trimer start per original sample. No physical-validity selection claim.
selected=[]
for s in range(5):
 r=min([v for v in allrows if v['source_adapter_sample']==s],key=lambda v:v['motif30_CA_RMSD_A']);allat,ax,motif=structures[s];R=np.array(r['rotation']);t=np.array(r['translation']);moved=ax@R+t
 # Verify one rigid transform preserves all trimer geometry to floating-point precision.
 sampleinds=np.arange(0,len(ax),max(1,len(ax)//40));before=np.linalg.norm(ax[sampleinds,None]-ax[sampleinds][None,:],axis=-1);after=np.linalg.norm(moved[sampleinds,None]-moved[sampleinds][None,:],axis=-1)
 assert np.max(np.abs(before-after))<1e-8
 path=ROOT/f'guided55_source{s}_UNRELAXED.pdb';lines=['REMARK 900 REFERENCE-GUIDED RIGID START ONLY; NOT RELAXED OR VALIDATED',f'REMARK 900 PENTON SOURCE 55AA SAMPLE {ps}; ADAPTER TRIMER SOURCE SAMPLE {s}',f'REMARK 900 FGH MATCHED TO REFERENCE PRIMARY PENTON SITES {r["guide_assignment_FGH"]}',f'REMARK 900 30-CA GUIDE RMSD {r["motif30_CA_RMSD_A"]:.3f} A',f'REMARK 900 CROSS-PROTEIN HEAVY ATOM PAIRS BELOW 2 A: {r["clash_audit"]["atom_pairs_lt2_A"]}', 'REMARK 900 B FACTORS SET TO ZERO: NO CONFIDENCE SCORES FOR THIS NEW POSE'];serial=1
 for chain in 'ABCDE':
  for a in penton[chain]:
   x,y,z=a['xyz'];name=a['name'];lines.append(f'ATOM  {serial:5d} {name:>4s} {inv[a["aa"]]:>3s} {chain}{a["pos"]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {name[0]:>2s}');serial+=1
  lines.append('TER')
 for chain in 'FGH':
  for (c,a),coord in zip(allat,moved):
   if c!=chain:continue
   x,y,z=coord;name=a['name'];lines.append(f'ATOM  {serial:5d} {name:>4s} {inv[a["aa"]]:>3s} {chain}{a["pos"]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}          {name[0]:>2s}');serial+=1
  lines.append('TER')
 lines.append('END');path.write_text('\n'.join(lines)+'\n');r['saved_file']=path.name;r['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();r['rigid_geometry_max_error_A']=float(np.max(np.abs(before-after)));selected.append(r)
(ROOT/'rigid_start_audit.json').write_text(json.dumps({'status':'UNRELAXED diagnostic reference-guided starts; no binding or disulfide validation','penton_source_sample':ps,'all_300_assignments':allrows,'saved_starts':selected,'selection':'lowest 30-CA reference-guide RMSD separately for each source trimer; NOT lowest energy','source_files':[{'path':str(BASE/f'sample_{s}.cif'),'sha256':hashlib.sha256((BASE/f'sample_{s}.cif').read_bytes()).hexdigest()} for s in range(5)]},indent=2))
