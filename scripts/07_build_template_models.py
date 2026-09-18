from pathlib import Path
import json
import pymol2
from pymol import editor
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('reconstruction')
OLD=(DATA/'reference52_original')
SEQ='GYIPEAPRDGQAYVRKDGEWVLLSTFLAKEVAKTKTKKVSRGTFDPVYPYDADNA'
records=[]
for s in range(5):
 source=next(OLD.rglob(f'*seed_64438_sample_{s}.cif'))
 with pymol2.PyMOL() as pm:
  cmd=pm.cmd;cmd.set('retain_order',0);cmd.load(str(source),'complex55');cmd.remove('hydro');before=cmd.get_coords('complex55 and chain A+B+C+D+E',1).copy()
  for chain in 'FGH':
   for pos,res in [(27,'LEU'),(45,'ASP'),(52,'ALA')]:
    cmd.wizard('mutagenesis');w=cmd.get_wizard();w.set_mode(res);cmd.select('mutation_target',f'complex55 and chain {chain} and resi {pos}');w.do_select('mutation_target');cmd.frame(1);w.apply();cmd.set_wizard()
   cmd.remove(f'complex55 and chain {chain} and resi 52 and name OXT')
   segment=cmd.get_model(f'complex55 and chain {chain} and resi 1 and name CA').atom[0].segi
   cmd.fab('ADNA','extension',resi=52,chain=chain,segi=segment,ss=2,hydro=0)
   fitargs=[]
   for atom in ['N','CA','C']:
    fitargs.extend([f'extension and resi 52 and name {atom}',f'complex55 and chain {chain} and resi 52 and name {atom}'])
   cmd.pair_fit(*fitargs)
   cmd.remove('extension and resi 52')
   cmd.create('complex55','complex55 or extension',1,1);cmd.delete('extension')
   cmd.bond(f'complex55 and chain {chain} and resi 52 and name C',f'complex55 and chain {chain} and resi 53 and name N')
   cmd.unpick()
   cmd.sort('complex55')
   got=''.join(cmd.get_fastastr(f'complex55 and chain {chain}').splitlines()[1:]);assert got==SEQ,(chain,got)
  cmd.remove('hydro');cmd.sort('complex55');out=ROOT/f'template55_source{s}_UNRELAXED.pdb';cmd.save(str(out),'complex55',state=1)
  content=out.read_text().splitlines();out.write_text('REMARK 900 TEMPLATE-DERIVED 55AA INITIAL MODEL, UNRELAXED, NOT VALIDATED\nREMARK 900 B FACTORS ZERO; NO NEW MODEL CONFIDENCE\n'+'\n'.join(l[:60]+'  0.00'+l[66:] if l.startswith('ATOM') else l for l in content)+'\n')
  atoms=cmd.get_model('complex55').atom
  records.append({'sample':s,'source':str(source),'output':str(out),'chain_lengths':{c:cmd.count_atoms(f'complex55 and chain {c} and name CA') for c in 'ABCDEFGH'},'sequence55_checked_all3':True,'source_of_new_residues':'PyMOL residue fragments; extended initial phi/psi (-139,135); not validated','mutations':['K27L','N45D','T52A'],'appended':'DNA'})
  print('built',s,records[-1]['chain_lengths'],flush=True)
(ROOT/'template_build.json').write_text(json.dumps(records,indent=2))
