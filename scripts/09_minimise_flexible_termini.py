from pathlib import Path
import sys,os,json,time,random
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('reconstruction')
os.environ['OPENMM_CPU_THREADS']='4'
import numpy as np
import openmm as mm
from openmm import app,unit
s=int(sys.argv[1]);random.seed(5500+s);np.random.seed(5500+s)
source=ROOT/f'template55_source{s}_LOCAL_MINIMIZED.pdb';fixed=ROOT/f'template55_source{s}_flex_prepared.pdb'
# Remove inherited source confidence values and explicit bonds; reconstruct standard protein bonds.
lines=[l for l in source.read_text().splitlines() if l.startswith('ATOM')]
rows=[]
for l in lines:
 if l[76:78].strip() in ('H','D'):continue
 rows.append(l[:60]+'  0.00'+l[66:])
# Add missing terminal oxygen on all eight free C termini using planar carbonyl geometry.
added=[];serial=max(int(l[6:11]) for l in rows)+1
for chain in 'ABCDEFGH':
 rr=[l for l in rows if l[21]==chain];last=max(int(l[22:26]) for l in rr);r=[l for l in rr if int(l[22:26])==last];m={l[12:16].strip():np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in r}
 if 'OXT' not in m:
  u=(m['CA']-m['C']);u/=np.linalg.norm(u);v=m['O']-m['C'];v/=np.linalg.norm(v);direction=-(u+v);direction/=np.linalg.norm(direction);x,y,z=m['C']+1.25*direction
  line=f'ATOM  {serial:5d}  OXT {r[0][17:20]} {chain}{last:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}           O';rows.append(line);serial+=1;added.append((chain,last))
ordered=[]
for chain in 'ABCDEFGH':ordered.extend([l for l in rows if l[21]==chain]);ordered.append('TER')
fixed.write_text('\n'.join(ordered+['END'])+'\n');pdb=app.PDBFile(str(fixed));ff=app.ForceField('amber14-all.xml','implicit/gbn2.xml');mod=app.Modeller(pdb.topology,pdb.positions)
print('sample',s,'adding hydrogens',flush=True)
mod.addHydrogens(ff,pH=7.0,platform=mm.Platform.getPlatformByName('CPU'))
print('sample',s,'hydrogens ready',flush=True)
positions=mod.positions.value_in_unit(unit.nanometer);atoms=list(mod.topology.atoms());coords=np.array(positions)
system=ff.createSystem(mod.topology,nonbondedMethod=app.CutoffNonPeriodic,nonbondedCutoff=2*unit.nanometer,constraints=None,removeCMMotion=False)
adidx=[a.index for a in atoms if a.residue.chain.id in 'FGH' and a.element!=app.element.hydrogen]
ax=coords[adidx];flexres=set()
for res in mod.topology.residues():
 if res.chain.id in 'FGH':continue
 at=[a for a in res.atoms() if a.element!=app.element.hydrogen];xyz=coords[[a.index for a in at]]
 if np.linalg.norm(xyz[:,None]-ax[None,:],axis=-1).min()<=.8:flexres.add(res.index)
restraint=mm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)');restraint.addGlobalParameter('k',1000);[restraint.addPerParticleParameter(k) for k in ['x0','y0','z0']]
movable=[];fixedindices=[]
backbone={'N','CA','C','O','OXT','H','H1','H2','H3','HA','HA2','HA3'}
for a in atoms:
 r=a.residue;isadapter=r.chain.id in 'FGH';pos=int(r.id)
 canmove=(isadapter and (pos>=28 or a.name not in backbone)) or (not isadapter and (pos<=48 or (r.index in flexres and a.name not in backbone)))
 if canmove:movable.append(a.index)
 else:system.setParticleMass(a.index,0);fixedindices.append(a.index)
 if isadapter and a.name=='CA' and 28<=pos<=52:restraint.addParticle(a.index,coords[a.index].tolist())
system.addForce(restraint);integrator=mm.VerletIntegrator(.001*unit.picoseconds);ctx=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'4'});ctx.setPositions(mod.positions)
e0=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);print('sample',s,'initial energy',e0,'movable',len(movable),flush=True)
start=time.time();mm.LocalEnergyMinimizer.minimize(ctx,10,500)
state=ctx.getState(getEnergy=True,getPositions=True,getForces=True);e1=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);final=np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer));forces=state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
assert np.isfinite(e1) and np.isfinite(final).all();assert np.max(np.linalg.norm(final[fixedindices]-coords[fixedindices],axis=1))<1e-7
out=ROOT/f'template55_source{s}_FLEX_NTERM_MINIMIZED.pdb'
with out.open('w') as f:
 f.write('REMARK 900 REFERENCE-GUIDED 55AA; BOUNDED LOCAL MINIMIZATION ONLY\nREMARK 900 NOT AN INDEPENDENT PREDICTION OR VALIDATED BINDING MODEL\n')
 app.PDBFile.writeFile(mod.topology,state.getPositions(),f,keepIds=True)
report={'source':str(source),'sample':s,'openmm_version':mm.__version__,'status':'bounded_local_minimization_completed_not_binding_validated','forcefield':['amber14-all.xml','implicit/gbn2.xml'],'nonbonded_cutoff_nm':2,'hydrogen_pH':7,'added_terminal_oxygens':added,'seed':5500+s,'fixed':'penton residues49-543 backbone and outside-shell sidechains; adapter backbone1-27; penton1-48 fully mobile','movable_atoms':len(movable),'adapter_CA28_52_position_restraint_k_kJ_mol_nm2':1000,'distance_contact_restraints':False,'max_iterations':500,'tolerance_kJ_mol_nm':10,'initial_energy_kJ_mol':e0,'final_energy_kJ_mol':e1,'movable_force_rms_kJ_mol_nm':float(np.sqrt(np.mean(np.sum(forces[movable]**2,axis=1)))),'convergence_note':'iteration limit bounded; do not assume force convergence merely from return','elapsed_seconds':time.time()-start,'output':str(out)}
(ROOT/f'relaxation_flex_{s}.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
