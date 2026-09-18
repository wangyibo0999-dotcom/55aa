from pathlib import Path
import json,shlex,itertools,hashlib
import numpy as np
from _paths import DATA, RESULTS, SCRIPTS, output_dir
ROOT=output_dir('interface')
BASE=RESULTS/'predictions'
REF=DATA/'reference_1X9T/1X9T.cif'
OLD=(DATA/'reference52_original')
aa=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def parse(path):
 headers=[];out={}
 for l in path.read_text().splitlines():
  if l.startswith('_atom_site.'):headers.append(l.strip().split()[0].split('.',1)[1])
  elif l.startswith('ATOM '):
   d=dict(zip(headers,shlex.split(l)))
   if d['type_symbol'] in ('H','D') or float(d.get('occupancy',1))<=0:continue
   c=d['label_asym_id'];r=int(d['label_seq_id']);v={'pos':r,'auth':int(d['auth_seq_id']),'name':d['label_atom_id'],'aa':aa[d['label_comp_id']],'b':float(d['B_iso_or_equiv']),'xyz':np.array([float(d['Cartn_'+a]) for a in 'xyz'])}
   out.setdefault(c,[]).append(v)
 return out
def ca(atoms):return {a['pos']:a for a in atoms if a['name']=='CA'}
def xyz(atoms):return np.array([a['xyz'] for a in atoms])
def fit(x,y):
 xc=x.mean(0);yc=y.mean(0);u,s,vt=np.linalg.svd((x-xc).T@(y-yc));m=np.eye(3);m[2,2]=np.linalg.det(u@vt);rot=u@m@vt;t=yc-xc@rot
 return rot,t,float(np.sqrt(np.mean(np.sum((x@rot+t-y)**2,axis=1))))
# Exact alignment supplied by the user from EMBOSS Needle, stored as paired blocks.
blocks=[
('--------------------------------------------TGGRNS','MMRRAYPEGPPPSYESVMQQAMAAAAAMQPPLEAPYVPPRYLAPTEGRNS'),
('IRYSELAPLFDTTRVYLVDNKSTDVASLNYQNDHSNFLTTVIQNNDYSPG','IRYSELSPLYDTTRLYLVDNKSADIASLNYQNDHSNFLTTVVQNNDFTPT'),
('EASTQTINLDDRSHWGGDLKTILHTNMPNVNEFMFTNKFKARVMVSRSLT','EASTQTINFDERSRWGGQLKTIMHTNMPNVNEFMYSNKFKARVMVSRKTP'),
('KDKQVE----------LKYEWVEFTLPEGNYSETMTIDLMNNAIVEHYLK','NGEFVTVTDGPGSQDILEYEWVEFELPEGNFSVTMTIDLMNNAIIDNYLA'),
('VGRQNGVLESDIGVKFDTRNFRLGFDPVTGLVMPGVYTNEAFHPDIILLP','VGRQNGVLESDIGVKFDTRNFRLGWDPVTELVMPGVYTNEAFHPDIVLLP'),
('GCGVDFTHSRLSNLLGIRKRQPFQEGFRITYDDLEGGNIPALLDVDAYQA','GCGVDFTESRLSNLLGIRKRQPFQEGFQIMYEDLEGGNIPALLDVDAYEK'),
('SLKDDTEQGGDGAGGGNNSGSGAEENSNAAAAAMQPVEDMNDHAIRGDTF','SKEE----------------SAAAARTAAVATASTEVD------VRGDNF'),
('ATRAEEKRAEAEAAAEAAAPAAQPEVEKPQKKPVIKPLTEDSKKRSYNLI','ASPA----AELVAAAEAA------ETES-SRKIVIQPVEKDSKDRSYNVL'),
('SNDSTFTQYRSWYLAYNYGDPQTGIRSWTLLCTPDVTCGSEQVYWSLPDM','P-DKINTAYRSWYLAYNYGDPEKGVRSWTLLTTSDVTCGVEQVYWSLPDM'),
('MQDPVTFRSTSQISNFPVVGAELLPVHSKSFYNDQAVYSQLIRQFTSLTH','MQDPVTFRSTRQVSNYPVVGAELLPVYSKSFFNEQAVYSQQLRAFTSLTH'),
('VFNRFPENQILARPPAPTITTVSENVPALTDHGTLPLRNSIGGVQRVTIT','VFNRFPENQILVRPPAPTITTVSENVPALTDHGTLPLRSSIRGVQRVTVT'),
('DARRRTCPYVYKALGIVSPRVLSSRTF','DARRRTCPYVYKALGIVAPRVLSSRTF')]
inp=json.loads((BASE/'inputs_exported.json').read_text())[0];seqp=inp['sequences'][0]['proteinChain']['sequence'];seqa=inp['sequences'][1]['proteinChain']['sequence']
al1=''.join(a for a,b in blocks);al2=''.join(b for a,b in blocks)
assert all(len(a)==len(b) for a,b in blocks)
assert al2.replace('-','')==seqp
fastaseq=(REF.parent/'rcsb_pdb_1X9T.fasta').read_text().splitlines()[1];assert al1.replace('-','')==fastaseq
mapping={};i=j=0
for a,b in zip(al1,al2):
 i+=a!='-';j+=b!='-'
 if a!='-' and b!='-':mapping[i]=(j,a==b)
ref=parse(REF);rca=ca(ref['A']);pepca=np.array([next(a['xyz'] for a in ref['B'] if a['auth']==r and a['name']=='CA') for r in range(10,20)])
assert mapping[197-48][0]==203 and mapping[494-48][0]==466 and mapping[227-48][0]==233
# Match every penton permutation, not just chain-letter/copy identity.
oldfile=next(OLD.rglob('*seed_64438_sample_0.cif'));old=parse(oldfile);oca={c:ca(old[c]) for c in 'ABCDE'}
core=[i for i in range(1,544) if min(oca[c][i]['b'] for c in 'ABCDE')>=80]
X=np.concatenate([np.array([oca[c][i]['xyz'] for i in core]) for c in 'ABCDE'])
results=[]
for sample in range(5):
 path=BASE/f'sample_{sample}.cif';chains=parse(path);conf=json.loads((BASE/f'confidence_{sample}.json').read_text());cas={c:ca(v) for c,v in chains.items()}
 assert set(chains)==set('ABCDEFGH')
 for c in chains:assert ''.join(v['aa'] for v in cas[c].values())==(seqp if c in 'ABCDE' else seqa)
 ring=[]
 for perm in itertools.permutations('ABCDE'):
  Y=np.concatenate([np.array([cas[c][i]['xyz'] for i in core]) for c in perm]);R,t,rmsd=fit(X,Y);ring.append((rmsd,''.join(perm),R,t))
 best=min(ring,key=lambda v:v[0]);monos=[];projections={}
 for c in 'ABCDE':
  pairs=[(ri,p) for ri,(p,identical) in mapping.items() if identical and ri in rca and p in cas[c] and cas[c][p]['b']>=70]
  x=np.array([rca[ri]['xyz'] for ri,p in pairs]);y=np.array([cas[c][p]['xyz'] for ri,p in pairs]);R,t,allr=fit(x,y);mask=np.ones(len(x),dtype=bool)
  for _ in range(10):
   new=np.linalg.norm(x@R+t-y,axis=1)<=2.5
   if new.sum()<30:break
   R,t,trimr=fit(x[new],y[new])
   if np.array_equal(new,mask):mask=new;break
   mask=new
  projections[c]=pepca@R+t
  monos.append({'chain':c,'initial_identical_CA_pairs':len(x),'initial_RMSD_A':allr,'core_pairs':int(mask.sum()),'core_RMSD_A':float(np.sqrt(np.mean(np.sum((x[mask]@R+t-y[mask])**2,axis=1))))})
 interfaces=[];adapters=[]
 for ac in 'FGH':
  aa0=chains[ac];xa=xyz(aa0);ar=np.array([a['pos'] for a in aa0]);targetpep=np.array([cas[ac][i]['xyz'] for i in range(43,53)]);placements=[]
  for pc in 'ABCDE':
   pp=chains[pc];xp=xyz(pp);pr=np.array([a['pos'] for a in pp]);d=np.linalg.norm(xa[:,None]-xp[None,:],axis=-1)
   ii,jj=np.where(d<=6);ti,tj=np.where((d<=6)&((ar>=43)&(ar<=52))[:,None]);i4,j4=np.where(d<=4)
   anchors={f'{ap}-{p}':float(d[np.ix_(ar==ap,pr==p)].min()) for ap,p in [(45,203),(48,466),(50,233)]}
   interfaces.append({'adapter':ac,'penton':pc,'min_heavy_A':float(d.min()),'all_residue_pairs_6A':len(set(zip(ar[ii].tolist(),pr[jj].tolist()))),'all_residue_pairs_4A':len(set(zip(ar[i4].tolist(),pr[j4].tolist()))),'motif43_52_pairs_6A':len(set(zip(ar[ti].tolist(),pr[tj].tolist()))),'contact_adapter_positions':sorted(set(ar[ii].tolist())),'contact_penton_positions':sorted(set(pr[jj].tolist())),'atom_pairs_lt2A':int((d<2).sum()),'pair_iptm':conf['chain_pair_iptm'][ord(ac)-65][ord(pc)-65],'anchors':anchors})
   proj=projections[pc];placements.append({'reference_primary_penton':pc,'motif_CA_RMSD_without_refitting_A':float(np.sqrt(np.mean(np.sum((proj-targetpep)**2,axis=1)))),'motif_centroid_distance_A':float(np.linalg.norm(proj.mean(0)-targetpep.mean(0)))})
  ints=[v for v in interfaces if v['adapter']==ac]
  adapters.append({'chain':ac,'motif_CA_plddt_mean':float(np.mean([cas[ac][i]['b'] for i in range(43,53)])),'nearest_reference_placement':min(placements,key=lambda v:v['motif_CA_RMSD_without_refitting_A']),'all_reference_placements':placements,'best_any_chain_anchors':{key:min(({'penton':v['penton'],'distance_A':v['anchors'][key]} for v in ints),key=lambda v:v['distance_A']) for key in ['45-203','48-466','50-233']}})
 row={'sample':sample,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sequence_identity_verified':True,'penton_vs_reference52':{'reference_sample':0,'all_120_permutations_tested':True,'best_map_referenceABCDE_to_prediction':best[1],'matched_CA_count':len(X),'RMSD_A':best[0]},'penton_vs_1X9T':monos,'interfaces':interfaces,'adapters':adapters}
 results.append(row)
 print('sample',sample,'penton ring rmsd',round(best[0],2),'tail rmsd',[round(a['nearest_reference_placement']['motif_CA_RMSD_without_refitting_A'],2) for a in adapters],flush=True)
out={'reference_1X9T':str(REF),'reference_1X9T_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),'reference52':str(oldfile),'reference52_sha256':hashlib.sha256(oldfile.read_bytes()).hexdigest(),'alignment_reference_gapped':al1,'alignment_target_gapped':al2,'samples':results}
(ROOT/'spatial_audit.json').write_text(json.dumps(out,indent=2))
