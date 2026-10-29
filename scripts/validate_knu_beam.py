"""Validate a closed fixed-generated-count electron smoke run."""
import json,sys,hashlib,re
from pathlib import Path
import numpy as np,uproot
p=Path(sys.argv[1]);inp=json.loads((p/'INPUT.json').read_text());n=inp['beam']['generated_events']
assert json.loads((p/'resources.json').read_text())['returncode']==0
log=(p/'run.log').read_text()
assert 'KNU_BEAM_PROCESS_COMPLETE' in log
assert not re.search(r'G4Exception-START|Stuck Track|GeomNav[0-9]|GeomVol[0-9]|Segmentation fault|COMMAND NOT FOUND',log)
with uproot.open(p/'result.root') as f:
 e=f['events'].arrays(library='np');c=f['channels'].arrays(library='np');l=f['loss_summary'].arrays(library='np')
 assert len(e['event'])==n and np.array_equal(e['event'],np.arange(n))
 for axis,value in zip('xyz',inp['beam']['position_mm']):assert np.allclose(e['primary_'+axis+'_mm'],value,rtol=0,atol=1e-10)
 for axis,value in zip('xyz',inp['beam']['direction']):assert np.allclose(e['primary_dir_'+axis],value,rtol=0,atol=1e-12)
 assert np.allclose(e['primary_energy_gev'],inp['beam']['kinetic_energy_GeV'],rtol=0,atol=1e-10)
 sel=(e['t1']>0)&(e['t2']>0)
 for k in ('t1','t2'):assert np.array_equal(e[k]>0,e[k+'_edep_mev']>=.05)
 assert np.array_equal(e['accepted_t1t2_cumulative'],np.cumsum(sel))
 assert len(c['event'])==n*128 and np.array_equal(c['event'].reshape(n,128),np.repeat(np.arange(n)[:,None],128,axis=1))
 assert np.array_equal(c['global_channel'].reshape(n,128),np.tile(np.arange(128),(n,1)))
 hit=(c['final_adc']>c['threshold_adc']).reshape(n,4,32).any(axis=2)
 for i,key in enumerate(('round_x','round_y','square_x','square_y')):assert np.array_equal(hit[:,i],e['park_'+key+'_hit']>0)
 assert np.array_equal(hit[:,0]&hit[:,1],e['park_round_xy_hit']>0)
 assert np.array_equal(hit[:,2]&hit[:,3],e['park_square_xy_hit']>0)
 assert len(l['event'])==4*n
 for k in ('unclassified','duplicate_terminal_marks','terminal_minus_registered','registered_minus_legacy_produced'):assert not np.any(l[k])
 for k in ('edep_concrete_mev','edep_blackout_mev'):assert not np.any(e[k])
 for k in ('measured_noise_enabled','waveform_digitization_enabled','scintillation_yield_scale','pde_scale','coupling_efficiency'):assert np.all(e[k]==1)
 assert f['photons'].num_entries==f['photon_audit'].num_entries==0
 r={'status':'PASS','generated_electrons':n,'T1T2':int(sel.sum()),'all_four':int(hit.all(axis=1).sum()),'T1T2_and_all_four':int((sel&hit.all(axis=1)).sum()),'produced':int(e['produced'].sum()),'detected':int(e['detected'].sum()),'root_sha256':hashlib.sha256((p/'result.root').read_bytes()).hexdigest(),'scope':'fixed pencil-beam smoke; not an efficiency estimate or final KEK geometry'}
(p/'VALIDATION.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
