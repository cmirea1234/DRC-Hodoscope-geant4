"""Validate a closed batch result; runtime closure is not physical validation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from uproot.source.file import MemmapSource

parser=argparse.ArgumentParser();parser.add_argument('job',type=Path);args=parser.parse_args()
job=args.job.resolve();state=json.loads((job/'STATUS.json').read_text())
assert state['status']=='PROCESS_COMPLETE_ROOT_VALIDATION_REQUIRED','Do not read active or failed ROOT'
try:
    with uproot.open(job/'result.root',handler=MemmapSource) as root:
        e=root['events'].arrays(library='np');sel=(e['t1']>0)&(e['t2']>0);n=int(sel.sum())
        assert n==state['target'] and np.array_equal(e['event'],np.arange(len(sel)))
        assert np.array_equal(e['accepted_t1t2_cumulative'],np.cumsum(sel))
        assert np.all(e['primary_z_mm']==150)
        assert not np.any(e['edep_concrete_mev']) and not np.any(e['edep_blackout_mev'])
        for key in ['measured_noise_enabled','waveform_digitization_enabled','scintillation_yield_scale','pde_scale','coupling_efficiency']:
            assert np.all(e[key]==1),key
        cfg=root['v21_configuration'].arrays(library='np')
        assert cfg['h2_square_fibre_count'][0]==29
        ch=root['channels'].arrays(library='np')
        assert len(ch['event'])==128*n
        assert np.all(ch['event'].reshape(-1,128)==e['event'][sel,None])
        assert np.all(ch['global_channel'].reshape(-1,128)==np.arange(128)[None,:])
        hit=(ch['final_adc']>ch['threshold_adc']).reshape(-1,4,32).any(axis=2)
        for i,k in enumerate(['round_x','round_y','square_x','square_y']):
            assert np.array_equal(hit[:,i],e[f'park_{k}_hit'][sel]>0),k
        assert np.array_equal(hit[:,0]&hit[:,1],e['park_round_xy_hit'][sel]>0)
        assert np.array_equal(hit[:,2]&hit[:,3],e['park_square_xy_hit'][sel]>0)
        loss=root['loss_summary'].arrays(library='np');assert len(loss['event'])==4*n
        for key in ['unclassified','duplicate_terminal_marks','terminal_minus_registered','registered_minus_legacy_produced']:
            assert np.all(loss[key]==0),key
        assert root['photons'].num_entries==0 and root['photon_audit'].num_entries==0
        counts={name:int(h.sum()) for name,h in zip(['MID1','MID9','MID15','MID13'],hit.T)}
        h1=hit[:,0]&hit[:,1];h2=hit[:,2]&hit[:,3]
        counts.update(H1=int(h1.sum()),H2=int(h2.sum()),all_four=int((h1&h2).sum()))
        result=dict(status='PASS',generated=len(sel),T1T2=n,counts=counts,
                    produced=int(e['produced'][sel].sum()),detected=int(e['detected'][sel].sum()),
                    root_sha256=hashlib.sha256((job/'result.root').read_bytes()).hexdigest())
    (job/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
except BaseException as exc:
    (job/'VALIDATION.json').write_text(json.dumps(dict(status='FAIL',error=repr(exc)),indent=2)+'\n')
    raise
