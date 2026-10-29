"""Compare closed optical probes by actual recorded incident photon state."""
import json,sys
from pathlib import Path
import numpy as np
import uproot
folder=Path(sys.argv[1]); result={'cases':{},'input_comparisons':{}};data={}
for case in sorted((folder/'cases').iterdir()):
    validation=json.loads((case/'validation.json').read_text())
    with uproot.open(case/'result.root') as f:
        a=f['photon_audit'].arrays(library='np');order=np.argsort(a['event']);a={k:v[order] for k,v in a.items()}
        assert len(a['event'])==3200 and np.array_equal(a['event'],np.arange(3200))
        data[case.name]=a
        result['cases'][case.name]={'validation':validation,'resources':json.loads((case/'resources.json').read_text()),'mean_steps':float(a['steps'].mean()),'mean_track_length_mm':float(a['total_track_length_mm'].mean()),'terminal_reason_counts':{str(int(k)):int(v) for k,v in zip(*np.unique(a['terminal_reason'],return_counts=True))},'mean_tir_count':float(a['tir_count'].mean()),'fanout_core_entered':int((a['fanout_core_path_mm']>0).sum())}
keys=['event','origin_global_channel','wavelength_nm','origin_x_mm','origin_y_mm','origin_z_mm','initial_dir_x','initial_dir_y','initial_dir_z']
for source in ('active','terminal'):
    ref='original_'+source
    if ref not in data:continue
    for name,a in data.items():
        if not name.endswith('_'+source) or name==ref:continue
        same={k:bool(np.array_equal(a[k],data[ref][k])) for k in keys}
        result['input_comparisons'][name+'_vs_'+ref]=same
        assert all(same.values()),(name,same)
(folder/'probe_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
for name,v in result['cases'].items():
    print(name,'arrival',v['validation']['totals']['stage_sipm_boundary'],'detected',v['validation']['totals']['sipm_detected'],'steps',round(v['mean_steps'],2),'fanout_entered',v['fanout_core_entered'])
print('EXACT_INCIDENT_PHOTON_INPUT_MATCH_PASS')
