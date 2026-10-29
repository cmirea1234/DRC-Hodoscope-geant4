"""Read completed 500-trigger outputs; align by event ID before comparing."""
import json
from pathlib import Path
import numpy as np
import uproot
base=Path(__file__).resolve().parents[1]
source=base/'runs/knu_500_20261028_1643/job_2188630_0_uWeFDu'
out=base/'runs/knu_diagnostic_20261029/summary/loss_comparison.json'
events={};losses={}
for model in ('original','light32'):
    folder=source/(model+'_500')
    assert json.loads((folder/'VALIDATION.json').read_text())['status']=='PASS'
    with uproot.open(folder/'result.root') as root:
        keys=[k for k in root['events'].keys() if k in ('event','t1','t2') or any(x in k for x in ('photons_produced','core_edep','stage_','photons_detected'))]
        events[model]=root['events'].arrays(keys,library='np')
        losses[model]=root['loss_summary'].arrays(library='np')
selected={m: e['event'][(e['t1']>0)&(e['t2']>0)] for m,e in events.items()}
common=np.intersect1d(*selected.values())
report={'common_selected_count':len(common),'models':{},'notes':['Common selected events only; optical histories are not assumed identical.','Stage counts overlap and must not be summed as mutually exclusive losses.','Comparing only joint selection is conditional and does not replace each full sample.']}
for model,e in events.items():
    mask=np.isin(e['event'],common);l=losses[model]
    report['models'][model]={'common_event_totals':{k:float(v[mask].sum()) for k,v in e.items() if k not in ('event','t1','t2')},'planes':{}}
    for plane in range(4):
        sel=np.isin(l['event'],common)&(l['plane']==plane)
        report['models'][model]['planes'][plane]={k:int(v[sel].sum()) for k,v in l.items() if k not in ('event','plane')}
out.write_text(json.dumps(report,indent=2)+'\n')
for m,d in report['models'].items():
    print(m)
    for plane,v in d['planes'].items():
        print(plane,{k:v[k] for k in ('registered','stage_active_exit','stage_core_to_grease','stage_window_entry','stage_sipm_boundary','sipm_detected','bulk_core','bulk_inner_clad','world_escape','boundary_absorption','other_terminal')})
print('COMMON',len(common))
