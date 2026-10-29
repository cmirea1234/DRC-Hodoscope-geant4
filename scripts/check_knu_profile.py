"""Require timer-only instrumentation to preserve all numerical event outputs."""
import sys,json
from pathlib import Path
import numpy as np
import uproot
p=Path(sys.argv[1]); report={}
for model in ('original','light32'):
    plain=p/('baseline_'+model);timed=p/('instrumented_'+model)
    for folder in (plain,timed):
        assert json.loads((folder/'VALIDATION.json').read_text())['status']=='PASS'
    checks={}
    with uproot.open(plain/'result.root') as a,uproot.open(timed/'result.root') as b:
        for tree in ('events','channels','loss_summary','v21_configuration'):
            x=a[tree].arrays(library='np');y=b[tree].arrays(library='np')
            assert x.keys()==y.keys()
            bad=[k for k in x if not np.array_equal(x[k],y[k],equal_nan=True)]
            checks[tree]={'branches':len(x),'different':bad}
            assert not bad,(model,tree,bad)
    timings=json.loads((timed/'profile.json').read_text())
    total=sum(v['exclusive_wall_seconds'] for v in timings.values())
    baseline=json.loads((p/f'baseline_{model}.resources.json').read_text())
    instrumented=json.loads((p/f'instrumented_{model}.resources.json').read_text())
    report[model]={'checks':checks,'timings':timings,'sum_profile_wall_seconds':total,'baseline_resources':baseline,'instrumented_resources':instrumented}
    print(model,'NUMERIC_REGRESSION_PASS','wall baseline/instrumented',baseline['wall_seconds'],instrumented['wall_seconds'])
    for k,v in sorted(timings.items(),key=lambda kv:-kv[1]['exclusive_wall_seconds']):
        print(k,round(v['exclusive_wall_seconds'],3),round(v['exclusive_wall_seconds']/total*100,2),'%',v['calls'])
(p/'profile_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
