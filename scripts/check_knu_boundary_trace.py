import sys,json
from pathlib import Path
import numpy as np
import uproot
base=Path(__file__).resolve().parents[1];current=Path(sys.argv[1]);old=base/'runs/knu_diagnostic_20261029/disk_2188882_0_vbHOkv'
result={}
for folder in sorted((current/'cases').iterdir()):
 with uproot.open(folder/'result.root') as a,uproot.open(old/'cases'/folder.name/'result.root') as b:
  checks={}
  for tree in ('photon_audit','loss_summary'):
   x=a[tree].arrays(library='np');y=b[tree].arrays(library='np')
   bad=[k for k in x if not np.array_equal(x[k],y[k],equal_nan=True)]
   assert not bad,(folder.name,tree,bad)
   checks[tree]='ALL_NUMERIC_BRANCHES_IDENTICAL'
 result[folder.name]=checks
(current/'observer_regression.json').write_text(json.dumps(result,indent=2)+'\n')
print('BOUNDARY_OBSERVER_REGRESSION_PASS')
