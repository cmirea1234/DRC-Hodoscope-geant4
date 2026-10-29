"""Compare traced core-side encounters and independently check cylinder normals."""
from pathlib import Path
import math,json
import numpy as np,uproot
base=Path(__file__).resolve().parents[1];p=base/'runs/knu_diagnostic_20261029/boundary_2188883_0_2VJ8PE'
mask=(1<<64)-1
def uniform(event,stream):
 v=(261029101+0x9e3779b97f4a7c15*(event+1)+0xbf58476d1ce4e5b9*(stream+101))&mask
 v=((v^(v>>30))*0xbf58476d1ce4e5b9)&mask;v=((v^(v>>27))*0x94d049bb133111eb)&mask;v^=v>>31
 return (v>>11)/9007199254740992.
report={}
for case in sorted((p/'cases').iterdir()):
 records={}
 for line in (case/'run.log').read_text().splitlines():
  if not line.startswith('KNU_BOUNDARY '):continue
  fields=dict(s.split('=',1) for s in line.split()[1:]);event=int(fields['event'])
  if fields['pre_region'] not in ('1','2','12') or fields['candidate_region'] not in ('3','4','5','6','13','14') or float(fields['length_mm'])<=1e-9:continue
  v={k:float(x) for k,x in fields.items() if k not in ('pre','post')}
  dot=sum(v['d'+axis]*v['n'+axis] for axis in 'xyz')
  v['incidence_deg']=math.degrees(math.acos(min(1.,abs(dot))))
  v['critical_deg']=math.degrees(math.asin(v['n2']/v['n1']))
  records.setdefault(event,[]).append(v)
 with uproot.open(case/'result.root') as f:a=f['photon_audit'].arrays(library='np')
 if case.name=='original_active':
  checks=[]
  for e,rows in records.items():
   radius=.47*math.sqrt(uniform(e,3));phi=2*math.pi*uniform(e,4)
   y0=a['origin_y_mm'][e]-radius*math.cos(phi);z0=a['origin_z_mm'][e]-radius*math.sin(phi)
   # Straight active region only: independently known 0.4795 mm circular core.
   for row in rows:
    if row['pre_region']!=1 or row['x']>15.75:continue
    radial=np.array([0.,row['y']-y0,row['z']-z0]);r=np.linalg.norm(radial)
    normal=np.array([row['nx'],row['ny'],row['nz']])
    residual=float(np.linalg.norm(radial/r-normal))
    checks.append({'event':e,'step':int(row['step']),'radius_residual_mm':abs(r-.4795),'normal_residual':residual})
  assert max(v['radius_residual_mm'] for v in checks)<1e-9
  assert max(v['normal_residual'] for v in checks)<1e-8
  report['independent_cylinder_check']={'count':len(checks),'max_radius_residual_mm':max(v['radius_residual_mm'] for v in checks),'max_normal_residual':max(v['normal_residual'] for v in checks)}
 report[case.name]=records
 print(case.name)
 for e in (16,25,17,20):
  rows=records.get(e,[])
  first_non_tir=next((r for r in rows if not r['tir']),None)
  first_transmit=next((r for r in rows if not r['reflected']),None)
  print(e,'first_angles',[(int(r['step']),round(r['incidence_deg'],4),int(r['tir'])) for r in rows[:3]],'critical',round(rows[0]['critical_deg'],4),'first_non_tir_step',None if first_non_tir is None else int(first_non_tir['step']),'transmit_x',None if first_transmit is None else round(first_transmit['x'],6))
(p/'boundary_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(report['independent_cylinder_check'])
