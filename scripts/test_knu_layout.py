"""Compute-node tests: baseline event regression and transformed placements/navigation."""
import sys,json,subprocess,os,shutil
from pathlib import Path
import numpy as np,uproot
from layout_config import pose_text,rotation_matrix,ORDER
base=Path(__file__).resolve().parents[1];attempt=Path(sys.argv[1]);build=Path(sys.argv[2]);py=str(base/'.venv/bin/python')
layout=json.loads((base/'configs/layouts/modules_baseline.json').read_text())
pose=attempt/'baseline.pose';pose.write_text(pose_text(layout));env=dict(os.environ,KNU_LAYOUT_FILE=str(pose))
subprocess.run([py,str(base/'scripts/run.py'),'original','--target','5','--seed','261026304','--build-dir',str(build),'--output',str(attempt/'baseline'),'--timeout-hours','1'],env=env,check=True)
subprocess.run([py,str(base/'scripts/validate_root.py'),str(attempt/'baseline')],check=True)
old=base/'runs/knu_20261028_1512/job_2188629_0_bQsF6C/original_5'
with uproot.open(old/'result.root') as a,uproot.open(attempt/'baseline/result.root') as b:
    differences={}
    for tree in ('events','channels','loss_summary'):
        x=a[tree].arrays(library='np');y=b[tree].arrays(library='np')
        bad=[k for k in x if not np.array_equal(x[k],y[k],equal_nan=True)]
        differences[tree]={k:{'max_abs':float(np.nanmax(np.abs(x[k]-y[k])))} for k in bad}
        for k in x:
            if x[k].dtype.kind in 'iu':assert np.array_equal(x[k],y[k]),(tree,k)
            else:assert np.allclose(x[k],y[k],rtol=0,atol=1e-10,equal_nan=True),(tree,k)

(attempt/'baseline_regression.json').write_text(json.dumps(differences)+'\n')
print('LAYOUT_BASELINE_NUMERIC_PASS',flush=True)
R=np.array(rotation_matrix([7,11,13]));shift=np.array([20.,-10.,30.])
rigid=json.loads(json.dumps(layout));rigid['name']='test_global_rigid_transform';rigid['world_half_size_mm']=[400,400,400]
for p in rigid['modules'].values():p['position_mm']=(R@np.array(p['position_mm'])+shift).tolist();p['rotation_deg']=[7,11,13]
independent=json.loads(json.dumps(layout));independent['name']='test_independent_modules';independent['world_half_size_mm']=[400,400,600]
independent['modules']['T1']={'position_mm':[-25,0,200],'rotation_deg':[5,3,10]}
independent['modules']['T2']={'position_mm':[0,-25,140],'rotation_deg':[-4,2,-8]}
independent['modules']['H1']={'position_mm':[3,-2,0],'rotation_deg':[10,5,12]}
independent['modules']['H2']={'position_mm':[-4,1,-150],'rotation_deg':[-8,7,-9]}
for name,config in [('rigid',rigid),('independent',independent)]:
    folder=attempt/name;folder.mkdir();p=folder/'layout.pose';p.write_text(pose_text(config));(folder/'layout.json').write_text(json.dumps(config,indent=2)+'\n')
    macro=(base/'macros/original.mac').read_text()+'''
/hodo/cosmicBeam false
/hodo/fixedPhotonTest true
/hodo/fixedPhotonPlane h1x
/hodo/fixedPhotonChannel -1
/hodo/fixedPhotonSource terminal
/hodo/skipOpticalUnlessT1T2 false
/hodo/recordSelectedDiagnosticsOnly false
/hodo/gate4OpticalAudit true
/hodo/recordPhotonAuditRows true
/hodo/targetCoincidences 0
/hodo/cosmicSeed 261029333
/hodo/diagnosticRootFile result
/run/initialize
/process/inactivate Cerenkov
/run/beamOn 32
/control/echo LAYOUT_GEOMETRY_COMPLETE
'''
    (folder/'run.mac').write_text(macro)
    with (folder/'run.log').open('w') as log:subprocess.run([str(build/'hodoscope-original-layout'),'run.mac'],cwd=folder,env=dict(os.environ,KNU_LAYOUT_FILE=str(p)),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    assert 'MODULE_POSE_NAVIGATION_PASS' in (folder/'run.log').read_text()
    before=[json.loads(x) for x in (attempt/'baseline/layout_placements.jsonl').read_text().splitlines()]
    after=[json.loads(x) for x in (folder/'layout_placements.jsonl').read_text().splitlines()]
    assert len(before)==len(after)==554
    for a,b in zip(before,after):
        assert (a['name'],a['copy'],a['module'])==(b['name'],b['copy'],b['module'])
        module=ORDER[a['module']];oldp=np.array(layout['modules'][module]['position_mm']);newp=np.array(config['modules'][module]['position_mm']);r=np.array(rotation_matrix(config['modules'][module]['rotation_deg']))
        assert np.allclose(b['position_mm'],r@(np.array(a['position_mm'])-oldp)+newp,rtol=0,atol=1e-9)
        assert np.allclose(np.array(b['rotation']).reshape(3,3),r@np.array(a['rotation']).reshape(3,3),rtol=0,atol=1e-12)
    with uproot.open(folder/'result.root') as root:
        loss=root['loss_summary'].arrays(library='np')
        assert int(loss['registered'].sum())==32
        assert int(loss['terminal_total'].sum())==32
        assert int(loss['stage_sipm_boundary'].sum())>=28
        assert not np.any(loss['unclassified'])
    print('LAYOUT_TRANSFORM_PASS',name,'placements=554','optical_readout_pass',flush=True)
(attempt/'VALIDATION.json').write_text(json.dumps({'status':'PASS','baseline_integer_exact':True,'baseline_float_absolute_tolerance':1e-10,'baseline_bitwise_exact':not any(differences.values()),'rigid':True,'independent':True})+'\n')
print('LAYOUT_ALL_PASS',flush=True)
