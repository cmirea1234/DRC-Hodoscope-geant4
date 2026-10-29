"""Diagnostic optical probes, separate from cosmic physics comparison."""
import hashlib,json,os,re,shutil,subprocess,sys,time
from pathlib import Path
base=Path(__file__).resolve().parents[1]
out=Path(sys.argv[1]);out.mkdir(exist_ok=False)
assert os.environ.get('_CONDOR_JOB_AD')
py=str(base/'.venv/bin/python')
build=Path(sys.argv[2])
cases=[('original_active','original','active',None),('light32_active','light32','active',None),('light64_control_active','light32','active',64),]
for name,model,source,sides in cases:
    folder=out/name;folder.mkdir()
    for p in (base/'data/runtime').iterdir():
        if p.is_file():shutil.copy2(p,folder/p.name)
    macro=(base/'macros'/f'{model}.mac').read_text()+f'''
# Isolated diagnostic: one optical photon per event, cycle over all H1 X channels.
/hodo/cosmicBeam false
/hodo/fixedPhotonTest true
/hodo/fixedPhotonPlane h1x
/hodo/fixedPhotonChannel -1
/hodo/fixedPhotonSource {source}
/hodo/fixedPhotonAngularMode forward
/hodo/fixedPhotonSpectrum mono
/hodo/fixedPhotonWavelength 450 nm
/hodo/skipOpticalUnlessT1T2 false
/hodo/recordSelectedDiagnosticsOnly false
/hodo/gate4OpticalAudit true
/hodo/recordPhotonAuditRows true
/hodo/targetCoincidences 0
/random/setSeeds 261029101 261029118
/hodo/cosmicSeed 261029101
/hodo/diagnosticRootFile result
'''
    if sides:macro+=f'/hodo/roundTessellatedSides {sides}\n'
    macro+='''/run/initialize
/process/inactivate Cerenkov
/geometry/navigator/check_mode true
/geometry/navigator/push_notify true
/run/beamOn 3200
/control/echo KNU_PHOTON_PROBE_COMPLETE
'''
    (folder/'run.mac').write_text(macro)
    print('PROBE_START',name,flush=True)
    with (folder/'run.log').open('w') as log:
        subprocess.run([py,str(base/'scripts/measure_command.py'),'--output',str(folder/'resources.json'),'--',str(build/f'hodoscope-{model}'),'run.mac'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=1800,check=True)
    text=(folder/'run.log').read_text(errors='replace')
    assert 'KNU_PHOTON_PROBE_COMPLETE' in text
    assert not re.search(r'G4Exception-START|Stuck Track|GeomNav\d|GeomVol\d|Segmentation fault',text),name
    import uproot,numpy as np
    with uproot.open(folder/'result.root') as f:
        a=f['loss_summary'].arrays(library='np'); mask=a['plane']==0
        totals={k:int(v[mask].sum()) for k,v in a.items() if k not in ('event','plane')}
        assert totals['registered']==3200 and totals['terminal_total']==3200,totals
        assert totals['unclassified']==totals['duplicate_terminal_marks']==totals['terminal_minus_registered']==0
        result={'case':name,'diagnostic_only':True,'photons':3200,'emission_disk_radius_mm':0.47,'model':model,'source':source,'sides_override':sides,'totals':totals,'macro_sha256':hashlib.sha256(macro.encode()).hexdigest()}
        (folder/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
        print('PROBE_PASS',name,json.dumps(totals),flush=True)
print('ALL_PROBES_PASS',flush=True)
