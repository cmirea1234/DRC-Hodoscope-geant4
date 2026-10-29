"""Compose original response + named layout + beam input, preserving baseline payload."""
import argparse,json,math,shutil,hashlib
from pathlib import Path
from verify_bundle import verify,BASE
from layout_config import validate_layout,pose_text
p=argparse.ArgumentParser();p.add_argument('--layout',type=Path,required=True);p.add_argument('--beam',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
verify();layout=json.loads(a.layout.read_text());beam=json.loads(a.beam.read_text())
assert beam['particle']=='e-' and 0<beam['kinetic_energy_GeV']<=6.5
assert isinstance(beam['generated_events'],int) and 0<beam['generated_events']<=100
assert isinstance(beam['seed'],int) and 0<beam['seed']<2000000000
for k in ('position_mm','direction'):
 assert len(beam[k])==3 and all(isinstance(x,(int,float)) and math.isfinite(x) for x in beam[k])
module_layout=layout.get('schema_version')==2
if module_layout:validate_layout(layout)
half=layout['world_half_size_mm'] if module_layout else [250,250,250]
assert all(abs(x)<h for x,h in zip(beam['position_mm'],half)), 'Beam source is outside world'
assert abs(sum(x*x for x in beam['direction'])-1)<1e-10
units={k:('deg' if k.endswith('Yaw') else 'mm') for k in ('trigger1Z','trigger2Z','h1FanoutZ','h2FanoutZ','trigger1OffsetX','trigger1OffsetY','trigger2OffsetX','trigger2OffsetY','h1OffsetX','h1OffsetY','h2OffsetX','h2OffsetY','trigger1Yaw','trigger2Yaw','h1xYaw','h1yYaw','h2xYaw','h2yYaw')}
if not module_layout:assert layout['parameters'].keys()==units.keys()
macro=(BASE/'macros/original.mac').read_text()+'\n# Explicit layout overrides\n'
for k,d in ({} if module_layout else layout['parameters']).items():
 assert d['unit']==units[k] and isinstance(d['value'],(int,float)) and math.isfinite(d['value'])
 assert abs(d['value'])<=150, 'Large layout changes require world/containment validation first.'
 macro+=f"/hodo/{k} {d['value']} {d['unit']}\n"
seed=beam['seed'];n=beam['generated_events']
macro+=f'''# Beam run uses a fixed number of incident electrons, not a trigger target.
/hodo/cosmicBeam false
/hodo/useBuilding false
/hodo/fixedPhotonTest false
/hodo/skipOpticalUnlessT1T2 false
/hodo/recordSelectedDiagnosticsOnly false
/hodo/targetCoincidences 0
/hodo/cosmicSeed {seed}
/random/setSeeds {seed} {seed+17}
/hodo/diagnosticRootFile result
/run/initialize
# Preserve arbitrary beam coordinates through Geant4 unit conversion.
/control/useDoublePrecision true
/gun/particle {beam['particle']}
/gun/energy {beam['kinetic_energy_GeV']} GeV
/gun/position {' '.join(map(str,beam['position_mm']))} mm
/gun/direction {' '.join(map(str,beam['direction']))}
/process/inactivate Cerenkov
/geometry/navigator/check_mode true
/geometry/navigator/push_notify true
/run/beamOn {n}
/control/echo KNU_BEAM_PROCESS_COMPLETE
'''
out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
for f in (BASE/'data/runtime').iterdir():
 if f.is_file():shutil.copy2(f,out/f.name)
(out/'run.mac').write_text(macro)
if module_layout:(out/'layout.pose').write_text(pose_text(layout))
record={'model':'original','required_executable':'hodoscope-original-layout' if module_layout else 'hodoscope-original','status':'PREPARED','layout':layout,'beam':beam,'macro_sha256':hashlib.sha256(macro.encode()).hexdigest(),'layout_sha256':hashlib.sha256(a.layout.read_bytes()).hexdigest(),'beam_sha256':hashlib.sha256(a.beam.read_bytes()).hexdigest(),'full_optical_for_all_incident_events':True}
(out/'INPUT.json').write_text(json.dumps(record,indent=2)+'\n');print('BEAM_PREPARED',out)
