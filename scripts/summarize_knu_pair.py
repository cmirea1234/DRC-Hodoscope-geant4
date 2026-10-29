"""Compare only closed, validated original/light32 paired results; never tune inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from uproot.source.file import MemmapSource

parser = argparse.ArgumentParser()
parser.add_argument('attempt', type=Path)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--target', type=int, required=True)
parser.add_argument('--seed', type=int, required=True)
args = parser.parse_args()
base = Path(__file__).resolve().parents[1]
report = {'scope': f'{args.target}-trigger trial per model; one seed pair is not a general performance conclusion',
          'denominator': 'T1 AND T2 energy deposit >= 0.05 MeV each; not strict primary path',
          'attempt': str(args.attempt.resolve()), 'models': {}}
primary = {}
selected = {}
hit_masks = {}
for model in ('original', 'light32'):
    folder = args.attempt / (model + '_' + str(args.target))
    state = json.loads((folder / 'STATUS.json').read_text())
    validation = json.loads((folder / 'VALIDATION.json').read_text())
    assert state['status'] == 'PROCESS_COMPLETE_ROOT_VALIDATION_REQUIRED'
    assert validation['status'] == 'PASS'
    assert state['target'] == args.target and state['seed'] == args.seed
    root_hash = hashlib.sha256((folder / 'result.root').read_bytes()).hexdigest()
    assert root_hash == validation['root_sha256']
    expected_macro = (base / 'macros' / (model + '.mac')).read_text()
    macro = (folder / 'run.mac').read_text()
    assert macro.startswith(expected_macro)
    assert hashlib.sha256(macro.encode()).hexdigest() == state['macro_sha256']
    with uproot.open(folder / 'result.root', handler=MemmapSource) as root:
        events = root['events'].arrays(library='np')
        channels = root['channels'].arrays(library='np')
        cfg = root['v21_configuration'].arrays(library='np')
        sel = (events['t1'] > 0) & (events['t2'] > 0)
        assert np.array_equal(events['t1'] > 0, events['t1_edep_mev'] >= 0.05)
        assert np.array_equal(events['t2'] > 0, events['t2_edep_mev'] >= 0.05)
        assert np.all((events['primary_energy_gev'] >= 0.2) & (events['primary_energy_gev'] <= 1000))
        assert np.all(events['primary_zenith_deg'] <= 45 + 1e-9)
        assert cfg['h1_offset_x_mm'][0] == cfg['h1_offset_y_mm'][0] == -0.25
        assert np.array_equal(channels['plane'].reshape(-1,128), np.tile(np.repeat(np.arange(4),32),(int(sel.sum()),1)))
        assert np.array_equal(channels['local_channel'].reshape(-1,128), np.tile(np.arange(32), (int(sel.sum()),4)))
        names = ('MID1', 'MID9', 'MID15', 'MID13')
        planes = ('round_x', 'round_y', 'square_x', 'square_y')
        masks = {name: events['park_' + plane + '_hit'][sel] > 0 for name,plane in zip(names,planes)}
        masks.update(H1=masks['MID1'] & masks['MID9'], H2=masks['MID15'] & masks['MID13'])
        masks['all_four'] = masks['H1'] & masks['H2']
        hit_masks[model] = masks
        n = int(sel.sum())
        rates = {k: {'count': int(v.sum()), 'denominator': n, 'fraction': float(v.mean()),
                     'binomial_1sigma': float(np.sqrt(v.mean()*(1-v.mean())/n))} for k,v in masks.items()}
        paths = {}
        for category, suffix in [('core','_primary_core_track_mm'), ('active_core','_primary_active_core_track_mm')]:
            m = {name: events[name.lower()+suffix][sel] > 0 for name in names}
            m.update(H1=m['MID1'] & m['MID9'], H2=m['MID15'] & m['MID13'])
            m['all_four'] = m['H1'] & m['H2']
            paths[category] = {k:int(v.sum()) for k,v in m.items()}
        primary[model] = {k: events[k] for k in events if k.startswith('primary_') and k in (
            'primary_x_mm','primary_y_mm','primary_z_mm','primary_dir_x','primary_dir_y','primary_dir_z','primary_energy_gev')}
        selected[model] = events['event'][sel]
        optical = {k: int(events[k][sel].sum()) for k in (
            'produced','detected','h1_sipm_boundary_photons','h2_sipm_boundary_photons',
            'h1_pde_rejected_photons','h2_pde_rejected_photons')}
        peak = {plane: events[plane+'_adc_max'][sel].tolist() for plane in planes}
        report['models'][model] = dict(validation=validation,
            hit_rates=rates, primary_path_counts=paths, optical_totals=optical,
            selected_event_ids=selected[model].tolist(), peak_adc_values=peak,
            configuration={k: v.tolist() for k,v in cfg.items()},
            process_elapsed_seconds=state['elapsed_seconds'], seed=state['seed'],
            build=state['build'], resources=json.loads((args.attempt/(model+'.resources.txt')).read_text()),
            binary_sha256=state['binary_sha256'])
common = min(len(primary[m]['primary_x_mm']) for m in primary)
report['primary_comparison'] = {'common_prefix_length': common,
    'exactly_equal_fields': {k: bool(np.array_equal(primary['original'][k][:common],primary['light32'][k][:common])) for k in primary['original']},
    'selected_event_ids_equal': bool(np.array_equal(selected['original'],selected['light32']))}
matched = all(report['primary_comparison']['exactly_equal_fields'].values()) and report['primary_comparison']['selected_event_ids_equal']
report['paired_hit_difference'] = {}
if matched:
    for key in hit_masks['original']:
        diff = hit_masks['light32'][key].astype(float)-hit_masks['original'][key].astype(float)
        report['paired_hit_difference'][key] = {'light32_minus_original':float(diff.mean()),
            'paired_standard_error':float(diff.std(ddof=1)/np.sqrt(len(diff))) if len(diff)>1 else None}
report['notes'] = ['A single seed pair cannot establish model equivalence or a general speed advantage; cache and execution order are uncontrolled.',
                   'Wald binomial errors and paired errors can be zero for all-equal tiny samples; zero does not establish certainty.',
                   'Primary equality does not imply identical optical histories.',
                   'Both runs use the same node, build configuration, recording settings and Geant4 libraries.']
args.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({m:report['models'][m]['hit_rates'] for m in report['models']},indent=2))
print('PRIMARY_COMPARISON',report['primary_comparison'])
