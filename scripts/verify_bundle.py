"""Check payload hashes and the two intended effective configurations."""
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]

def verify():
    entries=json.loads((BASE/'PAYLOAD_SHA256.json').read_text())
    for relative,expected in entries.items():
        p=BASE/relative
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==expected, relative
    for model in ['original','light32']:
        text=(BASE/'macros'/f'{model}.mac').read_text()
        assert '/home/' not in text and '/control/execute' not in text
        commands={s.split()[0]:' '.join(s.split()[1:]) for s in text.splitlines() if s.startswith('/')}
        for key,value in {'useBuilding':'false','h2SquareFibreCount':'29',
                          'cosmicPositionModel':'target_backproject','cosmicAngularDistribution':'modified_gaisser_joint',
                          'cosmicMinEnergy':'0.2 GeV','cosmicMaxEnergy':'1000 GeV',
                          'cosmicTargetHalfX':'16 mm','cosmicTargetHalfY':'16 mm',
                          'maxZenith':'45 deg','scintillationYieldScale':'1',
                          'useWaveformDigitization':'true','useMeasuredNoise':'true',
                          'recordPhotonNtuple':'false','recordChannelNtuple':'true'}.items():
            assert commands['/hodo/'+key]==value,(model,key)
        assert commands['/hodo/roundSweepRepresentation']==('native_balanced_union' if model=='original' else 'tessellated')
        if model=='light32':assert commands['/hodo/roundTessellatedSides']=='32'
    for name in ['ActionInitialization.cc','main.cc']:
        assert (BASE/'models/original/src'/name).read_bytes()==(BASE/'models/light32/src'/name).read_bytes()
    print('BUNDLE_INTEGRITY_PASS',len(entries),'files')

if __name__=='__main__':verify()
