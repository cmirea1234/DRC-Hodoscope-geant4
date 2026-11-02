"""Read-only preflight for the supported original simulation environment."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from verify_bundle import BASE, verify

p = argparse.ArgumentParser()
p.add_argument('--patch-record',type=Path,required=True)
p.add_argument('--output',type=Path)
a = p.parse_args()
verify()
for command in ('cmake','g++','geant4-config'):
    assert shutil.which(command), 'Missing command: '+command
version=subprocess.check_output(['geant4-config','--version'],text=True).strip()
assert version=='11.4.1', 'Use the tested Geant4 11.4.1, got '+version
manifest=json.loads((BASE/'patches/manifest.json').read_text())
record=json.loads(a.patch_record.read_text())
assert set(record)==set(manifest), 'Unexpected Geant4 patch record'
for name,item in manifest.items():
    for key in ('original_sha256','patched_sha256'):
        assert record[name][key]==item[key], (name,key)
data=subprocess.check_output(['geant4-config','--check-datasets'],text=True)
lines=[line.split() for line in data.splitlines() if line.strip()]
assert len(lines)==12 and all(len(line)>=3 and line[1]=='INSTALLED' and Path(line[2]).is_dir() for line in lines), data
import numpy,uproot
assert tuple(map(int,numpy.__version__.split('.')[:2]))>=(1,24) and int(numpy.__version__.split('.')[0])<3
assert int(uproot.__version__.split('.')[0])==5
report={'status':'ENVIRONMENT_PASS','geant4':version,'python':sys.version.split()[0],
        'numpy':numpy.__version__,'uproot':uproot.__version__,
        'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
        'cmake':subprocess.check_output(['cmake','--version'],text=True).splitlines()[0],
        'datasets':data,'patch_record':str(a.patch_record.resolve()),
        'scope':'Checks environment and patch provenance record, not binary provenance or physical validation'}
if a.output:a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
