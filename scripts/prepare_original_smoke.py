"""Prepare a first-use KNU Condor build + original five-trigger run in any checkout."""
import argparse
import re
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
from verify_bundle import BASE,verify

p=argparse.ArgumentParser()
p.add_argument('--env-file',type=Path,required=True,help='Shell setup defining G4PREFIX and loading compiler/Geant4')
p.add_argument('--filesystem-domain',required=True,help='Shared filesystem domain advertised by compute nodes')
a=p.parse_args()
assert re.fullmatch(r'[A-Za-z0-9_.-]+',a.filesystem_domain), 'Invalid filesystem domain'
verify()
assert a.env_file.is_file()
assert not any(c.isspace() for c in str(BASE)), 'Use a checkout path without whitespace'
assert (BASE/'.venv/bin/python').is_file(), 'Prepare .venv first (docs/INSTALL_KO.md)'
root=BASE/'runs';root.mkdir(exist_ok=True)
out=Path(tempfile.mkdtemp(prefix='first_original_',dir=root))
shutil.copy2(a.env_file,out/'environment.sh')
q=shlex.quote
worker=f'''#!/usr/bin/env bash
set -eo pipefail
test -n "${{_CONDOR_JOB_AD:-}}"
repo={q(str(BASE))}
campaign={q(str(out))}
cd "$repo"
source "$campaign/environment.sh"
set -u
: "${{G4PREFIX:?Set G4PREFIX in the environment file}}"
unset KNU_LAYOUT_FILE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME CPLUS_INCLUDE_PATH C_INCLUDE_PATH
attempt=$(mktemp -d "$campaign/start_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
py="$repo/.venv/bin/python"
"$py" scripts/check_runtime.py --patch-record "$G4PREFIX/applied_patches.json" --output "$attempt/environment.json"
cmake -S "$repo" -B "$attempt/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$G4PREFIX" -DPython3_EXECUTABLE="$py"
cmake --build "$attempt/build" --target hodoscope-original -j2
ctest --test-dir "$attempt/build" --output-on-failure
"$py" scripts/measure_command.py --output "$attempt/resources.json" -- "$py" scripts/run.py original --target 5 --seed 261026304 --build-dir "$attempt/build" --output "$attempt/original_5" --timeout-hours 2
"$py" scripts/validate_root.py "$attempt/original_5"
"$py" - "$attempt" <<'CHECK'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);r=json.loads((p/'original_5/VALIDATION.json').read_text())
ref=json.loads(Path('reference/original_small.json').read_text())
for k in ('generated','T1T2','produced','detected'):assert r[k]==ref[k],(k,r[k],ref[k])
for k in ('H1','H2','all_four'):assert r['counts'][k]==ref['counts']['channel_threshold_hit'][k]
(p/'FIRST_RUN_PASS.json').write_text(json.dumps({{'status':'PASS','scope':'Five-trigger aggregate regression, not full bitwise/physical validation','validation':r}},indent=2)+'\\n')
CHECK
echo FIRST_ORIGINAL_PASS
'''
(out/'worker.sh').write_text(worker);(out/'worker.sh').chmod(0o755)
submit=f'''universe = vanilla
executable = {out}/worker.sh
initialdir = {BASE}
should_transfer_files = NO
getenv = False
request_cpus = 2
request_memory = 4 GB
request_disk = 4 GB
requirements = (OpSysAndVer == "AlmaLinux9") && (TARGET.FileSystemDomain == "{a.filesystem_domain}")
log = {out}/condor.$(ClusterId).log
output = {out}/condor.$(ClusterId).$(ProcId).out
error = {out}/condor.$(ClusterId).$(ProcId).err
periodic_hold = (JobStatus == 2) && ((CurrentTime - JobCurrentStartExecutingDate) > 10800)
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
notification = Never
queue 1
'''
(out/'job.sub').write_text(submit)
subprocess.run(['bash','-n',str(out/'worker.sh')],check=True)
subprocess.run(['condor_submit','-dry-run',str(out/'dryrun.classad'),str(out/'job.sub')],check=True)
print('PREPARED',out)
print('SUBMIT: condor_submit',out/'job.sub')
