"""Snapshot layout/beam JSON and submit one serial original-layout beam job."""
import argparse,json,subprocess,tempfile,shutil
from pathlib import Path
from layout_config import validate_layout
BASE=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--layout',type=Path,required=True);p.add_argument('--beam',type=Path,required=True);p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
validate_layout(json.loads(a.layout.read_text()))
binary=BASE/'build-knu/original-layout-validated/hodoscope-original-layout'
if not a.prepare_only and not binary.is_file():raise SystemExit('Validated layout executable not yet installed')
root=BASE/'runs/layout_beam';root.mkdir(exist_ok=True)
campaign=Path(tempfile.mkdtemp(prefix='submission_',dir=root))
for src,name in ((a.layout,'layout.json'),(a.beam,'beam.json')):shutil.copy2(src,campaign/name)
subprocess.run(['python3',str(BASE/'scripts/prepare_knu_beam.py'),'--layout',str(campaign/'layout.json'),'--beam',str(campaign/'beam.json'),'--output',str(campaign/'preflight')],check=True)
submit=f'''universe = vanilla
executable = {BASE}/scripts/run_knu_layout_beam.sh
arguments = {campaign} {binary.resolve()}
initialdir = {BASE}
should_transfer_files = NO
getenv = False
request_cpus = 1
request_memory = 4 GB
request_disk = 4 GB
requirements = (OpSysAndVer == "AlmaLinux9") && (TARGET.FileSystemDomain == "knu.ac.kr")
log = {campaign}/condor.$(ClusterId).log
output = {campaign}/condor.$(ClusterId).$(ProcId).out
error = {campaign}/condor.$(ClusterId).$(ProcId).err
notification = Never
periodic_hold = (JobStatus == 2) && ((CurrentTime - JobCurrentStartExecutingDate) > 46800)
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
queue 1
'''
(campaign/'job.sub').write_text(submit)
subprocess.run(['condor_submit','-dry-run',str(campaign/'dryrun.classad'),str(campaign/'job.sub')],check=True)
print('SNAPSHOT',campaign,flush=True)
if not a.prepare_only:
 r=subprocess.run(['condor_submit',str(campaign/'job.sub')],check=True,text=True,capture_output=True);(campaign/'submission.txt').write_text(r.stdout);print(r.stdout)
