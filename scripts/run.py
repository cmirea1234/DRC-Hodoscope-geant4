"""One isolated batch job; no scheduler assumptions and no silent retry."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from verify_bundle import verify, BASE

parser=argparse.ArgumentParser()
parser.add_argument('model',choices=['original','light32'])
parser.add_argument('--seed',type=int,default=261026304)
parser.add_argument('--target',type=int,default=5)
parser.add_argument('--build-dir',type=Path,default=BASE/'build')
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--prepare-only',action='store_true')
parser.add_argument('--timeout-hours',type=float,default=12)
args=parser.parse_args()
assert 0 < args.target <= 500 and 0 < args.seed < 2000000000, 'Use 1..500 triggers per job; scale with independent seed jobs'
verify()
out=args.output.resolve();build=args.build_dir.resolve()
assert not any(c.isspace() for c in str(out)), 'Use an output path without whitespace'
out.mkdir(parents=True,exist_ok=False)
for p in (BASE/'data/runtime').iterdir():
    if p.is_file():shutil.copy2(p,out/p.name)
macro=(BASE/'macros'/f'{args.model}.mac').read_text()+f'''
/random/setSeeds {args.seed} {args.seed+17}
/hodo/cosmicSeed {args.seed}
/hodo/targetCoincidences {args.target}
/hodo/diagnosticRootFile result
/run/initialize
/process/inactivate Cerenkov
/geometry/navigator/check_mode true
/geometry/navigator/push_notify true
/run/beamOn 20000000
/control/echo KNU_TWO_MODEL_JOB_COMPLETE
'''
(out/'run.mac').write_text(macro)
exe=build/f'hodoscope-{args.model}'
record=dict(model=args.model,seed=args.seed,target=args.target,status='PREPARED',
            macro_sha256=hashlib.sha256(macro.encode()).hexdigest(),
            root='result.root',started=time.time(),build=str(build))
try:
    record['git_commit']=subprocess.check_output(['git','-C',str(BASE),'rev-parse','HEAD'],stderr=subprocess.DEVNULL,text=True).strip()
except subprocess.CalledProcessError:record['git_commit']='not_committed'
def save():
    p=out/'STATUS.json.tmp';p.write_text(json.dumps(record,indent=2)+'\n');p.replace(out/'STATUS.json')
save()
if args.prepare_only:
    print('PREPARED',out);raise SystemExit(0)
proc=None
diagnostic=re.compile(r'G4Exception-START|Stuck Track|GeomNav\d|GeomVol\d|AddressSanitizer|UndefinedBehaviorSanitizer|runtime error:|Segmentation fault|COMMAND NOT FOUND|Illegal application state')
try:
    assert exe.is_file(), f'Build first: {exe}'
    record.update(status='RUNNING',binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    linkage=subprocess.run(['ldd',str(exe)],capture_output=True,text=True,check=True).stdout
    (out/'library_linkage.txt').write_text(linkage)
    with (out/'run.log').open('w') as log:
        proc=subprocess.Popen([str(exe),'run.mac'],cwd=out,stdout=log,stderr=subprocess.STDOUT)
        record['pid']=proc.pid;save()
        start=time.monotonic();suffix=''
        with (out/'run.log').open(errors='replace') as reader:
            while True:
                text=suffix+reader.read()
                assert not diagnostic.search(text), 'Runtime diagnostic: inspect run.log'
                suffix=text[-256:]
                code=proc.poll()
                if code is not None:
                    # Drain output once more after process exit.
                    tail=reader.read();assert not diagnostic.search(suffix+tail)
                    assert code==0, f'process exited {code}'
                    break
                assert time.monotonic()-start<args.timeout_hours*3600,'job timeout'
                time.sleep(1)
    assert 'KNU_TWO_MODEL_JOB_COMPLETE' in (out/'run.log').read_text(errors='replace')
    record.update(status='PROCESS_COMPLETE_ROOT_VALIDATION_REQUIRED',elapsed_seconds=time.time()-record['started'])
except BaseException as exc:
    if proc is not None and proc.poll() is None:
        proc.terminate()
        try:proc.wait(timeout=10)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()
    record.update(status='FAILED',error=repr(exc));raise
finally:save()
print(record['status'],out)
