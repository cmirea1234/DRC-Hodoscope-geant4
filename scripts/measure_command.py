"""Run a command once and record wall time and child rusage without GNU time."""
import argparse
import json
import resource
import subprocess
import time
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
p.add_argument('command', nargs=argparse.REMAINDER)
a = p.parse_args()
cmd = a.command[1:] if a.command and a.command[0] == '--' else a.command
assert cmd
assert not a.output.exists(), a.output
started = time.time()
start = time.monotonic()
returncode = None
error = None
try:
    returncode = subprocess.call(cmd)
except BaseException as exc:
    error = repr(exc)
    raise
finally:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    a.output.write_text(json.dumps(dict(command=cmd, started_unix=started,
        wall_seconds=time.monotonic()-start, user_cpu_seconds=usage.ru_utime,
        system_cpu_seconds=usage.ru_stime, max_child_rss_kib=usage.ru_maxrss,
        rss_note='Linux child rusage peak; not summed concurrent-process peak. Use HTCondor for whole-job memory.',
        returncode=returncode, error=error), indent=2)+'\n')
raise SystemExit(returncode if returncode >= 0 else 128-returncode)
