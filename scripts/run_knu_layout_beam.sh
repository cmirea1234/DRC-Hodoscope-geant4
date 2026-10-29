#!/usr/bin/env bash
set -eo pipefail
repo=/u/user/chjoung/DRC-Hodoscope-geant4
cd "$repo"
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
source scripts/knu_env.sh
set -u
campaign=$1
binary=$2
attempt=$(mktemp -d "$campaign/start_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
.venv/bin/python scripts/prepare_knu_beam.py --layout "$campaign/layout.json" --beam "$campaign/beam.json" --output "$attempt/run"
out="$attempt/run"
export KNU_LAYOUT_FILE="$out/layout.pose"
sha256sum "$binary" > "$attempt/binary_sha256.txt"
echo "LAYOUT_BEAM_START output=$out"
.venv/bin/python - "$out" "$binary" <<'PY'
from pathlib import Path
import subprocess,sys
base=Path.cwd();out=Path(sys.argv[1]);binary=sys.argv[2]
with (out/'run.log').open('w') as log:
 subprocess.run([str(base/'.venv/bin/python'),str(base/'scripts/measure_command.py'),'--output',str(out/'resources.json'),'--',binary,'run.mac'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=12*3600)
PY
.venv/bin/python scripts/validate_knu_beam.py "$out"
echo LAYOUT_BEAM_PASS
