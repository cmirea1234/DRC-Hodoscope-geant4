#!/usr/bin/env bash
set -eo pipefail
cd /u/user/chjoung/DRC-Hodoscope-geant4
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
source scripts/knu_env.sh
attempt=$(mktemp -d "$PWD/runs/knu_beam_20261029/job_${1}_${2}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" scripts/knu_beam_smoke.sh scripts/prepare_knu_beam.py scripts/validate_knu_beam.py "$attempt/"
.venv/bin/python scripts/prepare_knu_beam.py --layout configs/layouts/baseline.json --beam configs/beams/pencil_3gev_smoke.json --output "$attempt/original_electron_10"
out="$attempt/original_electron_10"
binary="$PWD/build-knu/job_2188629_0_bQsF6C/hodoscope-original"
sha256sum "$binary" > "$attempt/binary_sha256.txt"
echo "BEAM_START output=$out"
.venv/bin/python - "$out" "$binary" <<'PY'
import subprocess,sys
from pathlib import Path
out,binary=sys.argv[1:];base=Path.cwd()
with (Path(out)/'run.log').open('w') as log:
 subprocess.run([str(base/'.venv/bin/python'),str(base/'scripts/measure_command.py'),'--output',str(Path(out)/'resources.json'),'--',binary,'run.mac'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
PY
.venv/bin/python scripts/validate_knu_beam.py "$out"
echo BEAM_SMOKE_PASS
