#!/usr/bin/env bash
set -eo pipefail
cd /u/user/chjoung/DRC-Hodoscope-geant4
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
source scripts/knu_env.sh
attempt=$(mktemp -d "$PWD/runs/knu_layout_20261029/job_${1}_${2}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
build="$PWD/build-knu/$(basename "$attempt")"
cmake -S build-knu/original_layout_v2_20261029/source -B "$build" -DCMAKE_BUILD_TYPE=Release -DGeant4_DIR="$Geant4_DIR"
cmake --build "$build" -j2
ln -s hodoscope-original-layout "$build/hodoscope-original"
echo "$build" > "$attempt/build_path.txt"
.venv/bin/python scripts/test_knu_layout.py "$attempt" "$build"
