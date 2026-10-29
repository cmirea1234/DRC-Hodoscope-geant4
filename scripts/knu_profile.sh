#!/usr/bin/env bash
set -eo pipefail
repo=/u/user/chjoung/DRC-Hodoscope-geant4
cd "$repo"
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
campaign="$repo/runs/knu_diagnostic_20261029"
attempt=$(mktemp -d "$campaign/profile_${1}_${2}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"; echo "WORKFLOW_EXIT_CODE=$code"' EXIT
source scripts/knu_env.sh
set -u
py="$repo/.venv/bin/python"
source_dir="$repo/build-knu/profile_20261029/source"
build="$repo/build-knu/$(basename "$attempt")"
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
cp scripts/knu_profile.sh build-knu/profile_20261029/instrumentation.patch build-knu/profile_20261029/source_sha256.json "$attempt/"
hostname -f > "$attempt/host.txt"
"$py" scripts/verify_bundle.py
cmake -S "$source_dir" -B "$build" -DCMAKE_BUILD_TYPE=Release -DGeant4_DIR="$Geant4_DIR" -DPython3_EXECUTABLE="$py"
cmake --build "$build" -j2
sha256sum "$build"/hodoscope-* > "$attempt/binary_sha256.txt"
# Repeat a previously validated seed on the same worker with and without timers.
# Opposite model order partially checks order effects; these tiny runs are diagnostic.
for mode in baseline instrumented; do
 if [ "$mode" = baseline ]; then
  chosen="$repo/build-knu/job_2188629_0_bQsF6C"
  models='original light32'
 else
  chosen="$build"
  models='light32 original'
 fi
 for model in $models; do
  out="$attempt/${mode}_${model}"
  echo "PROFILE_START mode=$mode model=$model output=$out time=$(date -u --iso-8601=seconds)"
  "$py" scripts/measure_command.py --output "$attempt/${mode}_${model}.resources.json" -- \
   "$py" scripts/run.py "$model" --target 5 --seed 261026304 --build-dir "$chosen" --output "$out" --timeout-hours 1
  "$py" scripts/validate_root.py "$out" | tee "$attempt/${mode}_${model}.validation.json"
 done
done
"$py" scripts/check_knu_profile.py "$attempt" | tee "$attempt/comparison.txt"
echo PROFILE_PASS
