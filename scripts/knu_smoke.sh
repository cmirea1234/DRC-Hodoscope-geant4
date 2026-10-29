#!/usr/bin/env bash
set -eo pipefail
repo=/u/user/chjoung/DRC-Hodoscope-geant4
cd "$repo"
# This wrapper must be dispatched by HTCondor, never run on the login host.
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) echo 'Refusing UI host' >&2; exit 90;; esac
campaign="$repo/runs/knu_20261028_1512"
attempt="$campaign/job_${1}_${2}"
# A new start preserves prior artifacts even if the scheduler restarts this job.
attempt=$(mktemp -d "${attempt}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
cp "$repo/scripts/knu_env.sh" "$repo/scripts/knu_smoke.sh" "$attempt/"
source "$repo/scripts/knu_env.sh"
set -u
py="$repo/.venv/bin/python"
job_cpus=$(awk '$1 == "RequestCpus" {print $3;exit}' "$_CONDOR_JOB_AD")
test "$job_cpus" = 2
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
if [ -f "${_CONDOR_MACHINE_AD:-}" ]; then cp "$_CONDOR_MACHINE_AD" "$attempt/machine.ad"; fi
{
 date -u --iso-8601=seconds
 hostname -f
 id
 "$py" --version
 g++ --version
 cmake --version
 geant4-config --version
 geant4-config --datasets
 geant4-config --check-datasets
 git rev-parse HEAD
 git status --short
 "$py" -m pip freeze
} > "$attempt/environment.txt"
"$py" scripts/verify_bundle.py
"$py" -c 'import numpy,uproot; print("analysis dependencies",numpy.__version__,uproot.__version__)'
build="$repo/build-knu/$(basename "$attempt")"
cmake -S "$repo" -B "$build" -DCMAKE_BUILD_TYPE=Release -DGeant4_DIR="$Geant4_DIR" -DPython3_EXECUTABLE="$py" 2>&1 | tee "$attempt/configure.log"
"$py" "$repo/scripts/measure_command.py" --output "$attempt/build.resources.txt" -- cmake --build "$build" -j "$job_cpus" 2>&1 | tee "$attempt/build.log"
ctest --test-dir "$build" --output-on-failure 2>&1 | tee "$attempt/ctest.log"
cp "$build/CMakeCache.txt" "$attempt/CMakeCache.txt"
sha256sum "$build"/hodoscope-* > "$attempt/binary_sha256.txt"
cp "$knu_g4/applied_patches.json" "$attempt/applied_patches.json"
for exe in "$build"/hodoscope-*; do ldd "$exe"; done > "$attempt/build_linkage.txt"
if grep -q "not found" "$attempt/build_linkage.txt"; then exit 91; fi
result=0
for model in original light32; do
 out="$attempt/${model}_5"
 if "$py" "$repo/scripts/measure_command.py" --output "$attempt/${model}.resources.txt" -- "$py" scripts/run.py "$model" --target 5 --seed 261026304 --build-dir "$build" --output "$out" --timeout-hours 1; then
  "$py" scripts/validate_root.py "$out" > "$attempt/${model}.validation.txt" 2>&1 || result=1
 else
  result=1
 fi
done
printf '%s\n' "$result" > "$attempt/WORKFLOW_EXIT_CODE"
exit "$result"
