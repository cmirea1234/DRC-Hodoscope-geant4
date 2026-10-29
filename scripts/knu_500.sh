#!/usr/bin/env bash
set -eo pipefail
repo=/u/user/chjoung/DRC-Hodoscope-geant4
cd "$repo"
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) echo 'Refusing UI host' >&2; exit 90;; esac
test "$#" = 4
cluster=$1
process=$2
target=$3
seed=$4
test "$target" = 500
test "$seed" = 261028500
campaign="$repo/runs/knu_500_20261028_1643"
attempt=$(mktemp -d "$campaign/job_${cluster}_${process}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; printf "%s\n" "$code" > "$attempt/WORKFLOW_EXIT_CODE"; echo "WORKFLOW_EXIT_CODE=$code"' EXIT
printf 'RESULT_DIRECTORY=%s\n' "$attempt"
source "$repo/scripts/knu_env.sh"
set -u
py="$repo/.venv/bin/python"
job_cpus=$(awk '$1 == "RequestCpus" {print $3;exit}' "$_CONDOR_JOB_AD")
test "$job_cpus" = 1
cp "$_CONDOR_JOB_AD" "$attempt/job.ad"
if [ -f "${_CONDOR_MACHINE_AD:-}" ]; then cp "$_CONDOR_MACHINE_AD" "$attempt/machine.ad"; fi
build="$repo/build-knu/job_2188629_0_bQsF6C"
previous="$repo/runs/knu_20261028_1512/job_2188629_0_bQsF6C"
sha256sum -c "$previous/binary_sha256.txt"
cp "$previous/binary_sha256.txt" "$build/CMakeCache.txt" "$knu_g4/applied_patches.json" "$attempt/"
for file in knu_env.sh knu_500.sh run.py validate_root.py verify_bundle.py measure_command.py summarize_knu_pair.py; do
 cp "$repo/scripts/$file" "$attempt/"
done
sha256sum scripts/{knu_env.sh,knu_500.sh,run.py,validate_root.py,verify_bundle.py,measure_command.py,summarize_knu_pair.py} > "$attempt/scripts.sha256"
{
 date -u --iso-8601=seconds
 hostname -f
 "$py" --version
 g++ --version
 geant4-config --version
 geant4-config --datasets
 geant4-config --check-datasets
 git rev-parse HEAD
 git status --short
 "$py" -m pip freeze
} > "$attempt/environment.txt"
for model in light32 original; do ldd "$build/hodoscope-$model"; done > "$attempt/build_linkage.txt"
if grep -q 'not found' "$attempt/build_linkage.txt"; then exit 91; fi
# Loader addresses can change; compare resolved library paths instead.
diff <(sed -E 's/ \(0x[0-9a-f]+\)//g' "$previous/build_linkage.txt") <(sed -E 's/ \(0x[0-9a-f]+\)//g' "$attempt/build_linkage.txt")
"$py" scripts/verify_bundle.py
for model in original light32; do
 out="$attempt/${model}_${target}"
 echo "MODEL_START model=$model target=$target seed=$seed host=$(hostname -f) time=$(date -u --iso-8601=seconds)"
 echo "SIMULATION_LOG=$out/run.log"
 "$py" scripts/measure_command.py --output "$attempt/${model}.resources.txt" -- \
  "$py" scripts/run.py "$model" --target "$target" --seed "$seed" \
  --build-dir "$build" --output "$out" --timeout-hours 12
 echo "VALIDATION_START model=$model"
 "$py" scripts/validate_root.py "$out" 2>&1 | tee "$attempt/${model}.validation.txt"
 echo "MODEL_PASS model=$model time=$(date -u --iso-8601=seconds)"
done
sha256sum -c "$attempt/scripts.sha256"
"$py" scripts/summarize_knu_pair.py "$attempt" --target "$target" --seed "$seed" \
 --output "$campaign/summary/comparison_${cluster}_${process}_$(basename "$attempt").json"
echo 'PAIR_PASS: both models completed 500 triggers and ROOT validation'
