#!/usr/bin/env bash
set -eo pipefail
cd /u/user/chjoung/DRC-Hodoscope-geant4
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
source scripts/knu_env.sh
attempt=$(mktemp -d "$PWD/runs/knu_diagnostic_20261029/photon_${1}_${2}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" scripts/knu_photon_probe.py scripts/knu_photon_probe.sh "$attempt/"
.venv/bin/python scripts/verify_bundle.py
.venv/bin/python scripts/knu_photon_probe.py "$attempt/cases"
