#!/usr/bin/env bash
set -eo pipefail
cd /u/user/chjoung/DRC-Hodoscope-geant4
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
source scripts/knu_env.sh
attempt=$(mktemp -d "$PWD/runs/knu_diagnostic_20261029/disk_${1}_${2}_XXXXXX")
exec > >(tee "$attempt/workflow.log") 2>&1
trap 'code=$?; echo "$code" > "$attempt/WORKFLOW_EXIT_CODE"' EXIT
hostname -f > "$attempt/host.txt"
cp "$_CONDOR_JOB_AD" scripts/knu_photon_disk.py scripts/knu_photon_disk.sh "$attempt/"
.venv/bin/python scripts/verify_bundle.py
build="$PWD/build-knu/$(basename "$attempt")"
cp build-knu/photon_disk_20261029/disk_source.patch build-knu/photon_disk_20261029/source_sha256.json "$attempt/"
cmake -S build-knu/photon_disk_20261029/source -B "$build" -DCMAKE_BUILD_TYPE=Release -DGeant4_DIR="$Geant4_DIR"
cmake --build "$build" -j2
sha256sum "$build"/hodoscope-* > "$attempt/binary_sha256.txt"
.venv/bin/python scripts/knu_photon_disk.py "$attempt/cases" "$build"
