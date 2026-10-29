#!/usr/bin/env bash
set -eo pipefail
repo=/u/user/chjoung/DRC-Hodoscope-geant4
cd "$repo"
test -n "${_CONDOR_JOB_AD:-}"
case "$(hostname -s)" in cms|cms01|cms02|cms03) exit 90;; esac
record=$(mktemp -d "$repo/runs/knu_20261028_1512/geant4_${1}_${2}_XXXXXX")
exec > >(tee "$record/workflow.log") 2>&1
source /cvmfs/sft.cern.ch/lcg/views/LCG_108_geant4ext20250725/x86_64-el9-gcc11-opt/setup.sh
unset PYTHONPATH PYTHONHOME CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export PATH="/cvmfs/sft.cern.ch/lcg/releases/git/2.29.2-daa17/x86_64-el9-gcc11-opt/bin:$PATH"
set -u
py="$repo/.venv/bin/python"
work="$repo/build-knu/$(basename "$record")"
install="$repo/build-knu/geant4-11.4.1-patched"
test ! -e "$install"
mkdir -p "$work"
cp "$_CONDOR_JOB_AD" "$record/job.ad"
hostname -f > "$record/hostname.txt"
g++ --version > "$record/compiler.txt"
mkdir "$work/source"
printf '%s  %s\n' '99dcf5f9d4f806fb8c4fde85cb2674a42e4ca19833143464ff7efa55c1852140' "$repo/build-knu/geant4-v11.4.1.tar.gz" | sha256sum -c -
tar -xzf "$repo/build-knu/geant4-v11.4.1.tar.gz" -C "$work/source" --strip-components=1
"$py" - "$work/source" "$record" <<'PYCODE'
import hashlib,json,subprocess,sys
from pathlib import Path
source=Path(sys.argv[1]); record=Path(sys.argv[2]); manifest=json.loads(Path('patches/manifest.json').read_text()); results={}
for rel,item in manifest.items():
 p=source/rel; original=hashlib.sha256(p.read_bytes()).hexdigest(); assert original==item['original_sha256']
 patch=str((Path('patches')/item['patch']).resolve())
 subprocess.run(['patch','-p1','--batch','--dry-run','-i',patch],cwd=source,check=True)
 subprocess.run(['patch','-p1','--batch','-i',patch],cwd=source,check=True)
 changed=hashlib.sha256(p.read_bytes()).hexdigest(); assert changed==item['patched_sha256']
 results[rel]={'original_sha256':original,'patched_sha256':changed}
(record/'applied_patches.json').write_text(json.dumps(results,indent=2)+'\n')
PYCODE
cpus=$(awk '$1 == "RequestCpus" {print $3;exit}' "$_CONDOR_JOB_AD")
test "$cpus" = 8
cmake -S "$work/source" -B "$work/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX="$install" -DGEANT4_BUILD_MULTITHREADED=OFF -DGEANT4_USE_SYSTEM_EXPAT=OFF -DGEANT4_USE_SYSTEM_CLHEP=OFF -DGEANT4_USE_SYSTEM_ZLIB=OFF -DGEANT4_USE_GDML=OFF -DGEANT4_USE_HDF5=OFF -DGEANT4_USE_QT=OFF -DGEANT4_USE_OPENGL_X11=OFF -DGEANT4_USE_RAYTRACER_X11=OFF -DGEANT4_INSTALL_EXAMPLES=OFF -DGEANT4_INSTALL_DATA=OFF -DGEANT4_INSTALL_DATADIR=/cvmfs/geant4.cern.ch/share/data 2>&1 | tee "$record/configure.log"
"$py" "$repo/scripts/measure_command.py" --output "$record/build.resources.txt" -- cmake --build "$work/build" -j "$cpus" > "$record/build.log" 2>&1
cmake --install "$work/build" > "$record/install.log" 2>&1
cp "$work/build/CMakeCache.txt" "$record/CMakeCache.txt"
"$install/bin/geant4-config" --version > "$record/version.txt"
sha256sum "$install"/lib*/libG4*.so > "$record/library_sha256.txt"
cp "$record/applied_patches.json" "$install/applied_patches.json"
printf '%s\n' "$record" > "$install/BUILD_COMPLETE"
echo GEANT4_BUILD_COMPLETE
