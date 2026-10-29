#!/usr/bin/env bash
# Identical private Release Geant4 for both models; no system install changes.
source /cvmfs/sft.cern.ch/lcg/releases/gcc/11.3.0/x86_64-el9/setup.sh
export PATH="/cvmfs/sft.cern.ch/lcg/releases/CMake/3.30.6-34165/x86_64-el9-gcc11-opt/bin:/cvmfs/sft.cern.ch/lcg/releases/git/2.29.2-daa17/x86_64-el9-gcc11-opt/bin:$PATH"
knu_g4=/u/user/chjoung/DRC-Hodoscope-geant4/build-knu/geant4-11.4.1-patched
test -f "$knu_g4/BUILD_COMPLETE" || { echo "Private Geant4 build is incomplete" >&2; return 1; }
source "$knu_g4/bin/geant4.sh"
export Geant4_DIR="$knu_g4/lib/cmake/Geant4"
if [ ! -f "$Geant4_DIR/Geant4Config.cmake" ]; then export Geant4_DIR="$knu_g4/lib64/cmake/Geant4"; fi
export CMAKE_PREFIX_PATH="$knu_g4"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME CPLUS_INCLUDE_PATH C_INCLUDE_PATH
