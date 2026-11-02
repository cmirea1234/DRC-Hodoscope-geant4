# 처음 설치하고 original을 실행하기

이 안내는 Linux x86_64 기준입니다. 검증 환경은 AlmaLinux 9, GCC 11.3.0, CMake 3.30.6, Python 3.9.25, **Geant4 11.4.1 + 프로젝트 패치 3종**입니다. 다른 환경은 첫 시험부터 다시 확인합니다.

KNU 로그인 노드에서는 코드 받기·설정·환경 확인·제출만 합니다. 아래 Geant4/프로그램 **컴파일과 시뮬레이션은 계산 자원을 할당받은 노드에서** 실행합니다. 현재 작업 폴더에 이미 설치된 환경을 쓰는 경우에는 2절을 반복하지 말고 환경 파일을 준비한 뒤 3~4절로 갑니다.

## 1. 코드와 Python 준비

```bash
git clone --branch knu/initial-validation-20261028 https://github.com/cmirea1234/DRC-Hodoscope-geant4.git
cd DRC-Hodoscope-geant4
python3 scripts/verify_bundle.py
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-analysis-tested.txt
```

공백 없는 저장소 경로를 사용합니다. `requirements-analysis-tested.txt`는 검증한 직접 의존성 NumPy 2.0.2/uproot 5.6.9를 고정합니다. Python 3.9.25에서 검증했으며 모든 최신 Python에서 설치된다는 뜻은 아닙니다. 기존 `requirements-analysis.txt`는 허용 범위이고, 설치 시점마다 같은 버전을 고정하지는 않습니다.

컴파일러·CMake·Python/venv·Git·patch·tar가 필요합니다. KNU에서는 CVMFS GCC 11.3.0/CMake 3.30.6을 사용했습니다. 계정별 환경 모듈을 사용해도 되지만 Geant4와 프로그램을 호환되는 컴파일러/ABI로 빌드해야 합니다. `g++ --version`, `cmake --version`, `python3 --version`으로 실제 버전을 확인합니다.

## 2. 패치된 Geant4 준비 — 처음 설치할 때만

[공식 11.4.1 릴리스](https://geant4.web.cern.ch/news/2026-03-13-patch1141) 소스를 사용합니다. 이 저장소는 현재 최신 Geant4가 아니라 **검증한 11.4.1**에 고정돼 있습니다. 같은 버전 번호의 공용 설치본이라도 프로젝트 패치가 없으면 동일한 검증 환경이 아닙니다.

저장소 최상위에서 전용 소스/빌드/설치 디렉터리를 만듭니다. 공용 Geant4 설치에 패치를 덮어쓰지 않습니다.

```bash
mkdir -p build-deps
git clone --depth 1 --branch v11.4.1 https://gitlab.cern.ch/geant4/geant4.git build-deps/geant4-src

# 세 파일의 원본 해시와 패치 적용 가능 여부 확인
python3 scripts/apply_geant4_patches.py --source build-deps/geant4-src
# 확인한 전용 소스 사본에 적용
python3 scripts/apply_geant4_patches.py --source build-deps/geant4-src --apply
```

도구는 세 파일의 해시를 먼저 확인하고, 알 수 없는 버전이면 적용 전에 중단합니다. 적용 후 해시도 검사합니다. 세 패치는 G4UIparsing.cc, G4Voxelizer.cc, G4SurfBits.icc 수정이며 [manifest](../patches/manifest.json)가 원본/수정 해시를 정의합니다.

아래는 **계산 노드에서**, 할당 CPU가 2개 이상일 때의 빌드 예입니다. 컴파일러/CMake 환경을 먼저 설정합니다.

```bash
export G4PREFIX="$PWD/build-deps/geant4-install"
cmake -S build-deps/geant4-src -B build-deps/geant4-build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_INSTALL_PREFIX="$G4PREFIX" \
  -DGEANT4_BUILD_MULTITHREADED=OFF \
  -DGEANT4_USE_SYSTEM_EXPAT=OFF -DGEANT4_USE_SYSTEM_CLHEP=OFF \
  -DGEANT4_USE_SYSTEM_ZLIB=OFF -DGEANT4_USE_GDML=OFF \
  -DGEANT4_USE_HDF5=OFF -DGEANT4_USE_QT=OFF \
  -DGEANT4_USE_OPENGL_X11=OFF -DGEANT4_USE_RAYTRACER_X11=OFF \
  -DGEANT4_INSTALL_EXAMPLES=OFF -DGEANT4_INSTALL_DATA=ON
cmake --build build-deps/geant4-build -j2
cmake --install build-deps/geant4-build
cp build-deps/geant4-src/applied_patches.json "$G4PREFIX/applied_patches.json"
source "$G4PREFIX/bin/geant4.sh"
geant4-config --version
geant4-config --check-datasets
```

`GEANT4_INSTALL_DATA=ON`은 빌드에 필요한 데이터셋을 내려받아 설치하므로 네트워크/저장 공간이 필요합니다. 계산 노드에서 외부 다운로드가 안 된다면 데이터셋을 사전에 준비해야 합니다. KNU 검증에서는 `GEANT4_INSTALL_DATA=OFF`와 `GEANT4_INSTALL_DATADIR`로 공유 데이터셋 경로를 지정했습니다. 이 선택은 해당 경로에 아래 데이터가 실제로 존재하는 경우에만 사용합니다. [공식 빌드 옵션 안내](https://geant4.web.cern.ch/documentation/dev/ig_html/InstallationGuide/installguide.html)

검증한 데이터셋은 G4NDL 4.7.1, G4EMLOW 8.8, PhotonEvaporation 6.1.2, RadioactiveDecay 6.1.2, G4PARTICLEXS 4.2, G4PII 1.3, RealSurface 2.2, G4SAIDDATA 2.0, G4ABLA 3.3, G4INCL 1.3, G4ENSDFSTATE 3.0, G4CHANNELING 2.0입니다. `--check-datasets`에서 12종 모두 INSTALLED여야 합니다.

## 3. 환경 파일과 실행 전 검사

저장소의 `build-env.sh`에 **자신의 절대 경로**를 넣습니다. 이름이 `build*`이므로 Git에 포함되지 않습니다. 아래 `/절대경로/...`는 실제 설치 위치로 교체해야 합니다.

```bash
# build-env.sh 내용 예: 먼저 이 설치를 빌드한 컴파일러/CMake를 활성화
# source /자신의/컴파일러/setup.sh
export G4PREFIX=/절대경로/DRC-Hodoscope-geant4/build-deps/geant4-install
source "$G4PREFIX/bin/geant4.sh"
```

KNU 검증에서는 CVMFS 제공 GCC 11.3.0/CMake 3.30.6을 사용했습니다. 자신의 사이트가 제공하는 컴파일러 setup 또는 module 명령을 환경 파일 앞부분에 넣으세요. 사이트별 설치 경로는 관리자 안내를 사용합니다.

Geant4 설치를 만든 컴파일러 설정과 이 실행 환경을 맞춥니다. 새 셸/계산 노드에서도 환경 파일을 읽어야 합니다. [공식 설치 후 환경 설정](https://geant4.web.cern.ch/documentation/dev/ig_html/InstallationGuide/postinstall.html)

```bash
source ./build-env.sh
unset PYTHONPATH PYTHONHOME CPLUS_INCLUDE_PATH C_INCLUDE_PATH
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python scripts/check_runtime.py --patch-record "$G4PREFIX/applied_patches.json"
```

`ENVIRONMENT_PASS`는 버전, 패치 이력, 데이터셋, 분석 의존성 확인입니다. 라이브러리 자체가 해당 소스로 빌드됐다는 암호학적 증명이나 물리 검증은 아닙니다. 패치 이력은 위에서 실제 빌드한 설치에만 복사합니다.

## 4. KNU에서 첫 Condor 시험

이 새 진입점은 과거 날짜의 빌드나 ROOT 경로를 요구하지 않습니다. 현재 체크아웃과 `.venv`, `build-env.sh`를 사용하며 계산 노드에서 새로 빌드합니다. 저장소/설치/환경 파일은 계산 노드에서도 같은 경로로 접근할 수 있어야 합니다.

로그인 노드에서:

```bash
condor_status -af FileSystemDomain | sort -u
# 아래 YOUR_DOMAIN을 사용할 계산 노드의 공유 파일시스템 도메인으로 교체
python3 scripts/prepare_original_smoke.py --env-file ./build-env.sh --filesystem-domain YOUR_DOMAIN
```

화면에 출력된 `SUBMIT: condor_submit .../job.sub` 명령을 실행합니다. 준비 단계는 자동 제출하지 않습니다.

```bash
# 실제 출력된 경로와 번호를 사용
condor_submit runs/first_original_실제생성값/job.sub
condor_q 실제잡번호 -nobatch
```

잡은 CPU 2개/메모리 4 GiB를 요청해 original만 빌드하고, seed 261026304로 5트리거를 실행합니다. 초기화/해시/데이터셋 검사, ROOT 검증, 저장소 기준 집계 비교까지 통과해야 `FIRST_ORIGINAL_PASS`가 출력됩니다.

성공 기록은 `runs/first_original_*/start_*/FIRST_RUN_PASS.json`에 있습니다. 주요 기준은 **1,240개 생성, 5트리거, 생성광자 40,188, 검출 p.e. 360**입니다. 이 비교는 소규모 집계 회귀이며 ROOT 파일 전체의 바이트 비교가 아닙니다. 차이가 나면 고통계로 확대하지 말고 환경과 입력 차이를 조사합니다.

이 제출 파일은 KNU AlmaLinux9/shared filesystem 조건을 사용합니다. 다른 사이트에서는 스케줄러·파일 전송·환경 경로 정책을 맞추거나 아래 계산 노드용 명령을 사용합니다.

2026-11-02에 별도 체크아웃 경로에서 이 절차를 실제 시험했습니다. 새 original 빌드, 환경/데이터셋 검사, 5트리거 ROOT 검증과 기준 집계 비교가 통과했습니다. [시험 기록](../reference/first_use_original_20261102.json). 이 시험은 기존 Geant4/분석 환경을 재사용했으며 Geant4 전체 재설치·데이터 다운로드를 다시 수행한 것은 아닙니다.

## 5. 할당된 계산 노드에서 직접 첫 실행

KNU 로그인 노드에서 실행하는 명령이 아닙니다. `build-env.sh`를 불러오고 3절의 환경 검사를 통과한 뒤 실행합니다.

```bash
cmake -S . -B build-original -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$G4PREFIX" -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build build-original --target hodoscope-original -j2
ctest --test-dir build-original --output-on-failure
.venv/bin/python scripts/run.py original --target 5 --seed 261026304 \
  --build-dir build-original --output runs/first_original_5
.venv/bin/python scripts/validate_root.py runs/first_original_5
```

출력 폴더가 이미 있으면 덮어쓰지 않습니다. 재시험은 새 이름을 사용합니다. CTest는 기준 묶음 해시 검사이고, 실제 ROOT 검증은 마지막 명령입니다. 기준 수치는 [original_small.json](../reference/original_small.json)을 참고합니다.

## 6. 이후 통계 확대와 모듈 확장

뮤온 한 잡의 `--target`은 최대500입니다. 서로 다른 seed의 여러 잡으로 통계를 늘리고 `(seed,event)`를 사건 식별자로 씁니다. [고통계 실행 원칙](MUON_CAMPAIGNS_KO.md)을 따르세요. 처음 받은 환경에서는 첫 실행이 통과한 빌드를 고정하고 사이트의 제출 경로·자원·seed를 새로 정합니다. 기존 `knu_*` 캠페인 스크립트는 과거 서버 경로를 참조하는 이력용 도구이므로 그대로 재제출하지 않습니다.

모듈 회전 확장은 별도 소스를 생성합니다:

```bash
python3 scripts/build_original_layout_source.py --output build-layout-source
# 아래 빌드는 계산 노드에서
cmake -S build-layout-source -B build-layout -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$G4PREFIX"
cmake --build build-layout -j2
```

`hodoscope-original-layout`은 `KNU_LAYOUT_FILE`로 전달되는 pose 파일을 읽습니다. `prepare_knu_beam.py`는 schema2 JSON에서 pose/매크로를 준비할 수 있지만 시뮬레이션은 계산 노드에서 실행해야 합니다. [배치 사용법](../configs/README_KO.txt)의 기존 `knu_*` 검증/제출 도구는 과거 KNU 경로를 참조하므로 새 환경에서는 경로와 기준 실행 자료를 준비해야 합니다. 회전된 배치의 cosmic 타기팅은 지원하지 않습니다.

## 문제가 생기면

| 증상 | 먼저 확인 |
| --- | --- |
| Geant4Config.cmake를 못 찾음 / 버전 불일치 | `G4PREFIX`, `geant4-config --version`, CMAKE_PREFIX_PATH |
| 데이터셋 누락 | `source "$G4PREFIX/bin/geant4.sh"`, `geant4-config --check-datasets` |
| 공유 라이브러리 또는 GLIBCXX 오류 | Geant4/프로그램 컴파일러 호환성과 `ldd 실행파일` |
| 패치 해시 불일치 | 정확한 11.4.1 소스인지 확인. 오류를 무시하거나 manifest를 맞춰 바꾸지 않음 |
| Condor Held | `condor_q 잡번호 -af HoldReason`, 제출 폴더 `.err`/`.out`, `workflow.log` |
| 큐에서 사라짐 | `condor_history 잡번호 -af JobStatus ExitCode`와 결과 PASS 확인 |
| 준비만 했는데 ROOT가 없음 | prepare 단계는 매크로/제출 파일만 작성. 실제 제출/실행 여부 확인 |
