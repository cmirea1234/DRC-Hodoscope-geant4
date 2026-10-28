KNU 이전용: 경량화 전 원본 + 32각형 경량 모델

포함 모델은 original과 light32 두 개뿐이다.
현재 진행 중인 서버 계산을 옮기거나 중단하는 패키지가 아니다.
Git 저장소 주소는 아직 지정되지 않았으며 push는 하지 않았다.
비공개 저장소를 권장한다. data/runtime에는 실험 기반 gain/noise/threshold/template가 포함된다.
실험 raw, 대용량 ROOT, 과거 중간 결과, 실행 바이너리와 Geant4 설치본은 포함하지 않는다.

소스
models/original: 비교 시험에 사용 중인 경량화 전 원본 사본.
models/light32: 동결된 경량화 소스. 실행 macro는 32각형 하나만 제공.
광학 C++ 소스/헤더는 원본과 바이트 단위로 동일하다.
light32 소스에는 기존 다른 표현을 지원하는 분기가 남아 있지만 이를 삭제/변형하지 않았다.
StraightFibres.cc도 링크에 필요한 원본 코드이므로 포함한다. 직선 실행 모델을 제공하는 것은 아니다.
SOURCE_PROVENANCE.json에 출처와 해시를 기록한다.
CMake와 실행기는 이식용으로 새로 작성했고, 옛 서버의 systemd 실행기는 포함하지 않는다.
macros는 중첩된 설정을 풀어낸 것이며 출력 경로·seed·목표 수는 실행기가 별도로 지정한다.

환경
Geant4 11.4.1, C++17, CMake, Python3를 사용한다.
Geant4 데이터셋도 같은 버전에 맞춰 설치하고 geant4.sh 환경을 먼저 불러온다.
분석/검산은 requirements-analysis.txt의 numpy와 uproot를 사용한다.
시뮬레이션 ROOT 출력 자체는 Geant4 analysis가 담당하므로 CERN ROOT 설치가 필수는 아니다.
현재 서버는 별도 Geant4 수정 3종과 ASan/UBSan 진단 빌드를 사용 중이다.
patches/*.patch와 manifest.json은 그 수정 내역이다. 원본 소스 해시를 확인한 뒤
별도 Geant4 소스 사본에만 patch -p1 --dry-run으로 확인하고 적용/재빌드한다.
공용 Geant4 설치에 덮어쓰지 않는다. 이 패치가 모든 충돌을 해결했다는 뜻은 아니다.
KNU의 실제 compiler, Geant4, 배치 시스템은 아직 확인하지 않았다.

빌드 (저장소 최상위에서, 사용 가능한 Geant4 환경을 먼저 설정한 뒤)
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build -j2
  ctest --test-dir build --output-on-failure
CTest는 패키지 해시/설정 무결성 검사이지 전체 물리 Gate 재검증이 아니다.
이번 CMake는 Qt를 요구하지 않는다. 기존 main/physics는 그대로이며 batch macro로 실행한다.
GUI 실행기를 별도로 제공하지 않는다.

같은 입력으로 첫 소규모 실행 (허가된 계산 노드에서)
  python3 scripts/run.py original --target 5 --seed 261026304 --output runs/original_5
  python3 scripts/validate_root.py runs/original_5
  python3 scripts/run.py light32 --target 5 --seed 261026304 --output runs/light32_5
  python3 scripts/validate_root.py runs/light32_5
기존 output 폴더가 있으면 덮어쓰지 않고 실패한다. 재실행은 새 이름으로 원인 검토 후 실시한다.
run.py --prepare-only를 추가하면 물리 실행 없이 macro와 입력만 준비한다.
run.py는 프로세스 정상 종료만 표시하며 별도 validate_root.py 통과 전까지 검증 완료가 아니다.
작업당 한 Geant4 프로세스/한 코어이다. 병렬화는 다른 seed의 작업 여러 개를 배치 시스템에 제출한다.
한 작업은 최대 500개 트리거로 제한한다. 총 통계는 독립 seed 묶음을 추가해 늘린다.
KNU 로그인 노드에서 대규모 계산을 실행하지 않는다.

소규모 기준 (현재 서버, seed261026304, T1/T2 5개)
original: 생성 1240, 선택 사건 생성광자 40188, 검출 p.e. 360.
light32: 생성 1240, 선택 사건 생성광자 33509, 검출 p.e. 185.
reference/*.json은 실제 정상 종료 기록이다. 시간은 서로 다른 서버의 속도 비교 기준이 아니다.
컴파일러/빌드/플랫폼 변화로 차이가 나면 먼저 조사한다. 수치를 맞추는 보정은 하지 않는다.
그 후 모델별 500개와 공통 primary/광학 response 검산을 거쳐 통계를 확대한다.
KNU에서의 빌드와 실행은 아직 검증되지 않았다.

물리 조건
콘크리트/암막천 없음; source z=150 mm.
T1 높이의 목표 영역 반폭 16 mm에서 역투영하는 target_backproject.
modified_gaisser_joint, 0.2~1000 GeV, 천정각 <=45도.
H2 각 평면 사각 파이버 29개, H1 offset=(-0.25,-0.25) mm.
광생성/광수송/PDE/파형/노이즈/채널별 threshold, X/Y AND 유지.
T1/T2 각각 Edep>=0.05 MeV가 분모. 동일 primary의 직접 통과만 세는 strict 분모가 아니다.
raw 효율을 검출기 고유 효율로 단정하지 않는다.

Git
.gitignore는 build/runs/ROOT/로그/바이너리/키/환경 파일을 제외한다.
Git remote URL과 비공개 여부를 결정한 뒤 포함 파일을 검토하고 push한다.
Git commit 번호 외에 데이터/소스 해시와 실제 실행 라이브러리도 기록해야 한다.
PAYLOAD_SHA256.json은 처음 제공한 묶음 기준이다. 나중에 코드 수정 시에는 새 변경 기록과
검증을 남기고 manifest를 명시적으로 갱신해야 한다. 해시 오류를 무시하고 실행하지 않는다.
