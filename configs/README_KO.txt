original 모듈 배치 사용법 — 2026-10-29

1. 바꾸는 파일
configs/layouts/modules_baseline.json을 복사하여 사용한다.
T1, T2, H1, H2 각각의 position_mm와 rotation_deg만 바꾸면 된다.
위치는 [X,Y,Z] mm, 회전은 [X,Y,Z]도. 위치는 이동량이 아니라 회전 기준점의
최종 월드 좌표이다. 고정 월드축 X→Y→Z 순으로 추가 회전(Rz Ry Rx)을 적용한다.

기준점/기본 위치:
  T1: [-25,0,80]     트리거 중심
  T2: [0,-25,-75]    트리거 중심
  H1: [0,0,20]      H1 기준 평면 원점
  H2: [0,0,-20]     H2 기준 평면 원점
H1/H2는 X/Y 두 층과 fanout, cladding, grease, window, SiPM까지 같이 움직인다.
회전 [0,0,0]은 original의 기존 방향을 뜻한다. T2/H2의 기존 회전을 없애지 않는다.
예를 들어 H1의 position_mm=[0,0,40], rotation_deg=[10,0,0]이면
H1 기준점을 Z=40 mm로 옮기고 기준점 주위에서 X축으로 10도 기울인다.

world_half_size_mm은 월드 각 축의 반길이이다.
[500,500,2000]이면 전체 크기는 1 m × 1 m × 4 m.
빔 발생점과 모듈 전체가 월드 안에 있어야 한다. 부품 경계가 월드 밖으로 나가면
초기화에서 중단한다. 모듈 기준점만 안에 있는 것으로는 충분하지 않다.
빔은 월드 좌표이며 모듈을 옮겨도 자동으로 따라가지 않는다.

2. 빔과 제출
configs/beams/pencil_3gev_smoke.json은 기존 중앙 3 GeV 전자 10개 입력이다.
position_mm, direction(길이 1), generated_events, seed를 별도로 정한다.
생성 매크로는 /control/useDoublePrecision true로 소수 좌표 정밀도를 유지한다.
현재 제출 도구는 pencil 전자 1~100개를 지원한다. generated_events는 트리거 수가
아닌 생성 입자 수이다. 같은 seed 재실행은 독립 통계가 아니다.

로그인 노드에서 준비만:
  python3 scripts/submit_knu_layout_beam.py --layout configs/layouts/modules_baseline.json --beam configs/beams/pencil_3gev_smoke.json --prepare-only

Condor 제출:
  python3 scripts/submit_knu_layout_beam.py --layout configs/layouts/modules_baseline.json --beam configs/beams/pencil_3gev_smoke.json

터미널에 표시되는 cluster 번호로 condor_q 번호 -nobatch를 확인한다.
입력 사본/제출 파일/로그는 runs/layout_beam/submission_*/에 저장되며,
실행 결과는 그 안의 start_*/run/에 저장된다.
VALIDATION.json의 PASS와 resources.json으로 결과/시간을 확인한다.
prepare-only 후 실제 제출 명령을 실행하면 새 제출 폴더를 만든다.
계산과 빌드는 계산 노드에서만 수행한다.

3. 회전 예제
modules_rigid_example.json과 pencil_3gev_rigid_example.json은 모든 모듈과 빔을
함께 X=7,Y=11,Z=13도로 회전하고 [20,-10,30] mm 이동한 검증 예제이다.
실제 설치 예정 배치를 뜻하지 않는다. 두 파일을 짝으로 사용한다.
모듈별 독립 회전/이동도 검증했지만 실제 빔테스트 간격과 자석 조건은 미정이다.

4. 구현과 검증 범위
기준 models/original, macros와 payload 46개 파일은 보존한다.
extensions/layout/과 scripts/build_original_layout_source.py가 original의 별도
소스 사본을 만들고, 검증한 hodoscope-original-layout 실행 파일을 사용한다.
모듈 내부 형상·재료·광학 조건을 유지하며 최상위 부품의 배치만 변환한다.
중첩 부품은 부모와 함께 이동한다. 읽기 센서 좌표와 분석의 로컬 좌표/표면 통과
판정에도 같은 변환을 적용한다. 새 JSON을 바꿀 때마다 재컴파일할 필요는 없다.
기존 schema 1 baseline.json은 이전 실행 파일용이며 새 모듈 설정과 혼용하지 않는다.

검증 잡 2188886:
- 기본 5트리거에서 기존 original과 정수 집계/채널 신호 일치.
- 최종 빌드에서는 비교한 events/channels/loss_summary 열이 모두 정확히 일치.
  앞선 빌드는 일부 float 열에 약 10^-14의 차이가 있어 엄격 검사 실패를 기록했다.
- 전체 동시 변환 및 모듈별 독립 변환 각각 554개 부품 위치/회전을 별도 계산과 비교.
- 각 배치에서 코어/SiPM/트리거 246개 위치의 navigator 확인.
- 각 변환 배치에서 SiPM 직전 광자 32개 추적/손실 집계 확인.
현재 검사들은 모든 부품 표면의 겹침을 빠짐없이 검증하는 검사가 아니다.
새 설치 좌표를 정하면 모듈 간 간섭 및 실제 빔 경로를 별도 검증해야 한다.
변경 배치는 명시적 빔 또는 고정 광자 입력을 지원하며, 기존 cosmic 타기팅은
회전 배치에 대응하지 않아 차단한다. building/추가 트리거/lite 모델도 범위 밖이다.

5. 결과 해석과 개발
INPUT.json은 실제 배치/빔과 입력 해시를 보관한다.
layout_placements.jsonl은 부품별 최종 월드 위치/회전을 기록한다.
기존 ROOT의 일부 Z/yaw 설정 열은 original의 명목 설정값이다. 새 배치를 읽을 때
그 열만 사용하지 말고 INPUT.json 및 layout_placements.jsonl을 함께 사용한다.
validation 출력을 켠 경우 validation_geometry_map/validation_trigger_map은 변환을
반영한다. v21_geometry_map의 source 좌표/치수는 원래 모듈 내부의 명목 값이다.
Cerenkov 비활성, 기존 노이즈/PDE/threshold 조건을 유지한다.
모든 입사 사건을 기록하고 광학 수송하므로 입사 수와 트리거 수를 구별할 수 있다.
run.log의 [muon energy]는 기존 고정 라벨이다. 입자는 /gun/particle e-로 설정되며
현재 ROOT events에는 PDG 열이 없어 입자는 매크로/생성기 설정으로 확인한다.

구현을 바꿀 경우 scripts/build_original_layout_source.py --output 새_소스_경로로
새 소스를 만들고 scripts/knu_layout_check.sh의 소스 경로를 갱신한다.
condor/knu_layout_check.sub를 제출하여 회귀/변환 검증을 통과한 빌드만
build-knu/original-layout-validated 링크 대상으로 지정한다.
