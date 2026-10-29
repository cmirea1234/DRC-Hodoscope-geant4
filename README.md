# DRC Hodoscope Geant4 — 시작 안내

현재 빔테스트 기준 모델은 **original**입니다. T1·T2·H1·H2를 각각 이동·회전하는 설정을 추가했습니다. 기존 형상·재료·광학 조건은 보존하고, 모듈의 섬유·연결부·SiPM을 함께 배치합니다.

이 문서는 2026-10-29 KNU 검증 이후의 안내입니다. [README_KO.txt](README_KO.txt)는 최초 이전 당시 기록으로, 그 안의 “KNU에서 아직 미검증” 등 상태 설명은 당시 기준입니다.

## 먼저 볼 파일

| 목적 | 파일 |
| --- | --- |
| 좌표·회전 규칙, 실행 방법, 제한 사항 | [configs/README_KO.txt](configs/README_KO.txt) |
| 기본 모듈 배치 | [configs/layouts/modules_baseline.json](configs/layouts/modules_baseline.json) |
| 3 GeV 전자 10개 입력 | [configs/beams/pencil_3gev_smoke.json](configs/beams/pencil_3gev_smoke.json) |
| 전체 배치 회전 예제 | [configs/layouts/modules_rigid_example.json](configs/layouts/modules_rigid_example.json) |
| 위 회전 예제와 짝인 전자 1개 입력 | [configs/beams/pencil_3gev_rigid_example.json](configs/beams/pencil_3gev_rigid_example.json) |

## 위치와 각도 변경

기본 배치 JSON을 복사하고 원하는 모듈의 두 항목을 바꿉니다. 예:

```json
"H1": {
  "position_mm": [0, 0, 40],
  "rotation_deg": [10, 0, 0]
}
```

H1 기준점을 월드 좌표 Z=40 mm로 옮기고 X축으로 10도 기울입니다.

- `position_mm`: 이동량이 아니라 모듈 기준점의 **최종 절대 좌표** `[X,Y,Z]`, 단위 mm.
- `rotation_deg`: original의 기존 방향에 추가하는 회전. 고정 월드축 X→Y→Z 순서, 단위 도.
- `world_half_size_mm`: 월드 각 축의 반길이. 모듈 전체와 빔 발생점을 포함해야 합니다.
- 빔 좌표·방향은 별도 JSON에서 설정합니다. 모듈을 옮겨도 빔이 자동으로 따라가지 않습니다.

배치 JSON만 바꿀 때는 재컴파일할 필요가 없습니다. 예제 배치는 실제 설치 예정 좌표를 뜻하지 않습니다.

## 현재 KNU 작업 환경에서 실행

아래 명령은 저장소 최상위에서 실행합니다. 로그인 노드는 준비·제출·상태 확인에 사용하고, 빌드와 시뮬레이션은 Condor 계산 노드에서 수행합니다.

```bash
# 기존 원본 묶음의 해시/설정 확인 — 새 확장 전체의 물리 검증은 아님
python3 scripts/verify_bundle.py

# 입력과 제출 파일 준비만
python3 scripts/submit_knu_layout_beam.py \
  --layout configs/layouts/modules_baseline.json \
  --beam configs/beams/pencil_3gev_smoke.json \
  --prepare-only

# 실제 Condor 제출: 위 명령에서 --prepare-only를 제외
python3 scripts/submit_knu_layout_beam.py \
  --layout configs/layouts/modules_baseline.json \
  --beam configs/beams/pencil_3gev_smoke.json
```

준비/제출할 때마다 새 `runs/layout_beam/submission_*` 폴더를 만듭니다. 출력된 잡 번호로 `condor_q 잡번호 -nobatch`를 확인합니다. 큐에서 사라졌다는 것만으로 성공을 판단하지 않습니다.

결과는 해당 제출 폴더의 `start_*/run/` 아래에 있습니다.

- `INPUT.json`: 사용한 배치·빔 설정과 입력 해시
- `layout_placements.jsonl`: 부품별 최종 월드 위치·회전
- `run.log`, `result.root`: 실행 로그와 결과
- `VALIDATION.json`: 결과 검증 `PASS` 여부
- `resources.json`: 실행 시간과 종료 코드

현재 빔 제출 도구는 pencil 전자 1~100개를 지원합니다. 사건 수는 **생성 입자 수**이며 트리거 목표 수가 아닙니다.

## 다른 경로/서버에서 처음 받았을 때

코드만 받은 상태에서는 위 시뮬레이션 제출을 바로 실행할 수 없습니다. Git에는 `build*/`, `.venv/`, `runs/`가 포함되지 않아 설치 환경·검증된 바이너리·과거 ROOT 결과는 별도로 준비해야 합니다.

필요 환경은 Geant4 11.4.1과 데이터셋, 프로젝트의 [Geant4 패치 기록](patches/manifest.json), C++17/CMake 및 [분석 의존성](requirements-analysis.txt)입니다. 현재 KNU 도구에는 `/u/user/chjoung/DRC-Hodoscope-geant4` 경로와 사이트 조건이 들어 있습니다. 다른 경로에서는 `scripts/knu_env.sh`, 배치 실행 셸, `condor/*.sub`의 경로/조건을 검토해야 합니다.

모듈 확장은 다음 순서로 빌드·검증합니다.

1. `scripts/build_original_layout_source.py --output 새_소스_경로`로 original 기반 확장 소스를 생성합니다.
2. `scripts/knu_layout_check.sh`의 소스 경로를 맞추고 `condor/knu_layout_check.sub`로 계산 노드에서 빌드·검증합니다.
3. 이 회귀검사는 이전 original 5트리거 ROOT 결과를 참조합니다. 새 환경에서는 기준 실행 결과를 준비하고 `scripts/test_knu_layout.py`의 참조 경로를 맞춥니다.
4. 검증 통과한 빌드를 `build-knu/original-layout-validated` 링크로 연결한 뒤 빔 제출 도구를 사용합니다.

최상위 CMake만 빌드하면 기존 original/light32가 만들어집니다. 모듈 회전 확장 실행 파일 `hodoscope-original-layout`은 위 별도 소스 생성 절차가 필요합니다.

## 확인한 범위

KNU 잡 2188886에서 기본 5트리거 결과의 `events/channels/loss_summary` 비교 열이 기존 original과 정확히 일치했습니다. 전체 동시 변환과 모듈별 독립 변환 각각 554개 부품의 위치·회전, 246개 코어/센서/트리거 위치 및 광자 전달 시험을 통과했습니다.

잡 2188888에서는 배치와 빔을 함께 회전한 3 GeV 전자 1개를 새 제출 도구로 실행하여 ROOT 검증을 통과했습니다. 이는 동작 확인이며 효율이나 최종 빔테스트 물리를 검증한 것은 아닙니다.

새 배치는 명시적 빔/고정 광자 입력을 지원합니다. 회전한 배치의 기존 cosmic 타기팅, building, 추가 트리거는 지원 범위 밖입니다. 실제 설치 좌표를 정하면 모듈 간 겹침과 빔 경로를 별도 확인해야 합니다.

확장 코드는 [extensions/layout](extensions/layout)에 있고, 원본 `models/`, `macros/`는 보존했습니다. 상세 이력은 [WORK_PLAN_KO.txt](WORK_PLAN_KO.txt)에 있습니다.
