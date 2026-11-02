# DRC Hodoscope Geant4

**기준 모델: original · 기준 환경: Geant4 11.4.1 + 저장소 패치 3종**

2026-11-02 확인: original의 **뮤온 5,000 T1&T2 트리거(500 × 10잡)** 실행과 결과 검증을 완료했습니다. 현재는 original을 계산 기준으로 사용합니다. light32는 광학 전달 차이가 확인된 비교·개발용 모델입니다.

여기서 검증 완료는 실행 안정성, 기준 결과 재현, ROOT/채널/광자 집계의 일관성을 뜻합니다. **실제 검출기 데이터와의 물리적 일치나 모든 배치의 정확성을 인증한 것은 아닙니다.** 빔테스트 확대는 보류한 상태입니다.

## 처음 받았다면

1. [설치부터 첫 실행까지](docs/INSTALL_KO.md)를 따라 환경을 준비합니다.
2. `python3 scripts/verify_bundle.py`로 기준 소스 묶음을 확인합니다.
3. 계산 노드에서 original **5트리거**를 실행하고 결과를 검증합니다.
4. 첫 실행을 통과한 뒤 통계/배치를 변경합니다.

현재 변경 사항이 있는 브랜치를 명시해서 받으세요:

```bash
git clone --branch knu/initial-validation-20261028 https://github.com/cmirea1234/DRC-Hodoscope-geant4.git
cd DRC-Hodoscope-geant4
python3 scripts/verify_bundle.py
```

Git에는 소스·설정·검증 요약이 들어 있고, Geant4 설치본·바이너리·`.venv`·전체 ROOT 파일은 들어 있지 않습니다. 새 환경에서는 설치/빌드가 필요합니다. 최상위 CMake는 **Geant4 11.4.1을 정확히 요구**하며, 최신 버전으로 임의 교체하면 동일 검증 환경이 아닙니다.

## 무엇을 사용할 수 있나

| 기능 | 현재 상태 | 안내 |
| --- | --- | --- |
| original, 기존 배치의 우주선 뮤온 | 현재 기준. 5,000트리거/10잡 모두 정상 종료·ROOT 검증 통과 | [검증 범위](docs/VALIDATION_KO.md) |
| 독립 seed의 뮤온 Condor 잡 분산·합산 | 기존 환경에서 10잡 검증. 사이트별 자원/경로/seed 설정 필요 | [뮤온 실행 안내](docs/MUON_CAMPAIGNS_KO.md) |
| 새 체크아웃에서 original 첫 빌드/5트리거 | 별도 환경 파일로 KNU Condor 잡 준비. 과거 빌드/ROOT 경로 불필요 | [첫 실행 안내](docs/INSTALL_KO.md#4-knu에서-첫-condor-시험) |
| T1/T2/H1/H2 전체 모듈 이동·3D 회전 | original 기반 별도 확장. 기본 회귀 및 대표 변환 시험 통과 | [배치 설정](configs/README_KO.txt) |
| 명시적 3 GeV pencil 전자빔 | 10개 기본 배치 및 1개 회전 배치 동작 시험만 완료. 확대 보류 | [검증 범위](docs/VALIDATION_KO.md) |
| light32 | 원통형 original과 광학 전달이 달라 대체 기준으로 채택하지 않음 | [모델 선택 근거](docs/VALIDATION_KO.md) |

회전한 배치의 기존 cosmic 타기팅은 지원하지 않습니다. 실제 새 배치의 모듈 간 간섭, 빔 폭·발산, 자기장, 실측 데이터 비교는 추가 작업입니다.

## 확인된 실행 환경

| 항목 | 검증 환경 |
| --- | --- |
| 운영체제 | Linux x86_64, AlmaLinux 9 계산 노드 |
| Geant4 | **11.4.1**, [패치 3종](patches/manifest.json), Release, multithreading OFF |
| 컴파일러 / C++ | GCC 11.3.0 / C++17 |
| CMake | 3.30.6 |
| Python | 3.9.25 |
| 분석 | NumPy 2.0.2, uproot 5.6.9 — [직접 의존성 버전](requirements-analysis-tested.txt) |
| 데이터셋 | Geant4 11.4.1에 맞는 12종. 설치 안내의 검사로 확인 |
| 실행 방식 | 잡당 단일 Geant4 프로세스. 통계 확대는 독립 seed의 여러 잡 |

CERN ROOT나 Qt GUI는 이 batch 실행에 필요하지 않습니다. 다른 OS·컴파일러·라이브러리 조합은 동일하게 검증됐다고 가정하지 않습니다.

## 모듈 위치와 각도

[modules_baseline.json](configs/layouts/modules_baseline.json)을 복사해서 수정합니다.

```json
"H1": {
  "position_mm": [0, 0, 40],
  "rotation_deg": [10, 0, 0]
}
```

H1 기준점을 월드 좌표 Z=40 mm로 옮기고 X축으로 10도 기울입니다. 섬유·연결부·SiPM이 함께 움직입니다. 위치는 이동량이 아닌 절대 좌표이고, 회전은 기존 방향에 고정 월드축 X→Y→Z 순으로 추가합니다. 빔 좌표는 별도입니다. 자세한 기준점/월드 크기/지원 범위는 [배치 사용법](configs/README_KO.txt)을 보세요.

이 기능은 `hodoscope-original-layout`을 별도로 빌드해야 합니다. 기본 `hodoscope-original`에 JSON만 전달한다고 적용되지는 않습니다.

## 결과 읽기

5,000트리거 실행은 2026-10-29 17:21 KST에 완료됐습니다. 생성 뮤온 1,103,858개 중 목표 5,000트리거를 수집했고, 10개 잡 모두 종료 코드 0 및 ROOT 검증 PASS입니다. 네 평면 조건부 응답률은 35.18%, Wilson 95% 구간 33.87~36.51%입니다. 분모는 **T1&T2 트리거**이며 검출기 고유 효율이 아닙니다.

- [공유 가능한 검증 요약 JSON](reference/knu_muons_5000_summary.json): 독립 seed별 시간·집계·ROOT 해시
- [검증 범위와 남은 작업](docs/VALIDATION_KO.md)
- [뮤온 잡 확인과 합산 원칙](docs/MUON_CAMPAIGNS_KO.md)
- [상세 작업 이력](WORK_PLAN_KO.txt)

큐에서 사라진 것만으로 성공을 판단하지 않습니다. 종료 코드와 `VALIDATION.json`, 캠페인의 `SUCCESS.json`/`SUMMARY.json`을 함께 확인합니다. ROOT 이벤트는 `(seed,event)`로 구분하며 서로 다른 잡의 event 번호를 같은 사건으로 취급하지 않습니다.

[README_KO.txt](README_KO.txt)는 최초 이전 시점의 동결 기록입니다. 그 안의 “KNU 미검증” 등 상태 설명은 현재 상태가 아닙니다. 현재 사용 안내는 이 README와 `docs/`를 기준으로 합니다.
