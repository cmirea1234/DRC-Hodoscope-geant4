"""Consolidate completed diagnostics without claiming cosmic equivalence."""
from pathlib import Path
import json
base=Path(__file__).resolve().parents[1]
c=base/'runs/knu_diagnostic_20261029'
p=next(c.glob('profile_2188880_0_*/profile_comparison.json'))
t=json.loads(p.read_text())
lines=['원인 및 시간 진단 1차 보고 — 2026-10-29','',
'결론: primary 샘플링이 장시간의 원인은 아님. H1의 전달 차이는 단면 내',
'공통 광자 출발점 시험에서도 재현됨. 현재 light32를 동등한 빠른 대체물로',
'선택할 근거가 부족하며, 원통 기준 모델의 절대 정확성도 이 시험으로 확정되지 않음.','',
'1. 완료된 500트리거 결과의 공통 선택 473개',
'H1 생성/등록 광자: original 2,217,535 / light32 2,253,381.',
'H1 SiPM 도달: 114,275 / 71,181 (light32 약 37.7% 감소).',
'H1 도달/생성: 5.153% / 3.159%. 도달 후 검출 비율: 16.407% / 16.717%.',
'원인은 단순 생성량/PDE 변경보다 전달 단계 차이에 있다는 근거.',
'공통 생성 106,234개에서 T1&T2는 original 495, light32 500, 공통 473.',
'original의 나머지 선택 5개는 이 공통 생성 구간 뒤에 있음.',
'선택 조건 차이와 광학 응답 차이를 구분해야 함.','',
'2. 공통 입력 광자 대조 시험',
'H1X 32채널 순환, 450 nm, 전방 반구 방향, 같은 위치/방향/파장/채널을',
'ROOT에서 정확히 대조함. 각 조건 3200개. 실제 cosmic/beam 통계가 아닌 진단.',
'SiPM 바로 앞 출발: 원통/32면 모두 도달 3181, 검출 554로 일치.',
'중심선 출발: 도달 원통 187, 32면 193, 64면 191.',
'단면의 반경 0.47 mm 원판에 균일 출발: 도달 원통 322, 32면 199, 64면 227.',
'원통에서 도달하지만 32면에서 도달하지 못한 125개는 모두 world_escape.',
'이 중 56개는 active exit를 통과하지 못했고 69개는 통과 후 손실.',
'반대 방향(32면만 도달)은 2개. 대응 표본 도달률 차이 -3.844±0.346 pp (1sigma).',
'전체 광자 평균 step 수: 원통 81.79, 32면 11.02, 64면 15.08.',
'같은 입력에서도 전달/반사 경로와 생존 기간이 달라짐. 면수 증가로 일부',
'회복되지만 64면에서도 차이가 남아 있음. 입사 위치에 따른 표면 근사/경계',
'처리 민감성이 유력하며, 정확한 첫 분기 경계/법선 원인은 아직 미확정.',
'중심선 시험만으로 모델을 검증하면 이 차이를 놓칠 수 있음.','',
'3. 단계별 시간 (별도 timer-only 사본, 동일 노드, seed 261026304, 5트리거)',
'모든 카테고리는 배타적 wall time. 소규모 표본이며 500트리거 시간 분해가 아님.',
'process_other는 계측 범위 밖의 Geant4 초기화/프레임워크/관리/정리 등을 포함함.',
'nonoptical_tracking에는 에너지 침적과 섬광 secondary 생성 비용이 포함됨.',
'stacking_registration은 주로 생성 광자 lineage 등록과 스택 분류.',
'event_end_adc_and_rows에는 별도 계측 waveform/audit를 제외한 ADC/출력 작업 포함.',
'run_end_write_close만 전체 I/O 비용으로 보지 말 것.','']
keys=list(t['original']['timings'])
lines.append('항목 | original 초 | light32 초')
for k in keys:
 lines.append(f"{k} | {t['original']['timings'][k]['exclusive_wall_seconds']:.6f} | {t['light32']['timings'][k]['exclusive_wall_seconds']:.6f}")
for m,d in t.items():
 a=d['baseline_resources']['wall_seconds'];b=d['instrumented_resources']['wall_seconds']
 lines.append(f'{m}: 기준 wall {a:.3f}s, 계측 wall {b:.3f}s, 변화 {(b/a-1)*100:.2f}%.')
 lines.append(f'{m}: events/channels/loss_summary/v21_configuration 모든 수치 branch 동일 확인.')
lines += ['',
'현재 skipOpticalUnlessT1T2=true: 비트리거 광자는 추적 시작 시 종료됨.',
'다만 생성/등록/스택/종료 관리와 비트리거 ADC 작업은 남아 있음.',
'이전 대화에서 비트리거 사건도 전체 광학 수송을 한다고 읽힐 수 있던 설명은 정정.',
'',
'4. 해석 및 다음 검증',
'원통-메시 차이는 단순 동일 결과의 속도 최적화가 아닐 수 있음.',
'다음 정밀 진단은 대표 불일치 광자의 첫 경계 분기, 입사각/법선/전반사,',
'active 구간과 fanout의 차이 및 단면 면수 수렴을 확인하는 것이 적절함.',
'현재 단계에서는 파라미터를 조절해 광량을 맞추거나 모델을 교체하지 않음.',
'5트리거의 단계별 시간 비중을 500트리거 전체에 그대로 적용하지 않음.',
'비트리거 광자 등록/관리 비용 최적화는 난수 소비와 물리 결과 동일성 검증이 필요.',
'',
'재현 및 보존',
'원본 models/macros/PAYLOAD_SHA256.json은 변경하지 않음. 46파일 검증 통과.',
'계측 소스: build-knu/profile_20261029/source 및 instrumentation.patch.',
'원판 광자 소스: build-knu/photon_disk_20261029/source 및 disk_source.patch.',
'완료 잡: profile 2188880 / centre-terminal photon 2188881 / disk photon 2188882.',
'결과/입력 매크로/해시/로그는 각 잡의 고유 폴더에 보존.',
'light64는 단면 면수 영향만 보는 대조이며 새 production 모델로 채택한 것이 아님.']
(c/'summary/DIAGNOSTIC_REPORT_KO.txt').write_text('\n'.join(lines)+'\n')
print(c/'summary/DIAGNOSTIC_REPORT_KO.txt')
