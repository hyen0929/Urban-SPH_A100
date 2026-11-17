# 2026-09-10_leipzig-full-bcfix

> 2026-09-10 17:52:14 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Full run A (transient LES + LDM). CASE_ID 1002 경계 수정본: L3/L4=1307.5/1317.5(유출 버퍼 4면), L_left/L_right=-+648.75(y 주기, 주기 1297.5=519면x2.5), LDM 천장반사 L3->L_ceiling, 횡풍 초기속도 제거(INIT_CROSSFLOW_UY=0). 검증런에서 버퍼 8면이 주기 상과 |du|max=0.000e+00 일치, 양 측면 mean_ux 대칭(+9.9%), 유출단 연속, 클립 0개 확인. 설정: 스핀업 0~150 s(도메인 통과시간 153 s 1회분), 평균 150~300 s -> t=300 s 에 mean_field.bin(Run B 입력), 종료 300 s. LDM 13 소스 x 100,000 = 1,300,000, z=5 m, 150~300 s 연속방출(점당 666.7개/s). 소스별 p_type=101~113 - KERNEL_clc_concentration(functions_LDM.cuh:1078)이 p_type 일치 입자끼리만 농도를 합산하므로 소스별 분리에 필수. 유체 규약(0=경계,1=유체,>=1000=IBM)과 겹치지 않게 101부터 부여. 13개 전부 격자 정합·실외 확인(3D 벽거리 5.05~58.85 m). 예상 60,000 step x 4.8 s = 80 h

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_surface_s10` |
| 출력 | `runs/2026-09-10_leipzig-full-bcfix/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=300.0` `-DINIT_CROSSFLOW_UY=0.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-10_leipzig-full-bcfix/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-10_leipzig-full-bcfix/env/build.sh && bash runs/2026-09-10_leipzig-full-bcfix/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-10_leipzig-full-bcfix` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
