# 2026-09-30_leipzig-ens3-src29

> 2026-09-30 15:58:19 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 2

## 1. 목적

난류 실현 변동의 크기 측정 + 형태-확산 상관의 시간 수렴. 29 소스를 같은 위치에서 세 시간창(150~300, 300~450, 450~600 s)에 나눠 방출해(p_type 101~187) 한 런 안에서 3-멤버 앙상블을 만든다. 같은 소스가 서로 다른 난류 구간을 세 번 겪으므로 그 흩어짐이 곧 실현 변동이며, 이것 없이는 풍향별 상관 차이(A90 폭 0.46)가 진짜인지 우연인지 판정할 수 없다. 종료 750 s 는 세 창 모두 나이 300 s 까지 확보되는 최소 시각이다(창1 t=600, 창2 t=750, 창3 t=750). 소스당 50,000 개는 수렴성 검증(runs/lambdap/pconv)에서 확인된 안전선. 부수 효과: 기존 150 s 추적으로는 발자국이 계속 자라 지표가 시간에 따라 흔들렸는데(180도에서 60초에 -0.52→-0.92) 600 s 까지 따라가 이를 해소한다. GPU 2 (H100, sm_90).

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_ens3_src29` |
| 출력 | `runs/2026-09-30_leipzig-ens3-src29/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DCKPT_EVERY=2000` `-DSLICE_DUMP_EVERY=250` `-DFULL_DUMP_EVERY=10000` `-DVORT_DUMP_EVERY=2500` `-DLDM_DUMP_EVERY=5000` |
| 실행 | `CUDA_VISIBLE_DEVICES=2 runs/2026-09-30_leipzig-ens3-src29/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-30_leipzig-ens3-src29/env/build.sh && bash runs/2026-09-30_leipzig-ens3-src29/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-30_leipzig-ens3-src29` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
