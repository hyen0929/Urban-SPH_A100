# 2026-11-11_cedval-c2-puff-dx15-lust

> 2026-11-11 11:27:07 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

MUSCL(RIEMANN_EULERV=2)은 클립 0 이지만 분해 난류를 지워 플룸이 좁아짐 → 중앙차분 75% + 선형 upwind 25%(RIEMANN_EULERV=3, UPWIND_BLEND=0.25)로 격자 잡음만 줄여 클립을 없애면서 난류·플룸 폭을 유지하는지 확인. 비교: 2026-11-03 중앙차분, 2026-11-07 MUSCL

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CEDVAL_C2_dx15_ext` |
| 출력 | `runs/2026-11-11_cedval-c2-puff-dx15-lust/out` |
| 빌드 플래그 | `-DRIEMANN_EULERV=3` `-DUPWIND_BLEND=0.25` `-DCASE_ID=1021` `-DIBM_length=0.5` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=1.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_Z_BOTTOM=150.0` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_DISABLE=1` `-DABL_ALPHA=0.27` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=1179.0` `-DCKPT_EVERY=5000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=500` `-DFULL_DUMP_EVERY=50000` `-DVORT_DUMP_EVERY=10000` `-DTRANSIENT_DUMP_EVERY=200` `-DTRANSIENT_DUMP_ZMAX=50.0` `-DTRANSIENT_DUMP_START=0.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-11-11_cedval-c2-puff-dx15-lust/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-11-11_cedval-c2-puff-dx15-lust/env/build.sh && bash runs/2026-11-11_cedval-c2-puff-dx15-lust/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-11-11_cedval-c2-puff-dx15-lust` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
