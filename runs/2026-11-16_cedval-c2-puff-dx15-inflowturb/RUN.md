# 2026-11-16_cedval-c2-puff-dx15-inflowturb

> 2026-11-16 17:55:16 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 3

## 1. 목적

세 런 모두 퍼프 간 편차가 실험의 절반 이하이고 상류가 층류(유입 변동 0) → 유입 버퍼에 합성 난류(RFG, INFLOW_TURB=1, Iu=1/ln(z/z0) 상한 0.45, L=40 m)를 넣어 캐노피 난류·퍼프 간 편차·플룸 폭이 실험에 가까워지는지 확인. 이류는 중앙차분. 비교: 2026-11-03 중앙차분(유입 층류)

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CEDVAL_C2_dx15_ext` |
| 출력 | `runs/2026-11-16_cedval-c2-puff-dx15-inflowturb/out` |
| 빌드 플래그 | `-DRIEMANN_EULERV=0` `-DINFLOW_TURB=1` `-DCASE_ID=1021` `-DIBM_length=0.5` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=1.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_Z_BOTTOM=150.0` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_DISABLE=1` `-DABL_ALPHA=0.27` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=1179.0` `-DCKPT_EVERY=5000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=500` `-DFULL_DUMP_EVERY=50000` `-DVORT_DUMP_EVERY=10000` `-DTRANSIENT_DUMP_EVERY=200` `-DTRANSIENT_DUMP_ZMAX=50.0` `-DTRANSIENT_DUMP_START=0.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=3 runs/2026-11-16_cedval-c2-puff-dx15-inflowturb/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-11-16_cedval-c2-puff-dx15-inflowturb/env/build.sh && bash runs/2026-11-16_cedval-c2-puff-dx15-inflowturb/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-11-16_cedval-c2-puff-dx15-inflowturb` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
