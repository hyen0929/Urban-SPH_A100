# 2026-11-07_cedval-c2-puff-dx15-muscl

> 2026-11-07 03:21:15 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 3

## 1. 목적

중앙차분 모멘텀 이류에서 첫 두 줄·마지막 줄 건물 앞 지면층에 속도 상한(20 m/s) 입자 ~1700 개(MOST 끈 분기에서도 그대로) → RIEMANN_EULERV=2(MUSCL minmod upwind)로 t=0부터 다시 돌려 클립이 사라지는지, 캐노피 난류가 살아있는지, 퍼프 dosage 가 나아지는지 확인. 60 퍼프 완결 위해 t=1180

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CEDVAL_C2_dx15_ext` |
| 출력 | `runs/2026-11-07_cedval-c2-puff-dx15-muscl/out` |
| 빌드 플래그 | `-DCASE_ID=1021` `-DRIEMANN_EULERV=2` `-DIBM_length=0.5` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=1.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_Z_BOTTOM=150.0` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_DISABLE=1` `-DABL_ALPHA=0.27` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=1179.0` `-DCKPT_EVERY=5000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=500` `-DFULL_DUMP_EVERY=50000` `-DVORT_DUMP_EVERY=10000` `-DTRANSIENT_DUMP_EVERY=200` `-DTRANSIENT_DUMP_ZMAX=50.0` `-DTRANSIENT_DUMP_START=0.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=3 runs/2026-11-07_cedval-c2-puff-dx15-muscl/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-11-07_cedval-c2-puff-dx15-muscl/env/build.sh && bash runs/2026-11-07_cedval-c2-puff-dx15-muscl/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-11-07_cedval-c2-puff-dx15-muscl` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->

클립은 0 건(최대 |u| 10.8)이지만 분해 난류가 거의 사라졌다. 배열 안 TKE 는 중앙차분(2026-11-03)의 1/10–1/100 이다. 퍼프 7개 비교(같은 퍼프, 각 런별 상수 k):
- 감쇠비 117→477 m 는 6.7(실험 6.1), 225→477 m 는 3.5(실험 3.4)로 개선됐다.
- 플룸 가장자리(±90 m)가 실험보다 크게 낮다.
- 퍼프 간 CV 는 0.20(실험 0.54)이다.
- 전체 FAC2 0.64, NMSE 1.95 로 중앙차분 dx1.5(0.66 / 1.05)보다 나쁘다.

그림은 analysis/figs_cmp3/ 에 있다 (스크립트 figs3.py, turb_prof.py).
다음: 논문 본결과는 중앙차분 dx1.5 로 하고, MUSCL 은 수치감쇠 민감도로 쓴다. 혼합(upwind 비율 작게)은 선택 사항이다.
