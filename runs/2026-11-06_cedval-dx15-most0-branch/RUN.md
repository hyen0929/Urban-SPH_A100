# 2026-11-06_cedval-dx15-most0-branch

> 2026-11-06 20:12:01 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

본 런(dx1.5)에서 첫 두 줄·마지막 줄 건물 바람받이 앞 지면 첫 층(z=0.75)에 속도 상한(20 m/s) 입자 ~1700 개가 소스·45 m 측정선을 덮음 → t=1000 체크포인트에서 MOST 만 끄고 30 s 돌려 클립이 사라지는지 확인

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CEDVAL_C2_dx15_br` |
| 출력 | `runs/2026-11-06_cedval-dx15-most0-branch/out` |
| 빌드 플래그 | `-DCASE_ID=1021` `-DIBM_length=0.5` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=1.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_Z_BOTTOM=150.0` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_DISABLE=1` `-DABL_ALPHA=0.27` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=999.0` `-DCKPT_EVERY=5000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=500` `-DFULL_DUMP_EVERY=50000` `-DVORT_DUMP_EVERY=10000` `-DTRANSIENT_DUMP_EVERY=100` `-DTRANSIENT_DUMP_ZMAX=50.0` `-DTRANSIENT_DUMP_START=1000.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 SOPHIA_RESTART=runs/2026-11-03_cedval-c2-puff-dx15/out/restart_t1000.bin runs/2026-11-06_cedval-dx15-most0-branch/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-11-06_cedval-dx15-most0-branch/env/build.sh && bash runs/2026-11-06_cedval-dx15-most0-branch/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-11-06_cedval-dx15-most0-branch` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

MOST 는 클립 원인이 아니다 — MOST 끈 뒤에도 클립 입자 2,100–2,500 개(본 런 1,500–2,000), 소스 부근 470–800 개로 그대로. 지면층 |u| 중앙값만 0.55→0.65 m/s.
체크포인트가 t=950(count 95000)이라 분기는 t=950.6 부터. t=1023.6 에서 수동 종료(결론이 나서).
다음: 모멘텀 이류 중앙차분 불안정으로 보고 RIEMANN_EULERV=2(MUSCL) 를 구현해 t=0 부터 재계산 → runs/2026-11-07_cedval-c2-puff-dx15-muscl
