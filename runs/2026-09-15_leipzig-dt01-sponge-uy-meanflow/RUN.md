# 2026-09-15_leipzig-dt01-sponge-uy-meanflow

> 2026-09-15 13:42:48 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Run B (frozen mean-flow LDM). Run A 2026-09-15_leipzig-dt01-sponge-uy(dt 0.01, 스펀지 uz+uy, t=120 s 분기)가 t=299 s 에 덤프하는 mean_field.bin(평균창 150~299 s)으로 LDM 만 t=0~300 s 풀이. 입력·dt·LDM-frequency 동일 -> 방출 시각·step 일치, 같은 step 의 LDM_<step>stp.vtk 로 transient vs mean-flow 비교. env/chain.sh 가 Run A 정상 종료 + mean_field.bin 검증 후 GPU 4 에 자동 기동. 2026-09-15_leipzig-dt01-sponge-meanflow(미기동) 대체.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DMEAN_FIELD_PATH="./runs/2026-09-15_leipzig-dt01-sponge-uy/out/mean_field.bin"` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow/env/build.sh && bash runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS  exp/report.sh 가 생성 — 직접 고치지 말 것 (2026-09-16 15:22:36) -->

| 측정 | 값 |
|---|---|
| 상태 | **DONE** |
| 도달 step / 시뮬레이션 시간 | 30000 / 299.98 s |
| 계산 시간 (wall) | 3995.97 s (= 0.075 sim-s/wall-s) |
| 속도 tu_max (최대 / 최종) | 14.860 / 14.860 m/s |
| uz_max (최대) | 7.44 m/s |
| 온도 범위 tt_min–tt_max | 300.00 – 300.00 K |
| 출력 | 1872 vtk, 38G |
| 로그 이상 | 없음 |

<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
