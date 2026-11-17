# 2026-09-14_leipzig-dt01-meanflow

> 2026-09-14 16:04:11 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Run B (frozen mean-flow LDM). Run A 2026-09-14_leipzig-dt01-full 이 t=299 s 에 덤프하는 mean_field.bin(평균창 150~299 s)을 고정 유동장으로 읽어 LDM 만 t=0~300 s 풀이. Run A 와 같은 입력(dt 0.01, LDM-frequency 10, LDM 13 소스 t_release 150~300 s)이라 방출 시각·step 번호가 일치 -> 같은 step 의 LDM_<step>stp.vtk 로 transient vs mean-flow 플룸 비교. env/chain.sh 가 Run A 정상 종료 + mean_field.bin 검증 후 GPU 4 에 자동 기동. 이전 Run B(2026-09-10_leipzig-full-meanflow)는 Run A 정전 사망으로 미기동 -> 이것으로 대체.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-14_leipzig-dt01-meanflow/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DMEAN_FIELD_PATH="./runs/2026-09-14_leipzig-dt01-full/out/mean_field.bin"` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-14_leipzig-dt01-meanflow/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-14_leipzig-dt01-meanflow/env/build.sh && bash runs/2026-09-14_leipzig-dt01-meanflow/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-14_leipzig-dt01-meanflow` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**미기동·폐기** — 입력인 Run A(2026-09-14_leipzig-dt01-full)가 천장 불안정으로 t=139.5 s 에 중단. 2026-09-15_leipzig-dt01-sponge-meanflow 로 대체.
