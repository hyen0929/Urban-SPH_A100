# 2026-09-23_neutral-meanflow-localust

> 2026-09-23 17:25:36 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

중립 Run B 를 국소 u* 로 재실행 — 기존 중립 Run B 는 mean_field.bin 에 저장된 하드코딩 us_update=0.5 를 썼고 안정 Run B 는 0.2399 를 써서 두 런의 침적을 나란히 놓을 수 없었다. MOST 밴드를 열어 DD_local_most 가 각 런의 동결 평균장에서 국소 u* 를 계산하게 한다.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-23_neutral-meanflow-localust/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=2.5` `-DMEAN_FIELD_PATH="./runs/2026-09-15_leipzig-dt01-sponge-uy/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-23_neutral-meanflow-localust/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-23_neutral-meanflow-localust/env/build.sh && bash runs/2026-09-23_neutral-meanflow-localust/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-23_neutral-meanflow-localust` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
