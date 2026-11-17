# 2026-09-26_leipzig-rot90-meanflow

> 2026-09-26 13:07:57 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

풍향 민감도 90° 평균장 LDM(Run B). runs/2026-09-23_leipzig-rot90 의 mean_field.bin 을 얼려 LDM 만 재계산하여 전이 대비 비정상성 기여를 분리한다. LDM 출력 100 step. GPU 4 (A100, sm_80) — 180° Run A 는 H100 이었으나 Run B 는 두 각도 모두 A100 으로 맞춰 각도 간 비교를 일치시켰다.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_rot90` |
| 출력 | `runs/2026-09-26_leipzig-rot90-meanflow/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DLDM_DUMP_EVERY=100` `-DLDM_MEAN_MODE=1` `-DMEAN_FIELD_PATH="./runs/2026-09-23_leipzig-rot90/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-26_leipzig-rot90-meanflow/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-26_leipzig-rot90-meanflow/env/build.sh && bash runs/2026-09-26_leipzig-rot90-meanflow/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-26_leipzig-rot90-meanflow` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
