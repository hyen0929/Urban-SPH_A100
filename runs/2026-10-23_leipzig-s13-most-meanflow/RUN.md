# 2026-10-23_leipzig-s13-most-meanflow

> 2026-10-23 22:46:41 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 2

## 1. 목적

1편 4.4용 Run B — 2026-10-20_leipzig-s13-most-wall 의 평균장 위에서 LDM 재계산. 전이 대 평균장 비교와 u* 검증(대수법칙 적합)에 쓸 시간평균 유동장을 좌표와 함께 확보한다

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_wall` |
| 출력 | `runs/2026-10-23_leipzig-s13-most-meanflow/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=2.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_WALL_DEP=1` `-DLDM_MEAN_MODE=1` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DMEAN_FIELD_PATH="./runs/2026-10-20_leipzig-s13-most-wall/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=2 runs/2026-10-23_leipzig-s13-most-meanflow/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-23_leipzig-s13-most-meanflow/env/build.sh && bash runs/2026-10-23_leipzig-s13-most-meanflow/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-23_leipzig-s13-most-meanflow` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
