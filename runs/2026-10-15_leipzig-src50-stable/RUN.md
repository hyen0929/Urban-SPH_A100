# 2026-10-15_leipzig-src50-stable

> 2026-10-15 15:39:31 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 3

## 1. 목적

안정 성층에서도 소스 간 침적 편차가 유지되는지 — 중립 50소스(runs/2026-10-06_leipzig-src50-rot0)와 짝. u* 처리를 양쪽 모두 균일(폴백)로 두어 성층 효과만 남긴다

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_src50_stable_rot0` |
| 출력 | `runs/2026-10-15_leipzig-src50-stable/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DMEAN_FIELD_PATH="./runs/2026-09-21_leipzig-most-stable/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=3 runs/2026-10-15_leipzig-src50-stable/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-15_leipzig-src50-stable/env/build.sh && bash runs/2026-10-15_leipzig-src50-stable/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-15_leipzig-src50-stable` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
