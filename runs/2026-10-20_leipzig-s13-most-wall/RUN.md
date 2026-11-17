# 2026-10-20_leipzig-s13-most-wall

> 2026-10-20 16:14:58 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 1

## 1. 목적

1편 4장 기준 데이터셋 — 13소스 중립, MOST 지표 ON, 벽/지붕 침적 ON, 평균장 덤프. 기존 2026-09-15 런은 MOST 가 꺼져 있어 원고 2장 서술과 어긋나므로 이 런으로 대체한다. 4.2 regime / 4.3 침적경로 / 4.4 전이-평균장 / 4.5 보행자노출 전부 이 한 런에서 나온다

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_wall` |
| 출력 | `runs/2026-10-20_leipzig-s13-most-wall/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=2.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_WALL_DEP=1` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DCKPT_EVERY=1000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=1 runs/2026-10-20_leipzig-s13-most-wall/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-20_leipzig-s13-most-wall/env/build.sh && bash runs/2026-10-20_leipzig-s13-most-wall/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-20_leipzig-s13-most-wall` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
