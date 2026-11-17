# 2026-10-20_leipzig-s13-wall

> 2026-10-20 15:35:04 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 1

## 1. 목적

1편 4.3용 — 실도시 전이 런에서 지면/벽·지붕 침적 분리. 기존 13소스 Run A(2026-09-15_leipzig-dt01-sponge-uy)와 유동 비트동일하며 DD_WALL_DEP 만 추가. 평균장도 덤프해 4.4의 Run B 를 자체 재현

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_wall` |
| 출력 | `runs/2026-10-20_leipzig-s13-wall/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_WALL_DEP=1` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DCKPT_EVERY=1000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=1 runs/2026-10-20_leipzig-s13-wall/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-20_leipzig-s13-wall/env/build.sh && bash runs/2026-10-20_leipzig-s13-wall/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-20_leipzig-s13-wall` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
