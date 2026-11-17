# 2026-10-07_leipzig-src50w-rot270

> 2026-10-07 15:46:15 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

50 소스 × 풍향 270° 평균장 LDM, 벽면침적 ON. 2026-10-06_leipzig-src50-rot270(벽면 off)와 소스·유동장이 동일하고 DD_WALL_DEP=1 만 다르다 — 벽면침적 효과를 같은 소스집합에서 직접 비교하기 위함(09-17 벽면 런은 13소스였다). solv 에 buffer_type 출력을 추가해 0=공중/1=지면/2=벽면·지붕 을 직접 구분한다. DD_wall_dist 기본 2.0 m, 반사 임계 기본(0.5*dx_sp=1.25 m). 소스당 50,000 입자, LDM 출력 1 s.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_src50w_rot270` |
| 출력 | `runs/2026-10-07_leipzig-src50w-rot270/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DDD_WALL_DEP=1` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DLDM_DUMP_EVERY=100` `-DMEAN_FIELD_PATH="./runs/2026-09-26_leipzig-rot270/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-10-07_leipzig-src50w-rot270/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-07_leipzig-src50w-rot270/env/build.sh && bash runs/2026-10-07_leipzig-src50w-rot270/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-07_leipzig-src50w-rot270` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
