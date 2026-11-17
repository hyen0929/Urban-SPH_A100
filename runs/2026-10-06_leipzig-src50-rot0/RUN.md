# 2026-10-06_leipzig-src50-rot0

> 2026-10-06 17:36:03 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

50 소스 × 풍향 0° 평균장 LDM. 29소스(n=29)에서 소스를 50개로 더 늘려 검정력을 올린다. 선정: 회전중심 r<=410 m, 최소간격 80 m, 방출운(반경 15 m) 열린비율 >= 0.90, λp 8구간 균등 층화, 기존 13개 유지. 소스당 50,000 입자(총 250만), 질량 31.25 로 총 방출량 고정. 유동장은 각 방향 Run A 의 mean_field.bin 동결. LDM 출력 1 s(100 step). 주의: 간격 80 m 라 λp(75) 창이 35% 겹친다 — 명목 n=50 이지만 유효표본은 그보다 작다.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_src50_rot0` |
| 출력 | `runs/2026-10-06_leipzig-src50-rot0/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DLDM_DUMP_EVERY=100` `-DMEAN_FIELD_PATH="./runs/2026-09-15_leipzig-dt01-sponge-uy/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-10-06_leipzig-src50-rot0/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-06_leipzig-src50-rot0/env/build.sh && bash runs/2026-10-06_leipzig-src50-rot0/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-06_leipzig-src50-rot0` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
