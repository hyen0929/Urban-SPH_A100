# 2026-09-30_leipzig-src29-rot180

> 2026-09-30 15:34:33 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

29 소스 × 풍향 180° 평균장 LDM. 기존 13 소스의 n 부족(ρs 신뢰구간 폭 0.89)을 보완하려고 소스를 29개로 늘렸다. 선정 기준: 회전 중심에서 r<=410 m(네 방향 모두 풍하여유 246 m 확보), 최소간격 120 m(λp(75) 창 겹침 10%), λp 8구간 균등 층화, 기존 13개 유지. 소스당 50,000 입자 — 수렴성 검증에서 지표가 이미 수렴한 값(runs/lambdap/pconv). 총 방출량은 질량 31.25 로 고정. 유동장은 각 방향 Run A 의 mean_field.bin 을 동결해 쓰므로 방향 간 유동 차이만 남는다. 목적: 형태-확산 상관이 풍향을 넘어 살아남는지 n=29 로 검정.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_src29_rot180` |
| 출력 | `runs/2026-09-30_leipzig-src29-rot180/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DLDM_DUMP_EVERY=5000` `-DMEAN_FIELD_PATH="./runs/2026-09-23_leipzig-rot180/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-30_leipzig-src29-rot180/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-30_leipzig-src29-rot180/env/build.sh && bash runs/2026-09-30_leipzig-src29-rot180/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-30_leipzig-src29-rot180` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
