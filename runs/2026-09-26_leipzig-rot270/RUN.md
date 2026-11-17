# 2026-09-26_leipzig-rot270

> 2026-09-26 13:27:09 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

풍향 민감도 270°: 도시를 270° 회전시킨 중립 전이 LES. 건물(IBM 마커)만 도메인 중심 (657.5,0) 기준 회전하고 유체·지면·천장은 미변경. 평균창 150~299s 로 기준런 2026-09-15_leipzig-dt01-sponge-uy 및 90°·180° 와 동일. LDM 출력 100 step. GPU 4 (A100, sm_80). 목적: lambda_p-발자국 관계의 풍향 일반성 검증 (0/90/180/270 네 방향 완성).

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_rot270` |
| 출력 | `runs/2026-09-26_leipzig-rot270/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` `-DLDM_DUMP_EVERY=100` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DCKPT_EVERY=2000` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-26_leipzig-rot270/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-26_leipzig-rot270/env/build.sh && bash runs/2026-09-26_leipzig-rot270/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-26_leipzig-rot270` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
