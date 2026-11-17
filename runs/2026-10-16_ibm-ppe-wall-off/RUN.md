# 2026-10-16_ibm-ppe-wall-off

> 2026-10-16 13:33:16 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 1

## 1. 목적

건물 IBM 마커를 PPE·압력힘에서 제외한 현행 기준선. runs/2026-09-15_leipzig-dt01-sponge-uy 의 t=299.98s 체크포인트에서 분기해 20s(2000 step) 재주행. IBM_PPE_WALL=1 짝과 초기장·플래그·GPU 가 모두 같으므로 둘의 차이는 전부 벽 압력 BC 효과.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_wallbc` |
| 출력 | `runs/2026-10-16_ibm-ppe-wall-off/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=500` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=1 SOPHIA_RESTART=runs/2026-09-15_leipzig-dt01-sponge-uy/out/restart.bin runs/2026-10-16_ibm-ppe-wall-off/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-16_ibm-ppe-wall-off/env/build.sh && bash runs/2026-10-16_ibm-ppe-wall-off/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS  exp/report.sh 가 생성 — 직접 고치지 말 것 (2026-10-16 19:36:14) -->

| 측정 | 값 |
|---|---|
| 상태 | **NEVER_RAN** |
| 도달 step / 시뮬레이션 시간 | 32000 / 320.00 s |
| 계산 시간 (wall) | 9802.12 s (= 0.033 sim-s/wall-s) |
| 속도 tu_max (최대 / 최종) | 20.000 / 20.000 m/s |
| uz_max (최대) | 20.00 m/s |
| 온도 범위 tt_min–tt_max | 300.00 – 300.00 K |
| 출력 | 128 vtk, 14G |
| 로그 이상 | 없음 |

**PPE_diag_explicit.csv** (PPE 잔차 진단)

```
샘플 200개
bi_max      평균 120.3  최대 134.1  최종 111.9
resid_max   평균 81.3  최대 90.54  최종 62.97
div_un1_avg 평균 0.06061  최대 0.06189  최종 0.06024
```

<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
