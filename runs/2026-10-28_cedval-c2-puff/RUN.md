# 2026-10-28_cedval-c2-puff

> 2026-10-28 12:06:10 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 2

## 1. 목적

확산 검증: CEDVAL Complexity 2(BL3-2-1 엇갈림 큐브 배열) 지면 퍼프 60개를 풍동 퍼프 163회 앙상블과 비교 — 평균 dosage·도착/최대농도 분포. 1편 확산 검증 공백 채움

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CEDVAL_C2` |
| 출력 | `runs/2026-10-28_cedval-c2-puff/out` |
| 빌드 플래그 | `-DCASE_ID=1020` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=2.5` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_Z_BOTTOM=150.0` `-DSPONGE_DAMP_UY=1` `-DINIT_CROSSFLOW_UY=0.0` `-DDD_DISABLE=1` `-DABL_ALPHA=0.27` `-DAVG_START_TIME=200.0` `-DAVG_END_TIME=999.0` `-DCKPT_EVERY=5000` `-DLDM_DUMP_EVERY=100` `-DSLICE_DUMP_EVERY=500` `-DFULL_DUMP_EVERY=50000` `-DVORT_DUMP_EVERY=10000` |
| 실행 | `CUDA_VISIBLE_DEVICES=2 runs/2026-10-28_cedval-c2-puff/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-10-28_cedval-c2-puff/env/build.sh && bash runs/2026-10-28_cedval-c2-puff/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-10-28_cedval-c2-puff` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
