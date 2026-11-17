# 2026-09-10_interior-leak-diag

> 2026-09-10 20:11:58 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 0

## 1. 목적

Run A 의 건물 내부 미표시 유체입자 폭주(t=6.5 s 부터 tu_max=20 클립 포화, hot-u 가 (330,-70,z=6~21) 에 고정 — 그 컬럼의 IBM 높이장 z_top=29 m 이고 2026-09-08 높이장 판정에서도 내부로 분류된 지점)가 건물 밖 유동을 오염시키는지 판정한다. 본런은 count_avg<0 동안 y in [0,7.5) 슬래브(fluid_condition)만 덤프하는데 문제 지점이 y=-70 이라 t=150 s 전에는 볼 방법이 없다. 그래서 AVG_START_TIME=6 으로 평균 게이트를 일찍 열어 전체장·vorticity 덤프를 켜고 t=10 s 에서 조기 종료(LDM_DURATION_AFTER_AVG=4). FULL/VORT_DUMP_EVERY=200 -> t=6~10 s 전체장 10 장. 판정: 건물 외피 밖 반경별 |u| 와 vorticity 가 시간에 따라 커지면 오염 확산 -> 내부유체 p_type=99 처리를 끝내고 본런 재기동해야 하고, 건물 내부에 갇혀 있으면 본런 속행. 입력/플래그는 Run A(2026-09-10_leipzig-full-bcfix_2)와 동일.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_surface_s10` |
| 출력 | `runs/2026-09-10_interior-leak-diag/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DINIT_CROSSFLOW_UY=0.0` `-DAVG_START_TIME=6.0` `-DAVG_END_TIME=10.0` `-DLDM_DURATION_AFTER_AVG=4.0` `-DFULL_DUMP_EVERY=200` `-DVORT_DUMP_EVERY=200` |
| 실행 | `CUDA_VISIBLE_DEVICES=0 runs/2026-09-10_interior-leak-diag/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-10_interior-leak-diag/env/build.sh && bash runs/2026-09-10_interior-leak-diag/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-10_interior-leak-diag` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
