# 2026-09-14_ckpt-restart-verify

> 2026-09-14 15:23:25 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

체크포인트/재개(save_checkpoint/load_checkpoint, CKPT_EVERY, SOPHIA_RESTART) 구현 검증. Leipzig full 입력(43M)+LDM 방출을 t=0.2 s 로 당긴 CkptTest 케이스, dt=0.01, 80 step. (1) 연속 실행(count=40 에서 체크포인트 저장) vs (2) count=40 체크포인트에서 재개해 80 step 까지 — PPE_diag 행, stdout tu_max/Σk/Σ½|u|², LDM·슬라이스·전체장 VTK, mean_field.bin(평균창 0.2~0.6 s, 재개 뒤 덤프)이 비트 단위로 같으면 통과. 평균누적·LDM 방출·RNG seed·QS 상태 복원을 한 번에 검사.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `CkptTest` |
| 출력 | `runs/2026-09-14_ckpt-restart-verify/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DAVG_START_TIME=0.2` `-DAVG_END_TIME=0.6` `-DINIT_CROSSFLOW_UY=0.0` `-DCKPT_EVERY=40` `-DSLICE_DUMP_EVERY=25` `-DFULL_DUMP_EVERY=75` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-14_ckpt-restart-verify/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-14_ckpt-restart-verify/env/build.sh && bash runs/2026-09-14_ckpt-restart-verify/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-14_ckpt-restart-verify` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

**절차**: 같은 바이너리로 (1) 연속 실행 0→80 step (`run_cont.log`, `out_cont/`, count=40·80 에서 체크포인트 저장)
→ (2) `SOPHIA_RESTART=out_cont/restart.bin.prev`(count=40) 로 재개해 80 step 까지 (`run_restart.log`, `out/`).
비교는 `analysis/compare_restart.sh`.

| 비교 대상 | 결과 |
|---|---|
| PPE_diag step 40~80 (5 행) | 같음 |
| stdout step 50 진단 (tu_max / hot-u / cold-spot / Σk·Σ½\|u\|²) | 같음 |
| fluid_condition_50/75, LDM_50/75 (재개 뒤 방출·Langevin·RNG) | 비트 동일 |
| fluid_75stp.vtk 전체장 (2.24 GB, 43,046,441 입자) | 비트 동일 |
| mean_field.bin (평균창 0.2~0.6 s, 재개 뒤 t=0.6 s 덤프, avg_count=41) | 비트 동일 |
| count=80 체크포인트 (P1·P2·셀목록·LDM·장치전역·호스트상태 전체) | 비트 동일 |

→ 체크포인트가 다음 스텝에 영향을 주는 상태를 빠짐없이 담는다. 체크포인트 한 벌 10.88 GB, 쓰기 57.7~65.3 s
(fsync 포함), 재개 시 로드 수 분(입력 텍스트 읽기 포함).

첫 비교에서 PPE_diag 가 "다름"으로 나온 건 스크립트 결함이었다 — 재개 런이 새로 만든 CSV 엔 헤더가 없어
`NR>1` 필터가 첫 데이터 행을 버렸다(값은 동일). 스크립트를 고치고, 코드도 **파일이 비어 있을 때만 헤더를
쓰도록**(`ckpt_log_is_new`) 바꿨다. 이 수정은 로그 헤더 조건만 건드리므로 위 수치 검증은 그대로 유효하다
(검증 바이너리 = 수정 직전 소스).

용량 정리: 검증 뒤 restart.bin 3벌(33 GB)과 fluid_75stp.vtk 2벌(4.5 GB)은 삭제.

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**통과** — count=40 체크포인트에서 재개한 결과가 연속 실행과 비트 단위로 동일(PPE_diag·stdout·슬라이스·LDM·전체장·mean_field.bin·count=80 전체 상태). 체크포인트 10.9 GB / 쓰기 ~60 s.
