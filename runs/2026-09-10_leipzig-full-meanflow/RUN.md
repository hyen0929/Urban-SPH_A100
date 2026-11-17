# 2026-09-10_leipzig-full-meanflow

> 2026-09-10 18:07:22 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Run B (frozen mean-flow LDM). Run A 2026-09-10_leipzig-full-bcfix_2 가 t=299 s 에 덤프하는 mean_field.bin(150~299 s 평균 ux/uy/uz/k/e, 43,046,441 rec x 20 B = 861 MB)을 첫 NNPS 정렬 직후 P1 에 덮어쓰고 SPH 블록 전체를 skip, LDM 만 푼다. 목적: transient LES 위 플룸(Run A)과 시간평균 정지장 위 플룸(Run B)의 차이 정량화 - 도시 캐노피에서 순간 난류구조가 확산에 기여하는 몫. 입력/LDM 덱/solv 는 Run A 와 동일(Leipzig_surface_s10, 13 소스 x 100,000, t_release 150~300 s 절대시각)이라 같은 step 번호끼리 1:1 대조 가능. 예상 60,000 step x 0.2 s = 3.3 h. env/chain.sh 가 Run A 종료 + mean_field.bin 완성을 감지해 자동 기동한다.

**기동 방식**: 수동 실행 금지. `env/chain.sh` (2026-09-10 18:09 기동, 백그라운드)가
Run A [2026-09-10_leipzig-full-bcfix_2](../2026-09-10_leipzig-full-bcfix_2/RUN.md) 의
프로세스 종료 → `out/mean_field.bin` 크기(≥800 MB)·magic(`SOPHIA_MEAN_v1`) 검증 →
GPU 4 반납 확인 순으로 통과해야 띄운다. 진행/실패 사유는 `env/chain.log`.
서버가 재부팅되면 chain 도 죽으므로 수동 재기동 필요:
`setsid nohup runs/2026-09-10_leipzig-full-meanflow/env/chain.sh >/dev/null 2>&1 &`

**성공 판정**
1. 기동 로그에 `[mean] loaded 43046441 records, avg_count_at_dump=29800 ...` (실패 시 `[mean] FATAL` 후 무의미한 런이 됨).
2. LDM VTK 가 Run A 와 같은 step 번호로 나옴 → 같은 step 끼리 1:1 대조 가능.
3. 플룸이 정지장 위에서 이동 — 속도장이 0 으로 덮이는 입자(경계·IBM·버퍼)에 갇히지 않는지 확인.

**분석 시 주의**: `mean_field.bin` 은 fluid(`p_type==1`, `i_type<=2`) 입자만 누적된다.
경계·IBM·버퍼 입자의 평균값은 0 이라 Run B 에서 그 입자들의 속도는 0 이다.
Run B 의 LDM `velocity` 출력은 Langevin 변동 성분만 담는다(평균풍 미저장).

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_surface_s10` |
| 출력 | `runs/2026-09-10_leipzig-full-meanflow/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DINIT_CROSSFLOW_UY=0.0` `-DLDM_MEAN_MODE=1` `-DMEAN_FIELD_PATH="./runs/2026-09-10_leipzig-full-bcfix_2/out/mean_field.bin"` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-10_leipzig-full-meanflow/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-10_leipzig-full-meanflow/env/build.sh && bash runs/2026-09-10_leipzig-full-meanflow/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-10_leipzig-full-meanflow` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**미기동·폐기** — 입력인 Run A(2026-09-10_leipzig-full-bcfix_2)가 정전으로 mean_field.bin 없이 사망. 2026-09-14_leipzig-dt01-meanflow 로 대체(dt 0.01 입력, 새 체인).
