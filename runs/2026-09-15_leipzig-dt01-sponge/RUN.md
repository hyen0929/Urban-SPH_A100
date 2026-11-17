# 2026-09-15_leipzig-dt01-sponge

> 2026-09-15 10:48:27 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Full run A (dt=0.01) 를 천장 스펀지 ON 으로 t=110 s 부터 분기 재개. 2026-09-14_leipzig-dt01-full(SPONGE_GAMMA_MAX=0)에서 천장(z>=250 m) k 가 10 s 마다 ~3배 기하급수 성장(Σk z>=250: t=90 36 -> 110 314 -> 130 1752, 천장 |uz|max 1.3~1.6 m/s)해 도메인 Σk 비가 t=128.5 s 에 1.254(WARN). dt=0.005 Run A 에도 같은 불안정이 있었으나 느렸음(t=236 에 125, ~18 s 배증) — 원인은 반사성 freeslip 천장 아래 격자규모 진동(9/12 천장군 클립 473개와 동일), dt=0.01 은 ~130 s 앞당김. 조치: 이전 런의 count=11000 체크포인트(out/restart_t110.bin, Σk z>=250 = 314)에서 SOPHIA_RESTART 로 재개, 빌드 플래그만 -D SPONGE_GAMMA_MAX=0.2 (Z_BOTTOM 90 m, ramp^2, uz 만 감쇠, 천장 e-folding 10 s). 나머지 플래그·입력·dt 동일. 성공 판정: (1) 천장 Σk(z>=250, 슬래브)가 t=130 s 까지 증가를 멈추고 감소, (2) 도메인 Σk 비(대 dt=0.005 Run A)가 1.25 아래로 복귀, (3) 캐노피 지표·KE·발산 비 유지, (4) t=299 s mean_field.bin -> Run B(2026-09-15_leipzig-dt01-sponge-meanflow). 주의: 천장 성장 e-folding(~9 s)이 스펀지 최대 감쇠(10 s)와 비슷해 0.2 로 못 누를 수도 있음 — 그 경우 보고.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-15_leipzig-dt01-sponge/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DINIT_CROSSFLOW_UY=0.0` `-DCKPT_EVERY=1000` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-15_leipzig-dt01-sponge/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-15_leipzig-dt01-sponge/env/build.sh && bash runs/2026-09-15_leipzig-dt01-sponge/env/launch.sh`

### 운영 — 분기 재개·감시 (2026-09-15)

- **시작점**: `../2026-09-14_leipzig-dt01-full/out/restart_t110.bin`(count=11000, t=110.0051) 에서 `SOPHIA_RESTART` 로 기동(10:47, GPU 4).
  step 11000 PPE_diag 행이 원런과 완전히 같음 확인. **원런의 restart_t110.bin 은 이 런의 원점이니 지우지 말 것.**
- **체크포인트**: 1000 step 마다 `out/restart.bin`(+`.prev`). 첫 저장은 t=120 s. 그 전에 죽으면
  `env/resume.sh ../2026-09-14_leipzig-dt01-full/out/restart_t110.bin` 처럼 경로를 주고 재개(경로는 리포 루트 기준).
- **정전·재부팅 뒤**: `runs/2026-09-15_leipzig-dt01-sponge/env/resume.sh`.
- **감시**: `env/watchdog.sh`(5 min), `env/groundcheck.py`(30 min), `analysis/track_vs_runA.py`(대 dt=0.005 Run A, step*dt 정렬),
  `analysis/ceiling_k.py`(슬래브 z≥250 Σk·천장 |uz|max → `analysis/ceiling_k.csv`), Run B 체인 `../2026-09-15_leipzig-dt01-sponge-meanflow/env/chain.sh`.
- 이 런의 `PPE_diag`·`run.log` 는 t=110 s 부터만 있다. t<110 s 는 원런 파일을 볼 것.

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-15_leipzig-dt01-sponge` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

### 결과 — uz 는 눌렀으나 uy 전단 모드가 남음 (2026-09-15 13:30)

| z 340~372 m 슬래브 | std ux | std uy | std uz | k 평균 |
|---|---|---|---|---|
| Run A (dt 0.005) t=130 | 0.13 | 0.055 | 0.034 | 1.5e-4 |
| 스펀지 OFF t=130 | 0.22 | 0.45 | 0.35 | 3.6e-2 |
| 이 런 t=110 (분기점) | 0.21 | 0.31 | 0.42 | 1.5e-2 |
| 이 런 t=120 | 0.13 | 0.37 | 0.10 | 1.8e-2 |
| 이 런 t=130 | 0.11 | 0.49 | 0.028 | 2.1e-2 |

천장 Σk(z≥250 슬래브, `analysis/ceiling_k.csv`): 314(110) → 425(120) → 555(130), 2 s 증가율 3.3→7.8 % 로 가속.
uy 는 체커보드가 아니라(인접차/std 0.84~1.03) 고도에 따라 평균 uy 가 −0.09~+0.13 으로 바뀌는 횡방향 전단층.
120 m 이하·KE(0.997)·발산 비 정상. 사용자 승인 diff 로 `KERNEL_sponge_ceiling` 에 `if(SPONGE_DAMP_UY) P1[i].uy *= damp;`
+ `Parameters.cuh` `SPONGE_DAMP_UY`(기본 0) 추가, t=120 s 체크포인트(`out/restart_t120.bin` 하드링크 보존)에서
[2026-09-15_leipzig-dt01-sponge-uy](../2026-09-15_leipzig-dt01-sponge-uy/RUN.md) 로 분기. 분기 런 step 12000 PPE_diag 행이 이 런과 동일.
이 런은 t=131.2 s(step 13120)에서 kill. `out/restart.bin`(t=130)·`.prev`(=t=120, restart_t120.bin 과 같은 inode)는 원점 외엔 불필요.

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**중단(t=131.2 s)** — uz 전용 스펀지(0.2)가 천장 uz 는 Run A 수준으로 눌렀으나(std 0.42→0.028, |uz|max 1.39→0.12) uy 전단 모드가 계속 성장(std 0.31→0.37→0.49, Run A 0.055)해 천장 Σk 425→555(가속), 도메인 Σk 비 1.101→1.155. t=120 s 체크포인트(out/restart_t120.bin, 분기 원점이라 보존)에서 uz+uy 감쇠(SPONGE_DAMP_UY=1, 사용자 승인 코드수정)로 분기: 2026-09-15_leipzig-dt01-sponge-uy.
