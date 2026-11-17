# 2026-09-14_leipzig-dt01-full

> 2026-09-14 16:03:45 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Full run A 재실행 (transient LES + LDM, t=0~300 s). 이전 Run A(2026-09-10_leipzig-full-bcfix_2)가 9/13 정전으로 t=237.2 s 에서 죽어 mean_field.bin 미획득. 바뀐 점: dt 0.005->0.01(시험 2026-09-14_dt01-stability 통과, LDM-frequency 10 으로 dt_LDM=0.1 s 유지), 체크포인트 CKPT_EVERY=1000(10 sim-s, ~80 min wall; 검증 2026-09-14_ckpt-restart-verify 비트동일), 출력 간격은 스텝 기반이라 SLICE/FULL/VORT_DUMP_EVERY 를 절반으로 해 이전 Run A 와 같은 물리 간격(0.5/50/5 s). 나머지 플래그·입력은 이전 Run A 와 동일. 성공 판정: (1) t=299 s 에 mean_field.bin 덤프 -> Run B(2026-09-14_leipzig-dt01-meanflow) 체인 기동, (2) t<=237 s 에서 이전 Run A 대비 이탈 없음(analysis/track_vs_runA.py: div_avg 비 2 s 중앙값<=2, div_plain_max 비<=3, KE 비 +-10%). 이탈 FAIL 시 중단하고 dt=0.005 로 전환. 예상 ~40 h + 체크포인트 30회x1 min.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-14_leipzig-dt01-full/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DINIT_CROSSFLOW_UY=0.0` `-DCKPT_EVERY=1000` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-14_leipzig-dt01-full/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-14_leipzig-dt01-full/env/build.sh && bash runs/2026-09-14_leipzig-dt01-full/env/launch.sh`

### 운영 — 재개·감시 (2026-09-14)

- **체크포인트**: 1000 step(10 sim-s, ~80 min)마다 `out/restart.bin`(10.9 GB, 직전 한 벌은 `.prev`), 쓰기 ~1 min.
- **정전·재부팅·비정상 종료 뒤**: `runs/2026-09-14_leipzig-dt01-full/env/resume.sh` 하나만 실행.
  체크포인트 끝 magic 확인(깨졌으면 `.prev`) → GPU 4 비었는지 확인 → `SOPHIA_RESTART=<체크포인트>` 로 재기동,
  `run.log` 에 이어 씀 → watchdog·groundcheck·Run B 체인이 죽어 있으면 다시 띄움. 로드까지 수 분.
  재개 뒤 `PPE_diag`·`run.log` 에는 마지막 체크포인트~사망 사이 구간이 두 번 기록된다(분석 스크립트는 시각 키로 덮어씀).
- **감시** (전부 `setsid` 백그라운드, 결과는 `analysis/watchdog.log`):
  - `env/watchdog.sh` 5 min — NaN, 정체, 로그 에러·체크포인트 실패, 비정상 종료(재개 안내), 이탈판정 변화, 디스크, Run B 상태
  - `env/groundcheck.py` 30 min — 슬래브 지표층 속도 분포(`analysis/ground_layer.csv`)
  - `analysis/track_vs_runA.py` — 이전 Run A(dt=0.005) 대비 div_avg·div_plain_max(raw×dt, 2 s 중앙값)·Σ½|u|²·Σk_turb 비, t≤237 s
  - `../2026-09-14_leipzig-dt01-meanflow/env/chain.sh` — Run A 정상 종료 + mean_field.bin 검증 후 Run B 기동 (재개 추적)
- stdout 은 `stdbuf -oL` 로 줄 버퍼링(`env/launch.sh`).

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-14_leipzig-dt01-full` 로 채운다._
<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

### 천장 불안정 — 중단 (2026-09-15 10:45)

t≈100 s 부터 도메인 Σk_turb 비(대 dt=0.005 Run A)가 1.03 → 1.079(110) → 1.143(120) → 1.254(128.5, watchdog WARN)로 가속 증가.
슬래브 고도별 분해(t=115~120 s)에서 지표~120 m 는 Run A 와 같거나 낮고(비 0.66~1.00), 초과분은 z≥250 m 에 몰려 있었다.

| t [s] | z≥250 Σk: Run A(dt 0.005) | 이 런(dt 0.01) | 천장 z≥355 \|uz\|max: A | 이 런 |
|---|---|---|---|---|
| 50 | 7.3 | 7.8 | 0.02 | 0.09 |
| 90 | 7.7 | 35.9 | 0.06 | 1.17 |
| 110 | 8.2 | 314 | 0.10 | 1.40 |
| 130 | 9.2 | 1752 | 0.14 | 1.29 |
| 236 | 124.7 | — | 0.20 | — |

Run A 에도 같은 기하급수 성장이 있었다(t=200→236 에 ~18 s 배증; 9/12 전체장의 천장군 클립 473개). 원인은
`SPONGE_GAMMA_MAX=0` 인 반사성 freeslip 천장 아래 격자규모 진동이고, dt=0.01 은 이를 ~130 s 앞당겨 평균창(150 s~) 이전으로 끌어왔다.
캐노피·KE(비 0.998)·발산 비는 정상이었다. 사용자 결정으로 t=110 s 체크포인트(`restart.bin.prev` 를 하드링크로 `out/restart_t110.bin`
보존)에서 `-D SPONGE_GAMMA_MAX=0.2` 만 바꿔 분기 → [2026-09-15_leipzig-dt01-sponge](../2026-09-15_leipzig-dt01-sponge/RUN.md).
분기 런의 step 11000 PPE_diag 행이 이 런과 완전히 같음을 확인(진단은 스펀지 적용 전에 계산됨 → 상태 복원 확인).
이 런은 t=139.5 s(step 13950)에서 kill. `out/restart.bin`(t=130)·`.prev`(t=120) 는 필요 없으면 지워도 됨(22 GB).

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**중단(t=139.5 s)** — 천장(z≥250 m) k 가 10 s 마다 ~3배 기하급수 성장(Σk z≥250: 36@90 → 314@110 → 1752@130, 천장 |uz|max 1.3~1.6 m/s)해 도메인 Σk 비 1.254(WARN). dt=0.005 Run A 에도 있던 스펀지 OFF 천장 불안정을 dt=0.01 이 ~130 s 앞당김. t=110 s 체크포인트(out/restart_t110.bin, 분기 원점이라 보존)에서 SPONGE_GAMMA_MAX=0.2 로 분기: 2026-09-15_leipzig-dt01-sponge.
