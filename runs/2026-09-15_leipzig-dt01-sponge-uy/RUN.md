# 2026-09-15_leipzig-dt01-sponge-uy

> 2026-09-15 13:42:22 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

Full run A (dt=0.01) 를 스펀지 uz+uy 감쇠로 t=120 s 부터 분기 재개. 2026-09-15_leipzig-dt01-sponge(SPONGE_GAMMA_MAX=0.2, uz 만 감쇠)에서 천장(z=340~372 m) uz std 는 0.42 -> 0.028 로 Run A 수준이 됐으나 uy std 가 0.31(t=110) -> 0.37(120) -> 0.49(130) 로 계속 성장(Run A 0.055, 고도에 따라 평균 uy -0.09~+0.13 인 횡방향 전단층), 천장 Σk(z>=250 슬래브) 425 -> 555 로 증가율 3.3 -> 7.8%/2 s 가속, 도메인 Σk 비 1.101 -> 1.155. 조치(사용자 승인 diff): functions_BC.cuh KERNEL_sponge_ceiling 에 'if(SPONGE_DAMP_UY) P1[i].uy *= damp;' 추가 + Parameters.cuh SPONGE_DAMP_UY(기본 0). 이 케이스는 INIT_CROSSFLOW_UY=0(Run A 천장 평균 uy +0.005)이라 uy->0 감쇠가 평균장을 왜곡하지 않음. 시작점: 2026-09-15_leipzig-dt01-sponge/out/restart_t120.bin(count=12000). 나머지 플래그·입력·dt 동일. 성공 판정: 천장 uy std 가 Run A 수준(~0.05)으로 감소, 천장 Σk 증가 멈춤·감소, 도메인 Σk 비 하락, 캐노피·KE·발산 비 유지, t=299 s mean_field.bin -> Run B(2026-09-15_leipzig-dt01-sponge-uy-meanflow).

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01_full` |
| 출력 | `runs/2026-09-15_leipzig-dt01-sponge-uy/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DINIT_CROSSFLOW_UY=0.0` `-DCKPT_EVERY=1000` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-15_leipzig-dt01-sponge-uy/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-15_leipzig-dt01-sponge-uy/env/build.sh && bash runs/2026-09-15_leipzig-dt01-sponge-uy/env/launch.sh`

### 운영 — 분기 재개·감시 (2026-09-15)

- **계보**: 2026-09-14_leipzig-dt01-full(스펀지 OFF, t=0~110 s) → 2026-09-15_leipzig-dt01-sponge(uz 스펀지, t=110~120 s)
  → **이 런**(uz+uy 스펀지, t=120 s~). 시각별 기록은 각 런의 PPE_diag·run.log·슬래브에 나뉘어 있다.
- **시작점**: `../2026-09-15_leipzig-dt01-sponge/out/restart_t120.bin`(count=12000) 에서 `SOPHIA_RESTART` 로 기동(13:40, GPU 4).
  step 12000 PPE_diag 행이 원런과 완전히 같음 확인. **원점 체크포인트 두 개(restart_t110.bin, restart_t120.bin)는 지우지 말 것.**
- **코드**: 사용자 승인 diff — `functions_BC.cuh` `KERNEL_sponge_ceiling` 에 uy 감쇠(플래그), `Parameters.cuh` `SPONGE_DAMP_UY`(기본 0).
- **체크포인트**: 1000 step 마다 `out/restart.bin`(+`.prev`), 첫 저장 t=130 s. 그 전에 죽으면
  `env/resume.sh runs/2026-09-15_leipzig-dt01-sponge/out/restart_t120.bin`. 정전 뒤엔 `runs/2026-09-15_leipzig-dt01-sponge-uy/env/resume.sh`.
- **감시**: `env/watchdog.sh`, `env/groundcheck.py`, `analysis/track_vs_runA.py`, `analysis/ceiling_k.py`(uy std 열 추가),
  Run B 체인 `../2026-09-15_leipzig-dt01-sponge-uy-meanflow/env/chain.sh`.

## 3. 결과

<!-- AUTO:METRICS  exp/report.sh 가 생성 — 직접 고치지 말 것 (2026-09-16 14:12:28) -->

| 측정 | 값 |
|---|---|
| 상태 | **DONE** |
| 도달 step / 시뮬레이션 시간 | 30000 / 299.98 s |
| 계산 시간 (wall) | 86717.30 s (= 0.003 sim-s/wall-s) |
| 속도 tu_max (최대 / 최종) | 20.000 / 20.000 m/s |
| uz_max (최대) | 20.00 m/s |
| 온도 범위 tt_min–tt_max | 300.00 – 300.00 K |
| 출력 | 994 vtk, 52G |
| 로그 이상 | 없음 |

**PPE_diag_explicit.csv** (PPE 잔차 진단)

```
샘플 1801개
bi_max      평균 115.9  최대 183.3  최종 124.1
resid_max   평균 72.25  최대 95.15  최종 87.21
div_un1_avg 평균 0.05482  최대 0.06322  최종 0.05897
```

<!-- /AUTO:METRICS -->

### 해석

**성공 판정 4개 모두 충족.**

| 판정 | 결과 |
|---|---|
| (1) 천장 Σk(z≥250, 슬래브) 증가 중단·감소 | 425(t=120) → 336(130) → 187(150) → 34.4(300). uy std 0.366 → 0.0026 (Run A 0.055 보다 낮음), 천장 |uz|max 0.37 → 0.011 |
| (2) 도메인 Σk 비(대 dt=0.005 Run A) 1.25 미만 복귀 | 1.079(120) → 1.031(200) → **1.004(237, Run A 기록 끝)** |
| (3) 캐노피·KE·발산 비 유지 | t=237 까지 KE 비 0.994~1.011, div_avg 비 0.96~1.29, div_plain_max 비 ≤1.47 (기준 2/3배 이내). t>237 은 자기기준(t=220~237) 대비 div ×0.97~1.06, KE ×0.995~1.000, k ×1.001~1.013 |
| (4) t=299 s mean_field.bin | `SOPHIA_MEAN_v1`, 43,046,441 레코드(유체 40,548,434 = 94.2% 비-0), avg_count=14,901, 860,928,868 B, NaN 0, ⟨ux⟩=8.375, ⟨uy⟩=+0.009, ⟨uz⟩=+0.019 → 체인이 검증 후 Run B 기동(14:10) |

**t=200 s 전체장 대 정전 전 Run A 같은 시각** (같은 방식으로 재분석): `|u|>15` 가 **1,043 vs 1,206**, 천장(z≥355) **0 vs 254** —
스펀지가 천장군을 완전히 없앴다. 2.5~10 m 는 24 vs 4 로 늘었는데 15 개가 (335~352, −442~−427) 한 패치이고 소스에서 127~149 m 떨어져 있다.
t=250 s 에는 전체 780, 그 패치도 56 → 20 으로 줄어 간헐 현상으로 확인. 평균 |u| 연직 프로파일은 전 고도대에서 Run A 와 0.08 m/s 이내.
소스 r≤15 m 이상입자는 스냅샷마다 자리를 옮긴다(t=200 소스 113 에 5, t=250 소스 108 에 7, 모두 z=1.25 m 첫 유체층이고 소스고도 2.5~10 m 에는 0).

**LDM**: t=150 s 부터 방출, 종료 시 1,299,844 개(입력 1,300,000 중), LDM VTK 600 장.
**성능**: 30,000 step, wall 86,717 s(24.1 h) = 2.89 s/step. 체크포인트 30 회(각 ~60~85 s, 합 ~35 min 포함). dt=0.005 였다면 같은 구간에 ~48 h.
**계보**: t=0~110 s 는 [2026-09-14_leipzig-dt01-full](../2026-09-14_leipzig-dt01-full/RUN.md)(스펀지 OFF), t=110~120 s 는
[2026-09-15_leipzig-dt01-sponge](../2026-09-15_leipzig-dt01-sponge/RUN.md)(uz 스펀지) 기록에 있다. 이 런의 PPE_diag·run.log 는 t=120 s 부터.

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**성공** — dt=0.01 + 천장 스펀지(uz+uy)로 t=300 s 완주(24.1 h, 2.89 s/step). 천장 불안정 소멸(z≥250 Σk 425→34, uy std 0.0026 < Run A 0.055, t=200 전체장 천장 클립 0 vs Run A 254), 정전 전 Run A 와 t=237 s 까지 KE 비 0.994~1.011·k 비 1.004·발산 기준 이내, mean_field.bin(14,901 샘플) 정상 → Run B 자동 기동. 다음: Run B(2026-09-15_leipzig-dt01-sponge-uy-meanflow) 종료 후 같은 step 의 LDM VTK 로 transient vs mean-flow 플룸 비교.
