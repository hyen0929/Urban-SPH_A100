# 2026-09-14_dt01-stability

> 2026-09-14 11:05:02 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

dt=0.01 이 2.5 m 격자에서 안정한지 판정하는 4시간 A/B 시험. 정전(9/13 09:00)으로 Run A(2026-09-10_leipzig-full-bcfix_2)가 t=237.2 s(79%, wall 62.7 h)에서 죽어 mean_field.bin 을 못 받았고, save_restart 가 ISPH_Calc.cuh:1327 에 주석으로만 있고 함수 자체가 없어 처음부터 79 h 재실행이 필요한 상황. dt 를 2배로 키우면 40 h 로 줄지만 repo 전체가 dt/dx=0.002 s/m 로 맞춰져 있고(Kori 간격 5.0 m/dt 0.01, Leipzig_surface_s10 간격 2.5 m/dt 0.005) 2.5 m 격자에서 dt=0.01 전례가 없음. 이류 CFL 은 여유(클립 20 m/s 기준 C=0.04 -> 0.08)지만 EISPH 1-sweep Jacobi 라 스텝당 발산오차는 dt 에 비례해 커지는데 흡수능력은 그대로이고, solv 의 timestep update=NO 에 ISPH_Calc.cuh:764 CFL 블록도 주석처리라 안전장치가 전혀 없음(넘쳐도 20 m/s 클립이 가려 조용히 틀림). 솔버 소스는 Run A 빌드와 완전 동일(git diff HEAD 89,068 B 일치)이고 solv 에서 dt 0.005->0.01, end-time 300->30, LDM-frequency 20->10(dt_LDM=0.1 s 유지)만 바꿈. AVG_START_TIME=150 이라 t=30 s 까지 평균누적·LDM 방출 모두 꺼짐 - Run A 의 동일 구간과 순수 유동만 비교됨. 판정기준: Run A 로그의 tu_max/uz_max(t=1.25 15.44/12.32, 5.25 13.80/11.14, 7.25 20.00/17.00, 10.0 20.00/19.6)와 PPE_diag(div_avg 7.58e-1@0 -> 1.63e-1@5.25 -> 1.26e-1@10 -> 9.18e-2@20 단조감소, div_plain_max 1.31e2@20)를 같은 물리시각에서 대조. div_avg 가 감소추세 유지 + 기준 2배 이내, div_plain_max 가 기준 3배 이내면 통과 -> 전 구간 dt=0.01 로 40 h. 실패 시 4 h 만 버리고 dt=0.005 79 h 또는 dt=0.0075 53 h 중간안. 참고: 케이스명이 22자를 넘으면 SOPHIA_gpu.cu:142 char fn[64] 가 터짐(Leipzig_surface_s10_dt01=24자로 첫 시도 실패) - 그래서 이름이 Leipzig_s10_dt01.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_s10_dt01` |
| 출력 | `runs/2026-09-14_dt01-stability/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=0` `-DSPONGE_GAMMA_MAX=0` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DINIT_CROSSFLOW_UY=0.0` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-14_dt01-stability/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-14_dt01-stability/env/build.sh && bash runs/2026-09-14_dt01-stability/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS  exp/report.sh 가 생성 — 직접 고치지 말 것 (2026-09-14 15:06:35) -->

| 측정 | 값 |
|---|---|
| 상태 | **DONE** |
| 도달 step / 시뮬레이션 시간 | 3000 / 30.00 s |
| 계산 시간 (wall) | 14292.23 s (= 0.002 sim-s/wall-s) |
| 속도 tu_max (최대 / 최종) | 20.000 / 20.000 m/s |
| uz_max (최대) | 20.00 m/s |
| 온도 범위 tt_min–tt_max | 300.00 – 300.00 K |
| 출력 | 32 vtk, 447M |
| 로그 이상 | 없음 |

**PPE_diag_explicit.csv** (PPE 잔차 진단)

```
샘플 300개
bi_max      평균 85.24  최대 182.9  최종 109.3
resid_max   평균 16.21  최대 130.8  최종 17.15
div_un1_avg 평균 0.1052  최대 0.3789  최종 0.05568
```

<!-- /AUTO:METRICS -->

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

**판정기준 보정(2026-09-14 11:20)**: PPE_diag 의 `div_un1_*`, `bi_max` 는 `rho*div(u)/dt` 로 기록된다(ISPH_Calc.cuh:499). dt 가 2배면 같은 발산이 raw 값으로 절반이 되어 1절의 "2배/3배 이내" 기준이 raw 값 기준으론 사실상 4배/6배로 느슨해진다. 따라서 비교는 `raw*dt`(= rho·div u) 로 한다 — `analysis/compare.py`. 검증: t=0 에서 두 런의 raw*dt 가 avg 3.789e-3, plain_max 1.564 로 정확히 일치(비 1.00).


**최종 대조 (t=0~29.9 s, 300점, raw*dt 기준, 위반 0점)** — `analysis/compare.py`

| t [s] | 5 | 10 | 15 | 20 | 25.5 | 26.5 | 28 | 29.9 |
|---|---|---|---|---|---|---|---|---|
| div_avg 비 (한계 2) | 1.75 | 1.51 | 1.63 | 1.45 | 1.40 | **1.94** | 1.49 | 1.40 |
| div_plain_max 비 (한계 3) | 1.81 | 0.70 | 1.27 | 1.33 | **2.98** | 2.07 | 0.66 | 0.53 |
| Σ½\|u\|² 비 | 0.991 | 0.986 | 1.002 | 1.025 | 1.014(25) | 1.010(26) | 1.004 | 0.997(29.5) |
| Σk_turb 비 | 0.95 | 0.92 | 0.89 | 0.91 | 0.94 | 0.94 | 0.94 | 0.945 |

- div_avg 비는 t=6 s 이후 1.3~1.9 에서 추세 없이 진동. 1-sweep 이라 스텝당 잔여발산 ∝ dt 로 예상한 2배보다 작다.
- div_plain_max: 시험런 raw*dt 는 t=20~30 s 동안 0.9~1.6 범위에서 유계(단조증가 아님). t=25.5 s 의 2.98 은 단발 스파이크이고 곧 1.0 수준 복귀. t≥28 s 비가 0.5 대로 떨어진 건 Run A 쪽이 1.6~2.0 으로 튀었기 때문.
- t=22.5 s 까지 보였던 운동에너지 단조증가(×1.028)는 **위상 지연**이었다: 시험런은 t≈22.5 s(1.4645e9)에 정점 후 감소해 t=29.5 s 에 ×0.997. Run A 는 t≈10 s 정점.
- Σk_turb 는 t=15 s 최저 ×0.89 에서 ×0.945 로 회복 중.
- 두 런 모두 tu_max 는 t≥7 s 부터 20 m/s 클립 포화. 시험런 종료 무렵 hot-u 표본은 전부 z=1.25 m(첫 유체층, x=295~1010) — 알려진 근지표 클립 무리.
- 주의: `h-reset frequency` 는 코드에서 사용처가 전부 주석(functions_APR*.cuh) — step 2000 은 관문이 아니었다. dt 변경으로 물리 주기가 바뀌는 스텝 기반 설정: 필터(20 step, 0.1→0.2 s), 출력 SLICE/FULL/VORT/LDM 덤프 간격 2배, 평균 샘플 수 절반. 스칼라 서브사이클은 CASE_ID>2000 전용이라 무관.
- 여유가 넉넉하진 않다: 단일 점이지만 div_avg 1.94(한계 2), div_plain_max 2.98(한계 3)까지 접근. t=30~300 s(특히 LDM 방출 150 s 이후)는 미검증.

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
**통과(위반 0/300점)** — dt=0.01 은 t=0~30 s 에서 Run A 궤적을 따라감(div_avg ×1.3~1.9, KE ×0.997, k ×0.945). 여유는 빠듯(단발 1.94/2.98). 본 런을 dt=0.01 로 가되 Run A 로그와 t=237 s 까지 계속 대조하며 이탈 시 중단 권장.
