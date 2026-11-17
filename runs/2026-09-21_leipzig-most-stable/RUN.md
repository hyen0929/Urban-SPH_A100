# 2026-09-21_leipzig-most-stable

> 2026-09-21 18:19:42 · 브랜치 `consistent-ppe` @ `50f000e` +uncommitted · GPU 4

## 1. 목적

안정 성층 ABL: 원래 풍향(회전 없음)에서 MOST 지표경계를 켜고 지표온도 T_con=285K(T_ini=300K, ΔT=-15K)로 야간 역전층을 만든다. 기준런 2026-09-15_leipzig-dt01-sponge-uy(중립, MOST off, Conduction/Boussinesq NO)와 입력 입자·LDM·스펀지·평균창(150~299s)·dt 전부 동일하고 다른 것은 (1) solv 두 줄 Conduction/Boussinesq YES, (2) MOST_boundary=1, (3) MOST 밴드 z 0~2.5m, (4) T_con=285, (5) ABL_h=1e9 뿐이다. 밴드는 필수 조정: 기본값 0.4~0.6m 는 Leipzig 최저 유체 z=1.25m 를 못 잡아 MOST 가 전부 early-return 한다(무효). ABL_h=1e9 는 MOST_boundary 가 켜질 때 functions_TIME_ISPH.cuh:53 에서 전 유체에 붙는 지균풍 강제력 f_geo=us_ini^2/ABL_h=0.002 m/s^2 을 0 으로 만들기 위한 것 — CASE_ID 1002 에서 ABL_h 는 이 항에만 쓰인다(MOST_inlet3D/turbulent_viscosity3D 는 미런치 죽은 코드). 초기 온도장은 균일 300K 이고 유입 버퍼도 300K 고정(buffer_type!=0 이라 dtemp 미적분)이라 300K 공기가 냉각면 위를 지나며 안정 내부경계층이 자라는 fetch-limited 구성이다. 목적: 안정 성층에서 연직 혼합 억제가 지표 부근 갇힘·침적·보행자 노출을 얼마나 키우는지, 그리고 중립에서 얻은 lambda_p-발자국 관계(rho_s=-0.84 A90 / +0.92 갇힘)가 성층에 대해 견고한지. 성공 판정: (1) MOST 밴드 입자 수가 0 이 아닐 것, (2) 천장/발산 건전성이 중립 기준런과 동등, (3) t=299s mean_field.bin 확보. 예상 40h.

<!-- 무엇을 확인하면 이 실험이 성공/실패인가? 비교 대상 run 이 있으면 링크. -->

## 2. 설정

| 항목 | 값 |
|---|---|
| 입력 케이스 | `Leipzig_stable` |
| 출력 | `runs/2026-09-21_leipzig-most-stable/out` |
| 빌드 플래그 | `-DCASE_ID=1002` `-DDEACTIVATE_IBM=0` `-DMOST_boundary=1` `-DL_MOST_z_lo=0.0` `-DL_MOST_z_hi=2.5` `-DT_con=285.0` `-DABL_h=1.0e9` `-DSPONGE_GAMMA_MAX=0.2` `-DSPONGE_DAMP_UY=1` `-DAVG_START_TIME=150.0` `-DAVG_END_TIME=299.0` `-DINIT_CROSSFLOW_UY=0.0` `-DCKPT_EVERY=2000` `-DSLICE_DUMP_EVERY=50` `-DFULL_DUMP_EVERY=5000` `-DVORT_DUMP_EVERY=500` |
| 실행 | `CUDA_VISIBLE_DEVICES=4 runs/2026-09-21_leipzig-most-stable/env/SOPHIA_gpu 1` |

기본값과 다른 것은 위 "빌드 플래그"가 전부다. 나머지는 `env/Parameters.cuh`,
런타임 설정은 `env/solv.txt` 참조.

**재현**: `bash runs/2026-09-21_leipzig-most-stable/env/build.sh && bash runs/2026-09-21_leipzig-most-stable/env/launch.sh`

## 3. 결과

<!-- AUTO:METRICS -->
_아직 없음 — `exp/report.sh 2026-09-21_leipzig-most-stable` 로 채운다._
<!-- /AUTO:METRICS -->


### 2026-09-22 02:37 — 국소 u* 바이너리로 교체 (step 6000 체크포인트 재개)

t=60 s 지점에서 바이너리를 바꿔 이어 돌렸다. 유동은 영향 없음(아래 회귀 검증).
바뀐 것은 LDM 건성침적·재부유가 쓰는 마찰속도·Obukhov 길이의 산정 방식뿐이다.

**왜**: MOST 기준점이 건물 안이었다. `x_ref=0.5*(L_MOST_x0+L_MOST_xmax)=995` 인데
CASE_ID 1002 의 `L_MOST_xmax=1990` 은 2000 m 급 케이스(211~213)에서 복사된 값이고
실제 도메인은 x=1315 까지다. (995, 0, 1.25) 입자는 반경 2.5 m 안에 IBM 마커 78개,
z_top=18.38 m 로 건물 내부이며 IBM 이 붙들어 |U|=0.179 m/s. 그 결과
`us_update=0.0284` (밴드 중앙값 기준 기대치 0.38), `Ri_b=8.8` 로 임계 초과라
Newton 에 해가 없어 `OB_length_global` 이 -500 클램프에 박혔다 — 지표를 식히는데
부호가 불안정로 뒤집힌 값이다. 침적속도는 `Ra`,`Rb` 모두 u* 에 반비례하므로
`Vd ~ Vs` (중력침강만), 재부유는 `u*<u*_t` 라 항상 0 이 된다.

**무엇을**: `functions_DD.cuh` 에 `DD_local_most()` 추가 — LDM 입자 주변 27 셀의
MOST 밴드 유체 입자에서 수평풍속·온도·높이를 평균낸 뒤 Newton 을 한 번 돌려
그 위치의 u*·L 을 구한다. `DD_calc_Vd`/`DD_calc_Vd_wall`/`DD_calc_P_res` 는 전역
변수 대신 인자로 받는다. 이웃이 없을 때만 전역 밴드 평균으로 폴백하며, 그 전역값은
`KERNEL_MOST_boundary3D` 의 atomicAdd 누적 + `KERNEL_MOST_band_reduce` 평균으로
만든다(1/L 을 평균 — L 은 중립에서 발산하지만 1/L 은 0 을 지나며 연속).

**검증**:
- 구조체 크기 불변 → step 4000/6000 체크포인트 그대로 로드됨
- 회귀: 새 바이너리를 step 4000 체크포인트에서 재개시켜 계산한 step 4050 진단값이
  원본 바이너리의 연속 진행 값과 **6개 항목 전부 일치**
  (`tu_max=20.00 tt_max=300.00 tt_min=285.75 fbz_max=4.661e-01 uz_max=19.62 ftotal_max=4.449e+01`)
- step 4000 체크포인트로 오프라인 계산한 전역 폴백값:
  `us_update` 0.0284 → **0.3382 m/s**, `OB_length_global` -500 → **+48.31 m** (부호 교정)
  밴드 입자의 20.2 % 가 여전히 ±500 클램프지만 1/L 평균이라 기여가 ±0.002 로 작다

**보존**: 옛 바이너리 `env/SOPHIA_gpu.pre-localust`, 소스 diff `env/localust.diff`.
재개부터는 `stdbuf -oL` 로 띄워 로그가 줄 단위로 기록된다(종전 4 KB 블록 버퍼).

**주의 — 비교 가능성**: 기존 13개 소스 침적 분석은 전부 `us_update=0.5` 균일값 위에서
나온 결과다. 이 런과 나란히 놓으려면 중립 기준런의 Run B 를 새 바이너리로 다시 돌려야
한다(평균장이 이미 있으므로 10분이면 된다). Run A 재실행은 불필요.


### 2026-09-22 — 연직속도에 2Δz 격자 진동. 기대한 하향류는 없음

**관측**: t=140 s 체크포인트, 도메인 전체 수평평균 uz 프로파일이 격자 평면마다 부호 교대.
z<40 m 16개 평면에서 13회 전환 (1.25: -0.072 / 3.75: +0.061 / 6.25: -0.052 / 8.75: +0.049 ...).
인접 평면 쌍평균으로 모드를 걸러내면 남는 값은 -0.005 ~ +0.012 m/s 로 사실상 0.
평평한 지형 위 수평 균질층에서 기대하는 바 그대로이고, **찬공기 침강은 없다.**

**대조군**: 중립 기준런(2026-09-15_leipzig-dt01-sponge-uy, 부력 off) 체크포인트에 같은 진단.

| | 중립 | 안정 |
|---|---|---|
| 부호전환 (z<40 m) | 1/15 | 13/15 |
| 2Δz 진폭 | 0.0017 m/s | 0.0330 m/s (19배) |
| 걸러낸 실제 평균 uz | +0.0152 | +0.0046 m/s |
| 진동/실제신호 | 0.1배 | 5.3배 |

그림 `analysis/f_checkerboard.png`, 스크립트 `analysis/fig_checkerboard.py`,
`compare_uz.py` / `compare_uz_avg.py` / `uz_structure.py`.

**진단 — 부력항이 원인이 아니라 방아쇠다**: `fbz_max` 는 0.475 로 전 구간 안정적이다.
docs/Temperature_advection.md §11.7 이 2026-05-18 CaseAS 검증에서 똑같은 구도를 이미 기록했다 —
`fbz_max=0.489` 로 부력은 안정인데 `uz_max` 만 6.8 -> 15.00(clip) 으로 터졌고, 원인은
**`eulerx/y/z` 모멘텀 advection 이 여전히 central-diff + KGC 라서 생기는 비단조성**으로 결론났다.
온도(`eulert`)는 Vila 1999 보존형 1차 upwind(`RIEMANN_EULERT=1`)로 교체해 해결됐지만
모멘텀은 손대지 않은 채 남아 있다. §11.8 이 남긴 다음 할 일 1번이 `RIEMANN_EULERV` 인데
**코드에 아직 없다**(현재 `RIEMANN_EULERT` 만 존재).

부력은 연직 모멘텀에 지속적 강제력을 주어 그 취약성을 여기시키는 역할이다. 중립 런에는
연직 강제력이 없어 드러나지 않았을 뿐이다. 따라서 처방 순서는 **implicit PPE 보다
`RIEMANN_EULERV` 가 먼저**다 — 같은 계열에서 이미 통한 처방이고, 온도 구현이 참조로 남아 있다.

**이 런의 처분 (사용자 결정 2026-09-22)**: 중단하지 않고 완주시켜 **explicit 기준선**으로 남긴다.
나중에 implicit / RIEMANN_EULERV 런과 짝지어 "성층에서 1-sweep explicit 이 깨지는 정도"를
정량화하는 대조군으로 쓴다. **이 런의 침적·재부유 수치는 쓰지 않는다** — z=1.25 m 에
-0.072 m/s 의 지속적 하강 편향이 있어 LDM 입자가 인위적으로 지면에 밀린다.

### 해석

<!-- 위 수치가 뜻하는 것. 1절의 질문에 답이 됐는가. 그림은 analysis/ 에. -->

## 4. 결론

<!-- 한 줄 결론 — INDEX.md 에 그대로 들어간다. 그리고 다음에 뭘 할지. -->
