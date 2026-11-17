# DDF vs SSM SGS-model comparison — building-free convective ABL

작성: 2026-05-19 (자동 분석). 케이스: CaseAS building-free, CASE_ID 213.

## 1. 셋업

| 항목 | 값 |
|---|---|
| 도메인 | CaseAS, 건물 제거 (`DEACTIVATE_IBM 1` → IBM p_type=1000 입자 i_type=3 비활성) |
| 성층 | 지표 T_con=330 K, 상단 T_ini=300 K, ABL_h=125 m → **대류 불안정** (ground +30 K) |
| MOST | ON (`MOST_boundary 1`) — 지표가열 |
| 온도 advection | Riemann-SPH (`RIEMANN_EULERT 1`) |
| 적분 | dt=1e-3, 30 s = 30000 step |
| **유일한 차이** | `turbulence-model`: **DDF** (Deardorff 1-eq) vs **SSM** (Smagorinsky) |

두 런은 step 0 초기조건 완전히 동일 — controlled 비교.

## 2. 런 건전성

- 양쪽 다 step 30000 완주, **발산 없음**: tt_max≈329.9 K, tt_min=305.2 K 전 구간 안정, NaN/Inf 0.
  → Riemann `eulert`가 30000 step 동안 견고함 재확인.
- 양쪽 다 step~3000부터 속도 15 m/s clip 도달 (모멘텀 advection은 여전히 central-diff,
  강한 대류 fbz~0.98). **두 런 공통** → 상대 비교엔 영향 작으나 절대 속도크기는 clip-제한됨.

## 3. 핵심 결과 — SGS 모델이 결과를 바꾼다

### 3.1 온도장 (fig1, fig2)
- t ≤ 10 s: DDF·SSM 거의 동일 (max|Δθ| = 0.08 K).
- 대류가 발달하며 차이 증가: **max|Δθ| = 0.08 K (10s) → 0.50 K (20s) → 0.98 K (30s)**.
- t=30 s에서 DDF가 전 고도에서 SSM보다 차가움 (지표부 −0.3 K, 상공 z~100 −1.0 K).

### 3.2 난류·열수송 (fig4, fig5) — t=30 s
| 레이어 | DDF θ̄ | SSM θ̄ | DDF ⟨w'θ'⟩ | SSM ⟨w'θ'⟩ |
|---|---|---|---|---|
| z[5,30) | 327.10 | 327.40 | 0.329 | 0.151 (+118%) |
| z[30,60) | 322.74 | 323.16 | 1.982 | 1.574 (+26%) |
| z[60,100) | 316.90 | 317.68 | 5.12 | 3.77 (+36%) |

- 컬럼평균 resolved 수직열플럭스 ⟨w'θ'⟩: DDF 3.58 vs SSM 2.97 (**+21%**).
- 수직속도 분산 ⟨w'²⟩: DDF 7.9 vs SSM 6.9.
- → **DDF가 더 활발한 대류 난류**를 만든다.

### 3.3 해석
DDF의 SGS TKE 예산에는 부력생산항 `b = -(g/θ)(ν_t/Pr_t)(dθ/dz)`가 있어, 불안정 성층
(dθ/dz<0)에서 b>0 → SGS TKE 증폭 → resolved 대류 강화. SSM(Smagorinsky)은 strain만
보고 성층에 무지. 결과적으로 DDF는 수직 열수송이 더 효율적 → 하층 열을 위로/도메인 밖으로
더 잘 환기 → 대류층이 SSM보다 차갑고 더 잘 섞임. **부력을 통해 SGS 모델이 결과를 바꾼다는
것이 정량적으로 확인됨** — 단 차이는 대류가 발달하는 t≳10 s 이후에 드러남.

## 4. 그림
- `fig1_theta_profiles.png` — θ(z) 시간변화, DDF/SSM 각 패널
- `fig2_theta_diff.png` — θ(z) 중첩 + Δθ(z)=DDF−SSM
- `fig3_sgs_tke.png` — SGS TKE·dissipation 프로파일
- `fig4_heatflux.png` — resolved ⟨w'θ'⟩(z), ⟨w'²⟩(z)
- `fig5_timeseries.png` — 컬럼평균 시계열

## 5. 주의/한계
- **속도 15 m/s clip** (양쪽): 모멘텀 advection central-diff 미해결. ⟨w'²⟩ 절대크기는
  clip 영향. 상대 비교(DDF>SSM)는 유효.
- **SSM의 k_turb·e_turb는 live 진단 아님**: Smagorinsky는 대수모델이라 TKE 방정식이
  없음. 코드가 k_turb를 DDF 분기에서만 갱신 → fig3의 SSM 곡선은 초기 MOST 프로파일에
  동결된 값. SSM 난류는 vis_t(미출력)로만 작용. fig3는 DDF 진화 + 이 방법론적 caveat 확인용.
- 온도 필터링 매 2 step → θ 평활화.
- 30 s는 대류 turnover(~60–125 s)보다 짧음 — 더 길게 돌리면 Δ가 더 커질 것.
- 코드 비대칭: SGS *열*확산항(functions_PPE.cuh)은 DDF에서만 적용 (사용자 확인, 의도적 유지).
  → 이 비교는 "DDF(SGS 열플럭스 포함) vs SSM(SGS 열플럭스 없음)" 임을 감안.

## 6. 재현
`CaseAS_compare/analyze_sgs.py` — VTK 파싱·프로파일·플롯 생성 스크립트.
