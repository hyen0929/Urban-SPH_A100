# KGC & k-ε 모델 — SOPHIA_gpu 구현 수식 정리

본 문서는 현재 코드에서 실제로 사용 중인 **KGC(Kernel Gradient Correction)** 와
**표준 k-ε 난류모델**의 이산화 수식을 코드와 1:1로 정리한 것이다. 재현·검증용.

기호: 입자 `i`(대상)와 이웃 `j`, 위치 `x`, 속도 `u`, 질량 `m`, 밀도 `ρ`.
`r_ij = x_i − x_j`, `r = |r_ij|`. 커널 `W`(Wendland C6, `functions_KNL.cuh`),
비보정 기울기 `∇W_ij = (dW/dr)(r_ij / r)`.

---

## 1. KGC — Kernel Gradient Correction

구현: `KERNEL_clc_gradient_correction_3D` (`functions_KNL.cuh:118`),
동일 알고리즘 `KERNEL_clc_correction_KGC_3D` (`functions_PREP.cuh:213`),
적용 `apply_gradient_correction_3D` (`functions_KNL.cuh:105`).

### 1.1 보정행렬 B (renormalization tensor)

각 입자 `i`에 대해 3×3 대칭 텐서 `B_i` 를 이웃합으로 구성:

```
B_i = Σ_j (m_j / ρ_j) (x_j − x_i) ⊗ ∇W_ij
```

성분으로 (코드 `functions_PREP.cuh:272–289`, `txx … tzz`):

```
B_αβ = − Σ_j (m_j / ρ_j) (x_iα − x_jα) · (dW/dr)(x_iβ − x_jβ)/r
     =   Σ_j (m_j / ρ_j) (x_jα − x_iα) · ∂W_ij/∂x_β
```

- 합산 대상: 유체·경계 입자만 (`p_type ∈ {−1, 0, 1}`, 코드의 `incl` 플래그).
- B는 대칭(`B_xy = B_yx` 등)이므로 6개 성분만 계산.

### 1.2 행렬 역산 → 보정행렬 L

```
L_i = B_i^{-1}            (코드: P3[i].Cm[3][3] = inv_cm_*)
```

코드 `functions_PREP.cuh:298–312`: 행렬식 `det(B)` 를 구하고,
- `|det(B)| > Min_det (=1e-7)` → 여인수(cofactor)로 3×3 역행렬
- 그 외(이웃 부족·특이) → `L = I` (단위행렬, 보정 비활성 fallback)

### 1.3 보정된 커널 기울기

```
∇̃W_ij = L_i · ∇W_ij
```

성분 (`apply_gradient_correction_3D`, `k_kgc_solve == KGC` 일 때):

```
∂̃W/∂x = L_xx ∂W/∂x + L_xy ∂W/∂y + L_xz ∂W/∂z
∂̃W/∂y = L_xy ∂W/∂x + L_yy ∂W/∂y + L_yz ∂W/∂z
∂̃W/∂z = L_xz ∂W/∂x + L_yz ∂W/∂y + L_zz ∂W/∂z
```

(코드에서 `L` 은 `Cm[3][3]` 또는 6성분 `inv_cm_xx,xy,yy,yz,zx,zz` = `gc11,gc12,gc22,gc23,gc31,gc33` 로 저장.)

### 1.4 성질

보정 기울기는 **1차(선형) 일관성**을 보장한다. 즉 임의의 선형장 `f(x)=a+b·x` 에 대해

```
Σ_j (m_j/ρ_j)(f_j − f_i) ∇̃W_ij = ∇f = b   (정확)
```

비보정 SPH 기울기는 입자 무질서·경계 절단에서 이 성질이 깨지지만, KGC는
`B` 의 역산으로 이를 복원한다. 이후 모든 1차 미분(변형률, 온도구배, k-ε 확산
등)은 `∇W` 대신 `∇̃W` 를 사용.

---

## 2. 변형률 텐서 (k-ε 입력)

구현: `KERNEL_clc_prep3D` (`functions_PREP.cuh:140–172`). KGC 보정 기울기 사용.

대칭 변형률 텐서 `S` (유체 이웃 `p_type==1` 만 합산):

```
S_xx = Σ_j (m_j/ρ_j)(u_jx − u_ix) ∂̃W/∂x
S_yy = Σ_j (m_j/ρ_j)(u_jy − u_iy) ∂̃W/∂y
S_zz = Σ_j (m_j/ρ_j)(u_jz − u_iz) ∂̃W/∂z
S_xy = ½ Σ_j (m_j/ρ_j)[(u_jx−u_ix) ∂̃W/∂y + (u_jy−u_iy) ∂̃W/∂x]
S_yz = ½ Σ_j (m_j/ρ_j)[(u_jy−u_iy) ∂̃W/∂z + (u_jz−u_iz) ∂̃W/∂y]
S_zx = ½ Σ_j (m_j/ρ_j)[(u_jx−u_ix) ∂̃W/∂z + (u_jz−u_iz) ∂̃W/∂x]
```

변형률 크기 (`functions_PREP.cuh:170`, `P3[i].SR` 에 저장):

```
|S| = √( 2 ( S_xx² + S_yy² + S_zz² + 2S_xy² + 2S_yz² + 2S_zx² ) )
    = √( 2 S_ij S_ij )
```

---

## 3. 표준 k-ε 난류모델 (Launder–Spalding)

구현: `KERNEL_clc_ke_turb3D` (`functions_TURB.cuh:667`) → 증분 `dk/dt, dε/dt` 계산,
`KERNEL_update_turbulence` (`functions_TURB.cuh:33`) → 시간적분,
`KERNEL_turb_viscosity` (`functions_TURB.cuh:3`) → 와류점성.

상태량: SGS 난류운동에너지 `k` (`k_turb`), 소산율 `ε` (`e_turb`).

### 3.1 와류점성

```
ν_t = C_μ k² / ε              (ε > 0,  아니면 ν_t = 0)
```

운동량 방정식에는 `μ_t = ρ ν_t` 로 들어감.

### 3.2 난류 생산항 (realizability limiter 포함)

```
P = min( ν_t |S|² ,  √C_μ · k · |S| )
```

- `ν_t |S|² = C_μ (k²/ε) |S|²` — 표준 전단 생산항 (`P_k = 2ν_t S_ij S_ij = ν_t|S|²`).
- `√C_μ · k · |S|` — Kato–Launder / Durbin 전단 realizability 한계
  (정체점 등에서 생산 과대 방지). 둘 중 작은 값 채택.
- `k, ε` 는 0 division 방지를 위해 `k_safe = max(k,1e-10)`, `e_safe = max(ε,1e-10)` 로 보호.

### 3.3 난류량 확산항 (SPH, Cleary 조화평균형)

`k` 의 확산 (`functions_TURB.cuh:760–770`):

```
D_k,i = Σ_j  [ 4 m_j /(ρ_i ρ_j) ] · Γ^k_ij · (k_i − k_j) · (r_ij · ∇̃W_ij)/(r² )

  Γ^k_ij = Γ^k_i Γ^k_j / (Γ^k_i + Γ^k_j)        (조화평균)
  Γ^k    = μ + μ_t / σ_k                         (μ_t = ρ ν_t)
```

`ε` 의 확산도 동일 형태, `σ_k → σ_e`:

```
D_ε,i = Σ_j [ 4 m_j/(ρ_i ρ_j) ] · Γ^ε_ij · (ε_i − ε_j) · (r_ij·∇̃W_ij)/r²
  Γ^ε   = μ + μ_t / σ_e
```

(인자 `(r_ij·∇̃W_ij)/r²` 는 SPH 라플라시안의 Brookshaw/Cleary 이산화 핵심항.)

### 3.4 수송방정식 (증분)

```
dk/dt = P − ε + D_k
dε/dt = (ε/k)( C_ε1 P − C_ε2 ε ) + D_ε
```

코드: `dk_turb = tPi − e_turb + tmpk`, `de_turb = (e_safe/k_safe)(C_ε1·tPi − C_ε2·e_safe) + tmpe`.
**부력 생산항 없음** — 표준 전단구동 k-ε (성층 효과는 미포함; DDF 모델과의 차이점).

### 3.5 시간적분 (`KERNEL_update_turbulence`)

전진 오일러, 유체 입자(`p_type > 0`)만:

```
k^{n+1} = max( k^n + (dk/dt) Δt ,  1e-6 )
ε^{n+1} = max( ε^n + (dε/dt) Δt ,  1e-8 )
```

### 3.6 상수 (`Parameters.cuh`)

| 상수 | 기호 | 값 |
|---|---|---|
| `C_mu`   | C_μ  | 0.09 |
| `C_e1`   | C_ε1 | 1.44 |
| `C_e2`   | C_ε2 | 1.92 |
| `sigma_k`| σ_k  | 1.0  |
| `sigma_e`| σ_e  | 1.3  |
| `kappa_t`| κ    | 0.41 |

---

## 4. 호출 순서 (매 스텝, `ISPH_Calc.cuh`)

```
KERNEL_clc_gradient_correction_3D   → L_i (KGC 행렬)
KERNEL_clc_prep3D                   → S_ij, |S| (KGC 보정 기울기 사용)
KERNEL_clc_ke_turb3D                → dk/dt, dε/dt
KERNEL_update_turbulence            → k, ε 적분
KERNEL_turb_viscosity               → ν_t = C_μ k²/ε
```

---

## 5. 참고 / 한계

- KGC 역행렬이 특이하면(`|det B| < 1e-7`) 해당 입자는 `L = I` 로 폴백 → 그 입자만
  비보정 SPH 기울기 사용.
- 본 k-ε는 **전단 생산만** 포함(부력항 X). 성층 효과를 보려면 DDF(Deardorff)
  모델 사용 — `functions_PREP.cuh` 의 DDF 분기에는 부력 생산항
  `b = −(g/θ)(ν_t/Pr_t)(dθ/dz)` 가 있음.
- `functions_TURB.cuh` 의 2D·k-l-m·2D k-ε 커널은 `#if 0` 으로 비활성(현 파이프라인 미사용).
