# 표준 k-ε 난류 모델 구현 정리

SOPHIA_gpu(GPU ISPH 솔버)에 구현된 **표준 k-ε 난류 모델(Launder-Spalding)** 의
실제 구현 내용과 수식을 정리한 문서입니다. 이 문서 하나만 보고도 입력 → 초기화 →
매 시간스텝 계산 → 운동량 방정식 환류(feedback)까지 전체 흐름이 파악되도록 구성했습니다.

> 검증 목적: LES 계열 SGS 모델(DDF/Deardorff one-equation)과 비교하기 위한
> RANS 기준 모델로 사용됨. (`ISPH_Calc.cuh:159` 주석 참고)

---

## 1. 모델 선택과 활성화 경로

| 단계 | 위치 | 내용 |
|------|------|------|
| 런타임 입력 | `input/<CASE>/<CASE>_solv.txt` | 난류 모델 문자열 지정 |
| 문자열 파싱 | `function_init.cuh:279` | `"k-e"` → `turbulence_model = 2` |
| 상수 매핑 | `Parameters.cuh:47` | `#define K_E 2` |
| GPU 상수 | `function_init.cuh:104` | `k_turbulence_model = k_vii[18]` (`__constant__`) |
| 커널 분기 | `ISPH_Calc.cuh:160` | `if(turbulence_model==K_E){ ... }` |

다른 난류 모델과의 관계 (`Parameters.cuh:45-49`):

```
Laminar = 0   k-lm = 1   K_E = 2   HB = 3   SSM = 4   DDF = 5
```

k-ε는 **2-방정식 모델**로, 난류 운동에너지 `k`와 그 소산율 `ε`를 각각의 수송
방정식으로 풀어 와류 점성 `ν_t`를 결정합니다. (SSM/DDF는 SGS 모델로
`functions_PREP.cuh` 안에서 직접 처리되는 반면, k-ε는 전용 커널을 거칩니다.)

---

## 2. 모델 상수

`Parameters.cuh`에 정의된 표준 k-ε 계수입니다.

| 상수 | 값 | 의미 | 위치 |
|------|----|----|------|
| `C_mu`    | 0.09 | 와류 점성 계수 | `Parameters.cuh:147` |
| `C_e1`    | 1.44 | ε 생성항 계수 | `Parameters.cuh:148` |
| `C_e2`    | 1.92 | ε 소산항 계수 | `Parameters.cuh:149` |
| `sigma_k` | 1.0  | k 방정식 난류 Prandtl 수 | `Parameters.cuh:150` |
| `sigma_e` | 1.3  | ε 방정식 난류 Prandtl 수 | `Parameters.cuh:151` |
| `kappa_t` | 0.41 | von Kármán 상수 | `Parameters.cuh:152` |
| `Pr_t`    | 0.71 | 난류 Prandtl 수(스칼라 확산) | `Parameters.cuh:118` |

---

## 3. 지배 방정식

표준 k-ε 모델은 다음 두 개의 수송 방정식으로 구성됩니다.

**난류 운동에너지 k:**

```
Dk/Dt = P_k − ε + (1/ρ) ∇·[ (μ + μ_t/σ_k) ∇k ]
```

**소산율 ε:**

```
Dε/Dt = (ε/k)·( C_ε1·P_k − C_ε2·ε ) + (1/ρ) ∇·[ (μ + μ_t/σ_e) ∇ε ]
```

**와류 점성(닫힘 관계):**

```
ν_t = C_μ · k² / ε
```

여기서 `P_k`는 난류 생성항(production), `μ = ρν`는 분자 동점성, `μ_t = ρν_t`는
와류 점성입니다. 우변 첫째항(생성·소산)은 국소(local) 항이고, 셋째항(확산)만
이웃 입자 합산이 필요합니다.

---

## 4. 전체 워크플로 내 위치

매 시간스텝에서 k-ε는 **전처리(strain rate 계산) 직후, 경계조건/예측자 이전**에
호출됩니다. (`ISPH_Calc.cuh`)

```
시간스텝 t:
  1. NNPS               이웃 탐색 / 셀 정렬
  2. gradient_correction KGC 보정행렬
  3. KERNEL_clc_prep3D   ← 변형률 |S| 계산, P3[i].SR 에 저장        (PREP)
  ┌─────────────────────────────────────────────────────────────┐
  │ 4. k-ε 블록  if(turbulence_model==K_E)        ISPH_Calc.cuh:160│
  │   (a) KERNEL_clc_ke_turb3D   → dk_turb, de_turb 계산           │
  │   (b) KERNEL_update_turbulence → k, ε 시간적분 + 하한 클램프    │
  │   (c) KERNEL_turb_viscosity   → ν_t = C_μ k²/ε                 │
  └─────────────────────────────────────────────────────────────┘
  5. 경계조건 (no-slip, MOST)
  6. KERNEL_advection_force3D  ← ν_t 가 점성력에 사용됨            (PPE)
  7. PPE → 압력력 → corrector
  8. Diffusion (스칼라/에너지 — ν_t/Pr_t 사용)
  ...
```

핵심: **3번에서 계산된 변형률 `P3[i].SR`** 가 k-ε 생성항의 입력이 되고,
**4번에서 갱신된 `ν_t`** 가 6·8번의 점성력/스칼라 확산으로 되먹임됩니다.

---

## 5. 초기 조건

`functions_NNPS.cuh`의 입자 초기화 단계에서 `k`, `ε`의 연직 분포를 부여합니다.
중립 대기경계층(MOST) 케이스(`CASE_ID` 211~219, 3001)의 경우:

```cpp
HP1[i].k_turb = Interpolate_k_MOST_Neutral(HP1[i].z);   // functions_NNPS.cuh:39,57
HP1[i].e_turb = Interpolate_e_MOST_Neutral(HP1[i].z);
```

프로파일 정의 (`functions_initial_conditions.cuh:284-298`):

```
k(z)  = (u*²/√C_μ)·(1 − z/h)²  + 1e-4        z < ABL_h,  아니면 1e-4
ε(z)  = u*³ / (κ · z)                         z 는 [z0, ABL_h] 로 클램프
```

- `k(z)`: 지표에서 `u*²/√C_μ ≈ u*²/0.3`, 경계층 상단으로 갈수록 0으로 감쇠.
- `ε(z)`: 중립 표면층의 평형 가정 `ε = u*³/(κz)` — 로그층 해석해와 일치.

채널 검증(`CASE_ID==1500`)은 `k=ε=0`에서 시작해 브로드밴드 난수 속도장으로
천이를 유도합니다.

---

## 6. (a) 수송항 커널 — `KERNEL_clc_ke_turb3D`

`functions_TURB.cuh:667-792`. 입자별로 `dk_turb`, `de_turb`(우변 RHS)를 계산하여
`P3`에 저장합니다. 시간적분은 하지 않습니다(=(b) 단계가 담당).

### 6.1 입력 읽기

```cpp
tSR    = P3[i].SR;            // 변형률 크기 |S|, PREP에서 계산됨
vis_ti = P1[i].vis_t * rhoi;  // 동점성 μ_t = ρ·ν_t
k_turbi = P1[i].k_turb;
e_turbi = P1[i].e_turb;
visi   = viscosity(tempi,ptypei);   // 분자 동점성 μ
```

### 6.2 확산항 — SPH 이산화 (이웃 입자 합산)

3×3×3 셀을 순회하며 탐색반경 `κ·h` 안의 이웃 `j`에 대해 합산합니다.
표준 Cleary/Monaghan 형식의 SPH 라플라시안 연산자를 사용합니다.

```
(1/ρ)∇·(Γ∇k)|_i  ≈  Σ_j  4·(m_j/(ρ_i ρ_j)) · Γ̄_ij · (k_i − k_j) · (r_ij·∇W_ij)/|r_ij|²
```

여기서 유효 확산계수와 입자간 평균은:

```
Γ_k = μ + μ_t/σ_k          (k 방정식)
Γ_e = μ + μ_t/σ_e          (ε 방정식)
Γ̄_ij = (Γ_i·Γ_j)/(Γ_i+Γ_j)   ← 조화평균(harmonic mean)
```

코드 대응 (`functions_TURB.cuh:760-771`):

```cpp
tmprr = 1.0/(tdist*tdist+1e-7) * ((xi-xj)*tdwx + (yi-yj)*tdwy + (zi-zj)*tdwz);

tmpkk  = 4*mj/(rhoi*rhoj);
tmpkk *= (visi+vis_ti/sigma_k)*(visj+vis_tj/sigma_k)
         /(visi+vis_ti/sigma_k + visj+vis_tj/sigma_k);   // 조화평균 Γ̄
tmpkk *= (k_turbi - k_turbj);

tmpee  = 4*mj/(rhoi*rhoj);
tmpee *= (visi+vis_ti/sigma_e)*(visj+vis_tj/sigma_e)
         /(visi+vis_ti/sigma_e + visj+vis_tj/sigma_e);
tmpee *= (e_turbi - e_turbj);

tmpk += tmpkk*tmprr;       // → k 확산항 누적
tmpe += tmpee*tmprr;       // → ε 확산항 누적
```

조화평균을 쓰는 이유: 물성이 급변하는 계면(벽 근처, 유체/IBM 경계)에서도
플럭스 연속성을 유지하고 수치 안정성을 확보하기 위함입니다.

### 6.3 생성항 — 실현가능성 제한자(realizability limiter)

표준 생성항은 `P_k = ν_t·|S|²`이며 `ν_t = C_μ k²/ε`를 대입하면:

```
P_k = C_μ · k² / ε · |S|²
```

여기에 **Kato-Launder / Durbin 전단 실현가능성 제한자**를 적용합니다
(`functions_TURB.cuh:779-784`):

```cpp
Real e_safe = fmaxf(e_turbi, 1e-10f);   // ε 하한 — 0으로 나눔 방지(초기조건)
Real k_safe = fmaxf(k_turbi, 1e-10f);   // k 하한

tPi  = C_mu*k_safe*k_safe/e_safe * tSR*tSR;   // 표준 생성항  C_μ k²/ε ·|S|²
tPi2 = sqrt(C_mu)*k_safe*tSR;                 // 실현가능 상한  √C_μ·k·|S|
tPi  = min(tPi, tPi2);                        // P_k = min(둘 중 작은 값)
```

- `tPi`: 표준 k-ε 생성항.
- `tPi2`: 와류 점성이 무한정 커지지 못하도록 하는 상한.
  정체점(stagnation point)이나 큰 변형률 영역에서 `P_k`가 비물리적으로
  폭주하는 것을 막습니다(건물 모서리 주변에서 특히 중요).
- 둘 중 작은 값을 생성항으로 채택.

### 6.4 우변(RHS) 조립

```cpp
dk_turbi = tPi − e_turbi + tmpk;                              // Dk/Dt
de_turbi = e_safe/k_safe*(C_e1*tPi − C_e2*e_safe) + tmpe;     // Dε/Dt

P3[i].dk_turb = dk_turbi;
P3[i].de_turb = de_turbi;
```

수식 대응:

```
dk/dt = P_k − ε + Diff_k
dε/dt = (ε/k)·(C_ε1·P_k − C_ε2·ε) + Diff_ε
```

---

## 7. (b) 시간적분 — `KERNEL_update_turbulence`

`functions_TURB.cuh:33-46`. 명시적 오일러(explicit Euler) 적분과 하한 클램프.

```cpp
P1[i].k_turb += P3[i].dk_turb * tdt * (p_type_i>0);   // 유체 입자만(p_type>0)
P1[i].e_turb += P3[i].de_turb * tdt * (p_type_i>0);

P1[i].k_turb = fmaxf(P1[i].k_turb, 1e-6f);   // k > 0 보장
P1[i].e_turb = fmaxf(P1[i].e_turb, 1e-8f);   // ε > 0 보장
```

- `(p_type_i>0)`: 경계 입자(`p_type==0`)는 갱신하지 않음 — k, ε 고정.
- 하한 클램프: 음수 `k`/`ε`로 인한 `ν_t` 발산·NaN을 방지.
  주석대로 DDF 모델의 클램프 범위와 일관성을 맞춤.

---

## 8. (c) 와류 점성 갱신 — `KERNEL_turb_viscosity`

`functions_TURB.cuh:3-12`. 닫힘 관계식 적용.

```cpp
if(P1[i].e_turb > 0) P1[i].vis_t = C_mu * P1[i].k_turb*P1[i].k_turb / P1[i].e_turb;
else                 P1[i].vis_t = 0;
```

```
ν_t = C_μ · k² / ε
```

> `vis_t`는 **운동학적(kinematic) 와류 점성 ν_t** 로 저장됩니다.
> 운동량/확산 방정식에서 사용할 때 `ρ`를 곱해 동점성 `μ_t`로 변환합니다.
> (`functions_TURB.cuh:9` 주석 참고)

---

## 9. 운동량/스칼라 방정식으로의 되먹임

갱신된 `ν_t`는 이후 단계에서 다음과 같이 소비됩니다.

### 9.1 점성력 (운동량 방정식) — `functions_PPE.cuh`

`KERNEL_advection_force3D`에서 유효 점성을 분자+와류 점성의 합으로 사용:

```cpp
visi = viscosity2(tempi,ptypei) + P1[i].vis_t*rhoi;   // μ_eff = μ + μ_t   (PPE:262)
visj = viscosity2(tempj,ptypej) + P1[j].vis_t*rhoj;   //                   (PPE:375,395)
```

점성력은 Cleary 또는 Monaghan 형식으로 이산화되며, `μ_eff`의 조화평균을
입자쌍 점성으로 사용합니다(`functions_PPE.cuh:373-389`).

### 9.2 스칼라/열 확산 — `functions_PPE.cuh:451`

```cpp
tmp_Rc_diff += mj/rhoj * ((P1[i].vis_t+P1[j].vis_t)/Pr_t) * (tempi-tempj) * tdwij/tdist;
```

```
난류 열확산도 = ν_t / Pr_t        (Pr_t = 0.71)
```

즉 k-ε가 만든 `ν_t`는 운동량(점성력)뿐 아니라 온도/농도 스칼라 수송의
난류 확산까지 일관되게 구동합니다.

---

## 10. 변형률 |S| 의 출처 — `KERNEL_clc_prep3D`

생성항 입력인 `P3[i].SR`은 k-ε 커널 이전 단계인 PREP에서 계산됩니다
(`functions_PREP.cuh:166-172`).

```cpp
tmp_SR = 2*( st11² + st22² + st33² + 2·st12² + 2·st31² + 2·st23² );
tmp_SR = sqrt(tmp_SR);
P3[i].SR = tmp_SR;   // k-ε 커널에 |S| 노출
```

```
|S| = √( 2·S_ij·S_ij )
```

`S_ij`는 변형률 텐서로, SPH 속도 구배(KGC 보정 포함)로부터 계산됩니다.
이 값은 `k_turbulence_model != Laminar`일 때만 계산됩니다.

---

## 11. 데이터 흐름 요약 (한 시간스텝)

```
                    ┌──────────────┐
  PREP ───────────► │  P3[i].SR    │ 변형률 |S|
                    └──────┬───────┘
                           │
                           ▼
   ┌────────────────────────────────────────────┐
   │ KERNEL_clc_ke_turb3D                        │
   │   생성항  P_k = min(C_μk²/ε·|S|², √C_μ·k|S|) │
   │   확산항  Σ_j 4m_j/(ρρ)·Γ̄·(Δk)·(r·∇W)/r²    │
   │   → P3[i].dk_turb,  P3[i].de_turb           │
   └────────────────────┬───────────────────────┘
                         ▼
   ┌────────────────────────────────────────────┐
   │ KERNEL_update_turbulence                    │
   │   k += dk·dt,  ε += dε·dt   (오일러 적분)    │
   │   k ≥ 1e-6,  ε ≥ 1e-8       (하한 클램프)    │
   └────────────────────┬───────────────────────┘
                         ▼
   ┌────────────────────────────────────────────┐
   │ KERNEL_turb_viscosity                       │
   │   ν_t = C_μ · k² / ε   →  P1[i].vis_t        │
   └────────────────────┬───────────────────────┘
                         ▼
        운동량 점성력 (μ + μ_t)  +  스칼라 확산 (ν_t/Pr_t)
```

---

## 12. 안정화 장치 요약 (최근 보완 항목)

GPU 명시적 적분에서 k-ε는 강성(stiff)이라 폭주하기 쉽습니다. 이를 막기 위한
장치들:

| 장치 | 위치 | 목적 |
|------|------|------|
| 실현가능성 제한자 `min(tPi,tPi2)` | `functions_TURB.cuh:782-784` | 큰 변형률(건물 모서리)에서 P_k 폭주 차단 |
| `e_safe`, `k_safe` 하한 | `functions_TURB.cuh:780-781` | 초기조건 `ε≈0`에서 0 나눗셈 방지 |
| `tmprr`의 `+1e-7` | `functions_TURB.cuh:760` | 입자 근접 시 분모 보호 |
| 적분 후 하한 클램프 | `functions_TURB.cuh:44-45` | 음수 k/ε → NaN 방지 |
| 경계입자 갱신 제외 `(p_type>0)` | `functions_TURB.cuh:41-42` | 벽면 k/ε 고정 |
| MOST 평형 초기 프로파일 | `functions_initial_conditions.cuh:284-298` | 물리적으로 일관된 출발점 |

---

## 13. 관련 커널·파일 빠른 참조

| 항목 | 위치 |
|------|------|
| k-ε 분기 호출 | `ISPH_Calc.cuh:159-167` |
| 수송항 RHS 커널 | `KERNEL_clc_ke_turb3D` — `functions_TURB.cuh:667` |
| 시간적분 커널 | `KERNEL_update_turbulence` — `functions_TURB.cuh:33` |
| 와류 점성 커널 | `KERNEL_turb_viscosity` — `functions_TURB.cuh:3` |
| 변형률 |S| 계산 | `KERNEL_clc_prep3D` — `functions_PREP.cuh:166` |
| 초기 k/ε 프로파일 | `functions_initial_conditions.cuh:284` |
| 모델 상수 | `Parameters.cuh:147-152` |
| 입력 파싱(`"k-e"`) | `function_init.cuh:279` |
| 점성력에서 ν_t 사용 | `functions_PPE.cuh:262, 375, 451` |

> 참고: `functions_TURB.cuh`의 2D 버전·k-lm·SPS 스트레스 텐서 커널들은
> 현재 구조체 필드(`P2.SR`, `P3.Sxx` 등) 부재로 `#if 0` 처리되어 비활성
> 상태입니다. 3D k-ε 파이프라인(위 3개 커널)만 컴파일·실행됩니다.
