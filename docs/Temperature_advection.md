# Temperature Advection 폭주 문제 — 분석 및 처방 결정

> 세션 인계 문서. 다음 세션에서 이 파일만 읽으면 작업을 이어갈 수 있도록 작성.
> 마지막 갱신: 2026-05-18 (§11 구현 완료 추가)

---

## 1. 문제 현상

- **케이스**: SOPHIA_gpu, urban ABL flow, stable stratification
- **증상**: windward roof corner 부근에서 온도가 비정상적으로 cold(또는 hot)으로 폭주, 결국 시뮬레이션 발산
- **수치적 추적**: `T = ±10²¹ K` 수준까지 발산, 속도 clip(15 m/s) 도달 후 전체 cascade

## 2. 트리거 동정 (이전 세션에서 완료)

`functions_PPE.cuh`의 advection 항 bisection 결과:

```c
// functions_PPE.cuh (line 425~431, 현재 복구된 상태)
eulerx += (uxi*(uxj-uxi)*tdwx + uyi*(uxj-uxi)*tdwy + uzi*(uxj-uxi)*tdwz)*mj/rhoj*(1-P1[i].eli);
eulery += (uxi*(uyj-uyi)*tdwx + uyi*(uyj-uyi)*tdwy + uzi*(uyj-uyi)*tdwz)*mj/rhoj*(1-P1[i].eli);
eulerz += (uxi*(uzj-uzi)*tdwx + uyi*(uzj-uzi)*tdwy + uzi*(uzj-uzi)*tdwz)*mj/rhoj*(1-P1[i].eli);
if(k_con_solve && ptypej<1000 && P1[j].buffer_type!=99){
    eulert -= 2.0*mj/(rhoi+rhoj)*(tempj-tempi)*((uxi-uxj)*tdwx+(uyi-uyj)*tdwy+(uzi-uzj)*tdwz)*(1-P1[i].eli);
}
```

→ **`eulert` (central-diff SPH temperature advection)이 폭주의 직접 트리거**로 확정.
→ `eulerx/y/z` (모멘텀 advection)도 **같은 central-diff + KGC 구조** — 같은 취약성 보유.

## 3. 진단 결과 정리

### 3.1 Cold-spot 위치 (ArgMin 진단)
- **위치**: `(95.00, 0.00, 19.50)` → **windward roof corner**
- 건물: x ∈ [95.5, 104.5], roof z = 19
- `p_type = 1` (fluid), `buffer_type = 0` (정상 외부 유체)
- → **BFS 내부 마킹 결함 아님**

### 3.2 Linear-init 테스트 (옵션 A1)
- 의도: 초기 T 프로파일을 `T(0)=T_con=315`, `T(ABL_h)=T_ini=300`로 **선형 보간** (T_con 자체는 변경 X)
- 구현: `functions_initial_conditions.cuh`에 `Interpolate_T_Linear_Stable()` 추가, `functions_NNPS.cuh:40`에서 호출
- 결과:
  - step 4400 (baseline 트리거 시점)에서 tt_max 안정
  - **step 6700**: cold-spot **동일 위치 (95, 0, 19.5)에서 재발**
  - step 10000: T = -1.26e21 K로 발산
- → **초기조건은 원인이 아님. ~2300 step 지연 효과만**.

## 4. 근본 원인 분석

### 4.1 디스크리타이제이션 자체의 비단조성
중앙차분 SPH advection 형태 `(Tj - Ti)·(uij·∇W)`는:
- 부호 보존성/단조성(monotonicity) **없음**
- noise를 부호 없이 누적 → overshoot/undershoot 가능

### 4.2 코너에서의 KGC ill-conditioning
windward roof corner는:
- 위: 공기, 옆: IBM 벽, 아래: 고체 → **stencil 비대칭**
- particle disorder 큼 (정체점 근방)
- KGC 행렬이 ill-conditioned → `tdwx,tdwy,tdwz` (보정 gradient)가 noisy
- noisy gradient × central-diff (단조성 없음) = **수치 노이즈 amplification**

### 4.3 cascade 메커니즘
1. T가 한 번 overshoot/undershoot 시작
2. Boussinesq buoyancy가 비정상 T에 반응 → 비정상 수직력
3. 속도 증가 → momentum advection (eulerx/y/z, 같은 central-diff 구조)도 amplify
4. 더 큰 속도 → 더 큰 advection error → 발산

### 4.4 문헌 일관성
- **Nasar/Manchester (Eulerian ISPH)**: 정확히 같은 문제 인식, **upwind/MUSCL** 명시 도입
  - "High-order Eulerian incompressible SPH with transition to Lagrangian free-surface"
  - "High-order velocity and pressure wall BC in Eulerian incompressible SPH"
  - "Eulerian ISPH Method for Simulating Internal Flows"
- **Adami (2013)** transport-velocity, **Lind/Stansby** particle shifting, **Marrone** δ-SPH
- → 문헌 일치된 입장: **순수 central-diff SPH advection은 high-Re에서 본질적으로 불안정**, 안정화 기법 필수

## 5. 사용자 처방 기준 (이번 세션에서 결정)

1. ❌ 단순한 limiting / band-aid 싫음 — 근본 처방 원함
2. ✅ Fully Eulerian 사용 — particle shifting 같은 거 안 씀
3. ✅ 스킴 변경 가능. 단 기존 스킴 **갈아치우는 경우 함수 수정 X, 새 함수 생성** (비교용)

## 6. 결정된 처방: Riemann-SPH Advection

### 6.1 출처
- Vila 1999 ("On particle weighted methods and SPH")
- Inutsuka 2002
- Nasar et al. (Manchester Eulerian ISPH 시리즈)

### 6.2 형태 (1차 upwind, 스칼라 온도)
```
eulert -= 2·(mj/(ρi+ρj)) · T*_ij · (uij·∇W)
```
- `T*_ij`: inter-particle 방향 따라 **upwind state** 선택
- 정보 흐름 방향(`uij·∇W`의 부호)에 따라 `Ti` 또는 `Tj` 채택
- **monotone-by-construction** (1차 한정)
- limiter를 central-diff 위에 얹는 게 아니라, **이산화 자체가 Godunov framework로 교체됨** → "근본 처방" 기준 충족

### 6.3 모멘텀 확장 (나중에)
- HLL/Roe-type flux로 eulerx/y/z 동일 구조 적용 가능
- 일단 eulert만 검증 후 진행

## 7. 계산량 평가

| 항목 | 비용 |
|---|---|
| central-diff (현행) | ~10 FLOPs/pair |
| 1st-order Riemann upwind | ~11-12 FLOPs/pair (분기 1회) |
| 추가 메모리 read | 없음 (`Ti, Tj, uij, ∇W` 모두 이미 로드됨) |
| advection 항 자체 증가 | ~10% |
| **전체 step 시간 증가** | **<1%** (PPE가 지배적이라 advection 변경 영향 작음) |
| 2차 MUSCL-SPH로 확장 시 | +20 FLOPs/pair, KGC gradient 재활용, 여전히 무시 가능 |

GPU warp divergence 우려: 같은 warp 내 i들이 인접 이웃 j-loop를 돌므로 분기 패턴이 비슷해 영향 미미 예상.

## 8. 진행 계획 (다음 세션 시작점)

1. **새 함수 생성** — 기존 `eulert` 식은 그대로 둠
   - 후보 위치: `functions_PPE.cuh` 내 별도 함수, 또는 새 헤더 파일
   - 후보 이름: `KERNEL_clc_eulert_Riemann` 또는 advection 묶음
2. **1차 upwind부터 구현** — `eulert`만, eulerx/y/z는 일단 현행 central-diff 유지
3. **스위치 메커니즘** — runtime 또는 compile-time flag로 central-diff vs Riemann 비교 가능하게
4. **검증 케이스**: Transport_advection (현재 working dir의 input/Transport_advection)
   - baseline (central-diff): cold-spot at (95, 0, 19.5) 재현
   - Riemann: 코너 cold-spot 사라지는지 확인
5. **통과 시**:
   - eulerx/y/z도 같은 구조(HLL/Roe)로 확장
   - 확산 과한 경우 2차 MUSCL-SPH (Nasar 2019 스타일)로 업그레이드
6. **실패 시 대안**:
   - 2차 reconstruction을 처음부터
   - 다른 Eulerian 안정화 기법 (Riemann 외) 재검토

## 9. 현재 코드 상태 (이 세션 종료 시점)

### 수정/추가된 파일
| 파일 | 상태 |
|---|---|
| `functions_PPE.cuh` | `eulert` central-diff 복구 상태 (line 425-431) |
| `ISPH_Calc.cuh` | `cub::DeviceReduce::ArgMin` 진단 추가됨 (cold-spot idx/pos 출력) |
| `functions_initial_conditions.cuh` | `Interpolate_T_Linear_Stable()` 추가됨 |
| `functions_NNPS.cuh` (line ~40) | linear-init 진단을 위해 `Interpolate_T_Linear_Stable` 호출 중 |
| `Parameters.cuh` | `T_con = 315.0` 유지 (변경 X) |

### 다음 세션 시작 시 결정 필요
- Linear-init 진단 코드 (`Interpolate_T_Linear_Stable` 호출)를 그대로 둘지, 원래 MOST_T로 복구할지
- ArgMin 진단 출력 유지 여부 (검증 시 유용하므로 유지 권장)

## 10. 참고: 폐기된 가설 / 시도

| 가설 | 결과 |
|---|---|
| BFS 내부 마킹 결함 | ❌ cold-spot이 정상 외부 유체(buffer_type=0)에서 발생 |
| 초기 sharp T gradient가 원인 | ❌ linear-init이 폭주를 막지 못함 (지연만) |
| 단순 limiter 추가 | ❌ 사용자 기준 #1로 기각 |
| Particle shifting | ❌ 사용자 기준 #2로 기각 |
| 고차 KGC만으로 해결 | ⚠️ 근본 원인(central-diff 비단조성)은 그대로, 비추천 |

---

## 11. 구현 완료 (2026-05-18)

### 11.1 구현 전 발견: 기존 eulert 항의 구조적 결함
§4.1은 기존 `eulert`를 "central-diff SPH advection"으로 규정했으나, 구현 직전 코드 검증에서
**기존 항이 표준 advection이 아님**을 발견:

- 기존 항 `2mj/(ρi+ρj)(Tj-Ti)((ui-uj)·∇W)` 은 **상대속도** `(ui-uj)·∇W` 사용
  → **균일 유동(ui≈uj)에서 0이 됨**. 균일 유동에서도 `u·∇T≠0` 이므로 advection을 못 잡음.
- 같은 커널의 다른 advection 항은 전부 표준형 `ui·∇W` 사용:
  - `eulerx/y/z`: `Σmj/ρj(uxj-uxi)(ui·∇W)` = `ui·∇ux` ✓
  - 농도 `tmp_Rd`: `Σmj/ρj(cj-ci)(ui·∇W)` = `ui·∇c` ✓
  - 온도 `eulert`만 `(ui-uj)·∇W` — 혼자 다름
- git 첫 커밋(91d4588): 표준형 `ui·∇T`가 바로 윗줄에 **주석 처리**돼 있음.
  과거에 표준형 → 상대속도 대칭형으로 의도적으로 교체됨 (이유 불명).

### 11.2 결정: 기반 연산자 = Vila 1999 보존형
사용자 결정 — Riemann upwind를 **현행 `(ui-uj)·∇W` 위에 얹지 않고**, Vila 보존형으로 교체:
```
eulert = u·∇T_i = Σ_j  2mj/(ρi+ρj) · (T*_ij - Ti) · (u*_ij·∇W) · (1-eli)
   u*_ij = 0.5(ui+uj)               (계면 속도, Vila/Nasar 표준)
   T*_ij = (u*_ij·∇W <= 0) ? Tj : Ti  (1차 upwind / Godunov)
```
- consistency: 상수 T → T*=Ti → 0 ✓
- monotonicity: `dtemp += -eulert` 의 `(Tj-Ti)` 계수 = `-2mj/(ρi+ρj)·s ≥ 0` (s≤0일 때만 유지)
  → **monotone by construction** ✓
- central 환산: `T*=½(Ti+Tj)` 넣으면 `Σmj/(ρi+ρj)(Tj-Ti)(u*·∇W)` = 정상 `u·∇T`
- upwind 방향 검증: `s≤0` ⇔ u*가 j→i 방향 ⇔ j가 upwind ⇒ T*=Tj 채택 (정상 upwind) ✓

### 11.3 스위치 = 컴파일 타임 #define (vii 아님)
- `Parameters.cuh:104` 부근 `#define RIEMANN_EULERT` 추가
  - `0` = central-diff (현행 production, baseline 비교용) — **현재 기본값**
  - `1` = Vila 1차 upwind Riemann-SPH
- `functions_PPE.cuh` `KERNEL_advection_force3D` 내부 eulert 라인을 `#if RIEMANN_EULERT / #else / #endif`로 분기.
  → 기존 central-diff 라인은 `#else`에 **원형 그대로 보존** (비교용 기준 충족).
- 분기 오버헤드 0, 두 스킴 모두 코드에 상존.

### 11.4 컴파일 검증
- `RIEMANN_EULERT 1` 로 전체 빌드 성공 (nvcc, sm_80, exit 0).
- 기본값은 `0`으로 되돌려 둠 → 평소 `make` 시 동작 변화 없음.

### 11.5 검증 케이스 정정 (2026-05-18)
문서 §8 #4는 검증 케이스를 `Transport_advection`으로 적었으나, 코드 검증 결과 **모순 확인**:
- cold-spot `(95,0,19.5)` / 건물 `x∈[95.5,104.5]` / roof `z=19` 지오메트리는 **`CaseAS`만 만족**
  (CaseAS: x 0~99.5, z -0.5~99.5, p_type에 1000=IBM 건물 포함, DDF+Conduction+Boussinesq ON).
- 현재 `input/Transport_advection`는 9.9×0.9 작은 박스, 건물 없음(p_type {-1,1}) → 코너 폭주 재현 불가.
- 즉 폭주가 관측된 production 케이스 = **CaseAS**. Transport_advection은 스킴 sanity용 clean 박스.

→ 사용자 확정: 검증은 **CaseAS**로 진행.

### 11.6 검증 실행 (진행 중)

#### 빌드 설정 함정 (1차 런 무효 → 정정)
CaseAS 검증에는 `INPUT_CASE_NAME` 외에 **2개 컴파일 스위치를 더 맞춰야** 함:
- **`case_config.cuh:8` `CASE_ID`** — `INPUT_CASE_NAME`과 별개. CaseAS = `213`
  (`211/212/213`은 동일 case_config 블록, case명/output dir 구분용. 모두 `>210 && <220`).
  `CASE_ID 213`이어야 `functions_NNPS.cuh:37-41`의 MOST 속도/난류 프로파일 + 성층 T init이 실행됨.
  1차 런은 `CASE_ID 1500`(채널유동, RANS_LES_Validation용)인 채로 빌드 → MOST·성층 init 전부 skip,
  T가 균일 300K로 동결 → **CaseAS 검증으로 무효**.
- **`Parameters.cuh:142` `T_con`** — 지표온도 + 성층 초기 T의 ground 끝점.
  `T_con == T_ini == 300`이면 `Interpolate_T_Linear_Stable`가 z 무관 300 반환 (성층 0, MOST 열속 0).
  doc §9는 `T_con=315`라 했으나 코드는 300으로 되돌려져 있었음(문서 stale).
  CaseAS stable 검증값 = **`T_con 315.0`** (지표 315K → ABL_h=125m 상단 300K 선형 성층).

→ 정정 후 빌드(`CASE_ID 213`, `T_con 315`, `RIEMANN_EULERT 1`)에서 step 0 온도 성층 확인:
  지표 315.0K → z90~100 평균 303.6K → 선형 `315+(300-315)z/125`. 정상.

#### 실행 계획
1. **CaseAS / Riemann(`RIEMANN_EULERT 1`)** — 코너 cold-spot/발산이 사라지는지 (핵심 테스트). → §11.7 완료.
2. (선택) Transport_advection / Riemann — 단순 advection을 안 깨뜨리는지 sanity.

### 11.7 검증 결과 (2026-05-18, CaseAS / CASE_ID 213 / T_con 315 / RIEMANN_EULERT 1)

런 완주: step 10000 (t=10s), computing time 6949s. 로그 stdout 진단 + VTK temp/velocity 파싱.

**✅ 온도 폭주 = 해결.**
- `tt_max` 314.94→314.99 K, `tt_min` 302.58 불변, NaN/Inf 전 구간 0.
- windward roof corner box(x88-112, z14-26, |y|<12): T 311.99→312.04 K — 전 구간 물리적, cold-spot 미발생.
- ArgMin cold-spot 진단: 최저점이 도메인 상단 `(4,-76,103.5)`(성층 최저, 무해)로 고정 — doc `(95,0,19.5)` 코너 cold-spot 안 생김.
- baseline ±10²¹ K 발산 대비 결정적 개선. central-diff 비단조 overshoot 제거 확인. **1차 처방 성공.**

**⚠️ 단, 모멘텀 advection(eulerx/y/z, 여전히 central-diff)이 독자적으로 불안정.**
- `uz_max`: step ~5000 6.8 → 5500 12.9 → 6500 14.9 → 8500+ 15.00(clip 고정). step 9000에 |uz|>14 입자 158개.
- `ftotal_max` 25→206 증가. `fbz_max`(부력)는 0.489 안정 → 온도/부력 원인 아님.
- step 5000 최대 uz 위치 = `(94,0,19.5)` = 바로 그 windward roof corner. step 9000엔 x=400 도메인 outlet으로 확산.
- → §2·§4.3 예측대로: eulerx/y/z도 central-diff+KGC라 같은 비단조 취약성. clip이 hard 발산은 막지만 비물리적.

### 11.8 다음 세션 할 일
1. **eulerx/y/z 모멘텀 advection에 동일 Vila Riemann 구조 적용** (§8 #5).
   - 단 momentum은 스칼라 T와 달리 벡터 → upwind state를 HLL/Roe-type flux로. 1차 upwind부터.
   - 기존 central-diff는 `#if RIEMANN_EULERT`와 별도 플래그(예: `RIEMANN_EULERV`)로 보존/전환.
2. 통과 시 Transport_advection sanity, 이후 2차 MUSCL-SPH 검토.
3. ArgMin 진단·linear-init은 모멘텀까지 검증 끝나면 정리.

이후:
- 통과 시: `eulerx/y/z` 모멘텀 advection도 동일 Vila/HLL 구조로 확장 검토.
- 1차 upwind이 과확산이면 2차 MUSCL-SPH reconstruction으로 업그레이드.
- ArgMin 진단(`ISPH_Calc.cuh`), linear-init(`Interpolate_T_Linear_Stable`)은 검증 끝날 때까지 유지.

### 11.6 수정된 파일 (2026-05-18 세션)
| 파일 | 변경 |
|---|---|
| `Parameters.cuh` | `#define RIEMANN_EULERT 0` 추가 (energy 섹션 아래) |
| `functions_PPE.cuh` | `KERNEL_advection_force3D` eulert 라인 `#if/#else` 분기, Vila upwind 분기 추가 |
