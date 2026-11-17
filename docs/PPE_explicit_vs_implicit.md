# PPE Solver 비교: Explicit vs Implicit (그리고 연산자 일치)

SOPHIA_gpu의 압력 푸아송 방정식(PPE) solver를 **explicit(EISPH 1-sweep)** 과
**implicit(N-sweep)** 으로 비교하고, 그 과정에서 드러난 **연산자 일관성(operator
consistency)** 문제와 IISPH 해법을 정리한다.

---

## 1. 무엇을 비교하려는가

| Config | 플래그 | PPE 풀이 | 비고 |
|---|---|---|---|
| **Explicit (production)** | `IMPLICIT_PPE 0` | Jacobi **1-sweep** | 현 production EISPH. 압력을 매 step 한 번만 완화 |
| **Implicit** | `IMPLICIT_PPE 1` | Jacobi **N-sweep** (max 200, early-exit) | 압력을 수렴까지 반복 |
| **IISPH (일치 연산자)** | `CONSISTENT_PPE 1` | relaxed-Jacobi N-sweep | 행렬 A를 적용 압력힘의 발산으로 구성 |

세 config 모두 같은 케이스(예: CaseCS1)·같은 입력으로 돌려서 **비압축성(incompressibility)
만족도**와 **유동 안정성**을 비교한다.

---

## 2. 왜 필요한가

production은 **EISPH 1-sweep** — PPE를 매 step 단 한 번만 완화한다. 자연스러운 가설:

> "PPE를 더 수렴시키면(implicit) 비압축성이 좋아지고 wake/코너 속도 스파이크도 줄지 않을까?"

이 가설을 **정량적으로 검증**하는 것이 1차 목적이다. 결과적으로 이 비교는
"1-sweep이 충분한가, 아니면 implicit이 필요한가", 더 나아가 "이 솔버의 PPE 이산화가
애초에 일관적인가"라는 더 근본적인 질문으로 이어진다.

---

## 3. 비교 지표 (metrics)

PPE는 선형계 `A·p = b`를 푼다. `A` = 라플라시안 연산자, `b` = 소스(= ρ·div(u*)/dt).

```
u*            : predictor 중간속도 (압력 보정 전)
u^{n+1}       : 보정 후 최종속도  (= u* + dt·a_p(p))
a_p(p)        : 압력 가속도 (pressureforce3D)
```

### 3.1 bi_max — 입력 발산 (before pressure)
```
bi_max = max | ρ·div(u*)/dt |
```
predictor가 만든 중간속도 u*의 발산. 압력으로 지워야 할 "입력 divergence". **모든
config에서 같은 Brookshaw div 연산자로 계산** → 정의가 동일.

### 3.2 resid_max — PPE 잔차 (Brookshaw)
```
resid_max = max | A·p − b |     (A = Brookshaw 라플라시안)
```
푼 압력이 PPE를 얼마나 만족하는지. **연산자가 일치할 때만** `resid_max ≈ ρ/dt·div(u^{n+1})`
의 근사가 된다. → Brookshaw를 쓰는 explicit/implicit에선 유효, IISPH에선 무의미(아래 §6).

### 3.3 div(u^{n+1}) — 출력 발산 (after pressure) [미구현]
```
div(u^{n+1}) = 보정 후 최종속도의 SPH 발산  = 진짜 비압축성 위반
```
solver 무관하게 "최종 속도장의 발산"만 직접 재는 가장 엄밀한 지표. resid_max가 하려던
일의 **정확한 버전**. 현재 PPE_DIAG엔 없음(추가하려면 양쪽 재실행 필요 → §6).

### 3.4 avg_dens_err — IISPH 밀도오차
```
avg_dens_err = mean | s_i − Ap_i | / ρ0      (s_i = ρ0 − ρ*_i)
```
IISPH 전용. 투영 후 남은 평균 밀도 위반 = IISPH의 수렴/early-exit 판정 지표.

---

## 4. 진단 인프라 (PPE_DIAG) 와 산출 파일

`PPE_DIAG 1`이면 매 `PPE_DIAG_FREQ`(=10) step마다, **post-solve 압력으로 PPE_prep을 한 번
더 돌려** bi_max·resid_max를 계산해 CSV로 남긴다. (post-solve prep 없이는 잔차가 항등적 0.)

| 파일 | 생성 조건 | 컬럼 |
|---|---|---|
| `./<CASE>/PPE_diag_explicit.csv` | `IMPLICIT_PPE 0` | step, time, bi_max, resid_max |
| `./<CASE>/PPE_diag_implicit.csv` | `IMPLICIT_PPE 1` | step, time, bi_max, resid_max |
| `./<CASE>/PPE_iter_log.txt` | implicit | step, time, iter_used, delta_p, bi_max, resid_max, conv |
| `./<CASE>/IISPH_log.txt` | `CONSISTENT_PPE 1` | step, time, iter_used, avg_dens_err, conv |

> 주의: 파일명은 `IMPLICIT_PPE`로 정해지므로 IISPH run(IMPLICIT_PPE=0)도 `PPE_diag_explicit.csv`
> 로 나온다. 비교 시 디렉토리로 구분할 것.

**단일 GPU 전용** — 진단 버퍼가 file-scope static이라 multi-GPU면 GPU-0 subdomain만 기록.

---

## 5. 실행 방법

각 config는 컴파일 타임 플래그라 **빌드를 따로** 해야 한다. 출력이 안 겹치게
`OUTPUT_CASE_NAME`을 config별로 다르게 둔다(입력은 `INPUT_CASE_NAME` 그대로).

```bash
# (1) explicit
#   Parameters.cuh: IMPLICIT_PPE 0, CONSISTENT_PPE 0, OUTPUT_CASE_NAME "CaseCS1(Explicit)"
make && CUDA_VISIBLE_DEVICES=0 ./SOPHIA_gpu 1 > CaseCS1_explicit_run.log 2>&1 &

# (2) implicit (Brookshaw N-sweep)
#   Parameters.cuh: IMPLICIT_PPE 1, CONSISTENT_PPE 0, OUTPUT_CASE_NAME "CaseCS1_implicit"
make && CUDA_VISIBLE_DEVICES=1 ./SOPHIA_gpu 1 > CaseCS1_implicit_run.log 2>&1 &

# (3) IISPH (일치 연산자)
#   Parameters.cuh: CONSISTENT_PPE 1, OUTPUT_CASE_NAME "CaseCS1_iisph"
make && CUDA_VISIBLE_DEVICES=2 ./SOPHIA_gpu 1 > CaseCS1_iisph_run.log 2>&1 &
```
- `maximum-iteration`(solv) = 200, `PPE_TOLERANCE`(implicit)/`IISPH_ETA`(IISPH) = 1e-3.
- 반드시 **단일 GPU**로.

---

## 6. 어떻게 비교하는가 (지표별 호환성)

| 지표 | explicit↔implicit | ↔ IISPH | 재실행 필요? |
|---|---|---|---|
| **bi_max** | ✅ 같은 div 정의 | ✅ 같은 div 정의 | ❌ 기존 데이터로 OK |
| **resid_max** | ✅ 둘 다 Brookshaw | ❌ IISPH는 Brookshaw 방정식을 푼 적 없음 → 무의미 | — |
| **div(u^{n+1})** | ✅ (엄밀) | ✅ (엄밀) | ⚠️ 양쪽 재실행 필요(미구현) |

### 비교 절차
1. **시간축 매칭**: 세 CSV 모두 매 10 step(dt 동일) 같은 time에 찍히므로 `time`(또는 `step`)
   컬럼으로 join.
2. **1차: bi_max 곡선 겹치기** — explicit vs implicit vs IISPH의 bi_max를 같은 time축에 overlay.
   재실행 불필요(기존 explicit run에 이미 있음). 발전된 유동 구간(예: t=10~)의 평균/추세 비교.
3. **2차: 물리장(VTK)** — wake/건물 코너의 **속도 스파이크**, 압력장 노이즈를 ParaView로 직접 비교.
4. **(선택) 엄밀**: div(u^{n+1}) 컬럼을 추가해 양쪽 재실행하면 "보정 후 발산" 직접 비교 가능.

### 왜 bi_max로 충분한가
투영이 불완전하면 보정 후 잔류 발산이 **다음 step u*로 넘어가 bi_max를 키운다**(피드백).
즉 bi_max 추세 차이가 이미 투영 품질·안정성 차이를 반영한다. 실제로 안정성 결론은
bi_max로 명확히 났다(§7).

---

## 7. 현재까지의 결과 (핵심)

CaseCS1 (4.4M 입자), 같은 200-sweep 반복 기준:

| Config | bi_max 추세 (step ~90) | 거동 |
|---|---|---|
| explicit (1-sweep) | ~75 (안정) | production, 안정 |
| **implicit (Brookshaw N-sweep)** | **400 → 11,000~22,000** | **폭주**, 속도클램프 |
| **IISPH (일치 N-sweep)** | **400 → ~460 → 감소** | 안정, 클램프 0 |

**가설 검증 결과 — 반직관적:**
- "PPE를 더 수렴시키면 좋아진다"는 **거짓**. Brookshaw 연산자를 tight하게 풀면 **오히려 폭주**.
- 원인은 solver 속도가 아니라 **연산자 불일치**: PPE 행렬 `A`(Brookshaw 라플라시안)가
  실제 적용하는 압력힘의 발산과 **다르다**. 1-sweep(explicit)은 under-relaxation으로 이를
  마스킹하지만, tight하게 수렴시키면(implicit) 불일치가 드러나 발산이 피드백 폭증.
- PCG·직접법(역행렬)도 해결 못 함 — 전부 "어떻게 푸나"일 뿐 "무엇을 푸나(A가 옳은가)"는 안 고침.
- **해법 = 연산자 일치(IISPH)**: A를 적용 압력힘의 발산으로 구성(대칭형 힘 + plain ∇W).
  같은 200-sweep인데 폭주가 사라지고 bi_max가 안정 감소. early-exit로 iter≈3에 수렴(~60× 빠름).

자세한 구현·근거: 메모 `project_iisph_consistent_operator.md`, `project_implicit_ppe_pipeline.md`.

---

## 8. 관련 파일·플래그 요약

- `Parameters.cuh`: `IMPLICIT_PPE`, `CONSISTENT_PPE`, `PPE_DIAG`/`PPE_DIAG_FREQ`,
  `PPE_TOLERANCE`, `IISPH_OMEGA`, `IISPH_ETA`, `OUTPUT_CASE_NAME`
- `functions_PPE.cuh`: `KERNEL_PPE_prep3D`/`PPE_calc3D`(Brookshaw),
  `KERNEL_compute_ppe_residual`(진단), `KERNEL_iisph_prep`/`pass1`/`pass2`(IISPH),
  `KERNEL_pressureforce3D`(CONSISTENT시 대칭형 분기)
- `ISPH_Calc.cuh`: PPE 블록 (CONSISTENT_PPE / IMPLICIT_PPE / EISPH 분기), PPE_DIAG 로깅
- solv 파일: `maximum-iteration`(=200), `minimum-iteration`
- branch: `consistent-ppe`
