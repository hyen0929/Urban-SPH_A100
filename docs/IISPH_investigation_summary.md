# Implicit PPE / 연산자 일관성 / div(u) 측정 조사 정리

SOPHIA_gpu의 PPE solver를 explicit(EISPH 1-sweep) vs implicit(N-sweep)으로 비교하고,
그 과정에서 드러난 **연산자 일관성**, **IISPH 이식**, **div(u) 측정의 함정**까지 정리.

---

## 1. 전체 흐름

1. implicit PPE 검증 파이프라인이 스캐폴딩만 돼 있고(주석 처리, helper 커널 미구현) 작동 안 함 → 완성.
2. PPE 잔차 진단(`PPE_DIAG`) 추가: `bi_max`, `resid_max`를 CSV로.
3. **explicit vs implicit(Brookshaw N-sweep) 비교 → implicit 폭주** (bi_max 400 → 11,000~22,000, 속도클램프).
4. 원인 규명: **연산자 불일치** (PPE 행렬 A ≠ 적용 압력힘의 발산).
5. **IISPH 일관 연산자**(`CONSISTENT_PPE`) 이식 → 폭주 해결.
6. early-exit, 속도/정확도 비교, div(u^{n+1}) 직접 측정 추가.
7. 측정 연산자(KGC vs plain) 문제로 "공정 비교" 논쟁.

---

## 2. 핵심 발견: implicit 폭주의 원인 = 연산자 불일치

선형계 `A·p = b`에서:
- **A** (PPE 행렬) = prep이 만드는 라플라시안.
- **적용 압력힘** = `pressureforce3D`의 gradient.

production(Brookshaw)은 **A ≠ div(적용 힘)**. 1-sweep(explicit)은 under-relaxation으로 이 불일치를
마스킹하지만, **tight하게 수렴(implicit)시키면 불일치가 드러나 발산이 피드백 폭증.**

- PCG·직접법(역행렬)도 **해결 못 함** — 전부 "어떻게 푸나(질문2)"일 뿐, "무엇을 푸나(A가 옳은가, 질문1)"는 안 고침.
  직접법은 정확해라 오히려 불일치를 더 노출.
- Brookshaw는 안정(대각우세)하지만 불일치 / 정확한 div(grad)는 일관되지만 nested-gradient 불안정 → **IISPH가 정답**.

---

## 3. IISPH 일관 연산자 이식 (`CONSISTENT_PPE`)

- A를 **적용 대칭 압력힘의 발산**으로 구성: `a_p = -Σ m_j(p_i/ρ_i² + p_j/ρ_j²)∇W`, **plain ∇W (KGC off)** → 연산자==힘==대칭.
- 커널: `KERNEL_iisph_prep` / `pass1` / `pass2` (relaxed-Jacobi, ω=0.5).
- 결과: 같은 200-sweep인데 **Brookshaw는 폭주(2.2만), IISPH는 bi_max ~330~460에서 안정.** NaN0/클램프0.
- → "implicit 폭주는 solver 속도가 아니라 연산자 일관성 문제"가 확정.

---

## 4. 진단 지표와 그 관계 (★ 핵심 개념)

```
u*        : predictor 중간속도 (압력 보정 전)
u^{n+1}   : 보정 후 최종속도

bi_max       = max|ρ·div(u*)/dt|        ← 입력 발산 (보정 "전")
div_un1      = |ρ·div(u^{n+1})/dt|       ← 출력 발산 (보정 "후") = 진짜 비압축성 위반
resid_max    = |A_brookshaw·p − b|       ← Brookshaw 잔차 (Brookshaw run만 의미)
avg_dens_err = IISPH 잔차                 ← 아래 참고
```

### avg_dens_err ≡ div(u^{n+1}) (밀도오차 = 속도발산)
IISPH 잔차를 풀면:
```
r_i = s_i − Ap_i = ρ0 − ρ^new_i = −Δt·ρ·div(u^{n+1})
→ avg_dens_err = mean|r|/Δt² = mean|ρ·div(u^{n+1})/Δt|
```
**즉 avg_dens_err는 별개 변수가 아니라, div(u^{n+1})를 IISPH 자기 연산자(plain ∇W)로 잰 것.**

### 밀도가 상수인데 왜 error가 생기나
상태의 ρ는 ρ0 고정이 맞음. 그러나 `ρ*`는 **예측 밀도** — 발산하는 u\*를 한 step 적용하면 될 뻔한 밀도:
```
ρ* = ρ0 − Δt·ρ·div(u*)   (u*가 발산하니 ρ0에서 벗어남)
s_i = ρ0 − ρ* = −Δt·ρ·div(u*)   ← "밀도오차"는 사실 속도발산의 다른 이름
```
**constant-density 비압축 솔버에서 density error ≡ velocity divergence (×Δt·ρ).** 한 물건의 두 이름.

---

## 5. 비교 결과 (★ 결론이 여러 번 뒤집힘 — 시간순)

### (a) bi_max (입력 발산)
- explicit ~74, IISPH ~202 → IISPH가 2.7× 높음.
- **단 bi_max는 입력(보정 전)이라 "implicit이 줄이는 대상"이 아님.** 비교 지표로 부적절.

### (b) div_un1 (출력 발산, KGC로 측정) — 첫 비교
- IISPH가 explicit보다 **높음** (max ~1.3×, avg ~2.9×). "implicit이 더 비압축적"은 미성립.
- 원인 의심: ① IISPH가 density-error 게이트로 trivial 수렴(iter=3) → 발산 안 줄임, ② 측정 연산자(KGC) 문제.

### (c) early-exit를 발산(div) 단위로 교정
- 기존 게이트 `avg_dens_err < 1e-3`(밀도단위)는 발산으론 ~49까지 허용 = 너무 헐거워 iter=3에 멈춤.
- **게이트를 div 단위로 바꾸니**(IISPH_ETA): 더 돌릴수록 div가 **실제로 줄어듦** (0.3 → 0.08).
  → 이전 주장 "더 돌려도 소용없다"는 **틀렸음** (잔차 r = div이라 반복이 div를 직접 줄임).
- 단 ~0.08에서 **Jacobi floor** (O(M²) 저주파 한계). 0.05 게이트면 200 cap, 0.09 게이트면 iter~13(빠름).

### (d) 속도
- IISPH(early-exit iter~13) ≈ explicit 대비 ~1.6× 느림 (blind 200-sweep이면 ~60~200×였음).

---

## 6. ★ 가장 최근 논쟁: KGC vs plain 측정 (어느 게 진짜 div냐)

같은 u^{n+1}인데 측정 연산자에 따라 div가 다름:
```
avg_dens_err (plain ∇W, IISPH 연산자)  ≈ 0.09   ← IISPH가 최소화한 값
div_un1      (KGC 보정 연산자)          ≈ 0.18~0.24
```

- **IISPH는 plain ∇W 발산만 풀었지(0.09), KGC 발산은 푼 적 없음(0.2).** 그래서 KGC로 재면 안 줄어듦.
- explicit은 KGC로 풀고 KGC로 재니 KGC값이 낮음 → "같은 발산을 다른 자로 재서 IISPH가 커 보임".

### 어느 자가 옳은가 → 해석이 뒤집힘
- **KGC는 더 정확한 발산 추정**(선형장 정확, 1차 일관). plain ∇W는 truncation 오차로 발산을 일부 놓침.
- 따라서 **IISPH의 진짜 발산 ≈ KGC값 ~0.2**, avg_dens_err ~0.09는 자기 연산자(덜 정확)로 잰 **과소평가**일 가능성.
- 즉 "KGC가 IISPH한테 불공정"이라던 초기 framing이 틀렸을 수 있음. KGC는 편파가 아니라 더 정확한 자.

### explicit 증거 (KGC≈plain)
- explicit fair run에서 explicit의 KGC-div ≈ plain-div (0.30 ≈ 0.30, ~2% 차이).
- **explicit은 두 자가 일치 → 발산이 robust.** IISPH는 두 자가 어긋남(0.09 vs 0.2) → IISPH의 낮은 plain값이 연산자 덕일 의심.

### "den/div" 질문 정리 (최근 Q&A)
- Q. div(u) 보자더니 왜 avg_dens_err가 진짜 변수냐? → **같은 것.** avg_dens_err = div(u^{n+1})를 plain ∇W로 잰 값.
- Q. 밀도 상수인데 왜 error? → ρ\*는 예측 밀도, error = 속도발산의 다른 이름.
- Q. 여러 번 반복하는데 왜 implicit이 더 큰가? → IISPH는 plain 발산을 줄였지(0.09) KGC 발산(div_un1)은 안 푼 연산자라 ~0.2로 남음. 게다가 KGC가 더 정확 → 진짜 div는 ~0.2.

### 공정 비교를 위한 조치
- `KERNEL_compute_div_un1`에 `use_kgc` 플래그 추가 → div를 KGC·plain **둘 다** 측정.
- **explicit만 재실행**(implicit은 재실행 불가, 그 plain-div는 이미 avg_dens_err에 있음) → `SOPHIA_exp_fair`, OUTPUT `CaseCS1_exp_fair`, CSV 8컬럼.
- 미결: explicit_fair가 발전유동까지 가면 KGC끼리 / plain끼리 최종 비교.

---

## 7. explicit의 시간 의존성 (정정)
초기에 "explicit div_un1_avg ~0.16 안정"이라 했으나 **틀림** — 그건 t=20~30 딥 구간만 본 것.
전체:
```
t=  0~2  : 0.510   t=  2~20 : 0.179   t= 20~61 : 0.159(딥)
t= 61~100: 0.365   t=100~150: 0.556(max 0.69, div_un1_max 스파이크 t=135에 3375)
```
explicit은 유동 발전하며 **div가 0.12 → 0.7로 증가**(악화). 그래서 implicit 우위 가능성은 **후반(t>60)** 에서
explicit이 악화될 때 — 단 implicit도 그때까지 floor를 유지하는지 확인 필요(미확인).

---

## 8. 플래그 / 파일 / 바이너리 요약

**Parameters.cuh 플래그**
| 플래그 | 의미 |
|---|---|
| `IMPLICIT_PPE` | 0=EISPH 1-sweep, 1=N-sweep Jacobi |
| `CONSISTENT_PPE` | 1=IISPH 일관 연산자 (대칭힘+plain∇W) |
| `IISPH_OMEGA` | relaxed-Jacobi 계수 (0.5) |
| `IISPH_ETA` | early-exit 게이트(**발산 단위**). 0.05=floor까지(200 iter), 0.09=빠름(iter~13) |
| `IISPH_CONV_FREQ` | step내 iter별 div 곡선 로깅 빈도. 0=off |
| `PPE_DIAG`/`PPE_DIAG_FREQ` | 진단 CSV(매 10 step) |

**진단 산출 파일**
- `<CASE>/PPE_diag_explicit.csv` : step,time,bi_max,resid_max,div_un1_max,div_un1_avg[,div_un1_plain_max,div_un1_plain_avg]
- `<CASE>/IISPH_log.txt` : step,time,iter_used,avg_dens_err(div 단위),conv
- `<CASE>/IISPH_conv.txt` : step,iter,avg_div (step내 수렴곡선, conv 켰을 때)

**구현 커널** (functions_PPE.cuh)
- `KERNEL_iisph_prep/pass1/pass2` (IISPH), `KERNEL_compute_div_un1(...,use_kgc)`, `KERNEL_fluid_mask`,
  `KERNEL_compute_ppe_residual`, `KERNEL_snapshot_pres`, `KERNEL_compute_pres_delta`.

**run 디렉토리**
- `CaseCS1_exp_div` : explicit (div_un1 KGC만, 6컬럼). t=150까지.
- `CaseCS1_iisph_div4` : IISPH 0.09 게이트, conv off. t≈32까지 확인됨.
- `CaseCS1_exp_fair` : explicit 재실행, **8컬럼(KGC+plain)** — 공정 비교용. GPU2에서 진행.

branch: `consistent-ppe`

---

## 9. 핵심 결론 (현재까지)

1. **implicit 폭주의 진짜 원인 = 연산자 불일치** (PCG/직접법으로 안 고쳐짐). IISPH 일관 연산자로 해결.
2. **반복(iteration)은 발산을 실제로 줄인다** — 단 게이트를 밀도단위가 아니라 **발산단위**로.
3. **순수 relaxed-Jacobi는 ~0.08 floor**(O(M²)). 더 낮추려면 CG/multigrid 필요.
4. **measurement 연산자가 결론을 가른다**: IISPH는 plain(자기 연산자)으론 낮지만(0.09), 더 정확한 KGC로는 ~0.2.
   → "각자 자기 연산자로" 비교는 IISPH를 봐주는 것일 수 있음. **둘 다 KGC(정확한 자)로** 보는 게 정직.
5. **미결(다음 단계)**: explicit_fair가 발전유동(t≥50, 가능하면 t>60)까지 가서 KGC끼리/plain끼리 최종 비교 →
   "implicit이 더 비압축적인가"의 답.

---

## 10. 미해결 / 의심 포인트
- implicit(iisph4) 후반(t≈32)에서 iter가 다시 200 cap·conv=0으로 갔다는 신호 → 발전유동에서 floor가 0.09 위로
  올라가 implicit이 악화될 가능성(미확인).
- IISPH의 plain self-measure가 진짜 비압축성을 과소평가하는지 — div_un1(KGC) 비교로 확정 필요.
- 속도: IISPH는 floor까지 가려면 여전히 비쌈(0.05 게이트=200 iter). 실용 구간은 0.09(iter~13, ~1.6× 느림).
