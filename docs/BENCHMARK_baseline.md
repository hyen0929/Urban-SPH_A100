# SOPHIA_gpu 성능 기준선 (Baseline)

**측정일** 2026-09-02  **브랜치** `consistent-ppe`  **기준 커밋** `50f000e` + 미커밋 working tree

최적화 작업 **이전**의 성능을 기록해 둔 문서입니다. 최적화 후 §2의 프로토콜을 그대로 반복해서 이 숫자들과 비교하십시오.

---

## 1. 요약

| | CaseAS (6.9M 입자) | Kori (28.9M 입자) |
|---|---:|---:|
| 스텝당 시간 | 1.53 s (스텝 251–300) | 3.03 s (스텝 51–200, 일정) |
| GPU 시간 중 이웃탐색 커널 | 95.6 % | 90.5 % |
| **물리를 바꾸지 않고 제거 가능한 몫** | **약 31 %** | **약 34 %** |

제거 가능한 31%의 내역 (§7):

| 항목 | CaseAS | Kori | 성격 |
|---|---:|---:|---|
| `compute_vorticity3D` 를 매 스텝 실행 | 12.6 % | 11.7 % | 출력 전용 필드인데 100스텝마다만 소비됨 |
| KGC 보정행렬 커널 2개 (내용 동일 + 시간불변) | 11.9 % | 12.1 % | 한 번만 계산하면 되는 것을 매 스텝 2번 계산 |
| `PPE_DIAG` 진단 (`compute_div_un1` + 여분 `PPE_prep`) | 4.3 % | 4.4 % | `Parameters.cuh` 스위치 하나로 off |
| 출력 DtoH 중 버려지는 3/4 | 1.4 % | 4.5 % | 저장 안 하는데 매번 전량 내려받음 |

여기에 FP64 제거·점유율 개선(§7.5)이 나머지 이웃탐색 커널에 추가로 얹힙니다.

---

## 2. 재현 프로토콜

재측정 시 **반드시 동일하게** 맞춰야 하는 조건들입니다.

### 2.1 빌드

```bash
# 프로덕션 빌드 (Makefile 과 동일 옵션)
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr \
     SOPHIA_gpu.cu -o sophia_<case> -I./cub-1.8.0/ -lpthread

# CaseAS: Parameters.cuh 에서 INPUT_CASE_NAME "CaseAS", OUTPUT_CASE_NAME "BENCH_CaseAS"
#         + 컴파일 시 -DCASE_ID=213
# Kori  : Parameters.cuh 기본값 (CASE_ID 1010), OUTPUT_CASE_NAME "BENCH_Kori"
```

> 출력 디렉토리를 `BENCH_*` 로 분리한 이유: 기존 `Kori_vy1/` 등 실제 연구 결과를 덮어쓰지 않기 위함.

### 2.2 실행

```bash
# 반드시 전용 GPU 에서. 다른 사람 작업과 SM 을 공유하면 수치가 무의미해집니다.
nvidia-smi --query-gpu=index,utilization.gpu --format=csv,noheader   # 대상 GPU 가 0% 인지 확인
CUDA_VISIBLE_DEVICES=4 ./sophia_<case> 1 > run.log 2>&1
```

### 2.3 측정 창 — 중요

**스텝당 비용이 정상상태가 아닙니다.** CaseAS 에서 733 ms/step (51–100 스텝) → 1527 ms/step (251–300 스텝)으로 **2.1배** 증가합니다. 유동이 발달하면서 `KERNEL_advection_force3D` 내부의 조건부 항(인공점성 `uij_xij<0`, 인공열확산 `uij_xij_t<0 || flt_s<0.7`)이 점점 더 많은 이웃 쌍에서 활성화되기 때문으로 보입니다.

따라서 **단일 숫자로 비교하지 말고, 반드시 같은 스텝 구간을 비교**하십시오. 본 문서의 커널 표는 모두 **CaseAS 스텝 약 250–310 구간**입니다.

### 2.4 프로파일

```bash
# 커널별 GPU 시간 (경합 없는 전용 GPU 전제)
nsys profile --trace=cuda --delay=150 --duration=90 -o prof ./sophia_CaseAS 1
nsys stats --report cuda_gpu_kern_sum   --format csv -o - prof.nsys-rep
nsys stats --report cuda_gpu_mem_time_sum --format csv -o - prof.nsys-rep
nsys stats --report cuda_gpu_mem_size_sum --format csv -o - prof.nsys-rep
```

---

## 3. 측정 환경

| 항목 | 값 |
|---|---|
| GPU | NVIDIA A100 80GB PCIe (device 4, 전용 — 측정 중 타 프로세스 0건) |
| CUDA | 12.1 (V12.1.66), nsys 2023.1.2 |
| 컴파일 | `-w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr` |
| 호스트 | 96 core, 1007 GB RAM |
| 실행 | 단일 GPU (`./SOPHIA_gpu 1`) |

> **`ncu` (Nsight Compute) 는 사용하지 못했습니다** — `ERR_NVGPUCTRPERM`, 일반 사용자에게 GPU 성능 카운터 접근 권한이 없습니다. 점유율·달성 대역폭·런타임 FP64 명령수 같은 하드웨어 카운터가 필요하면 관리자에게 `NVreg_RestrictProfilingToAdminUsers=0` 설정을 요청해야 합니다. 그 전까지는 §7.5의 **정적 SASS 분석**이 대체 근거입니다.

---

## 4. 벤치마크 케이스

| | CaseAS | Kori |
|---|---:|---:|
| CASE_ID | 213 | 1010 |
| 입자 수 | 6,908,875 | 28,867,981 |
| LDM 입자 수 | 250,000 | 400,000 |
| dt | 1e-3 s | 1e-2 s |
| freq_LDM | 1 (매 스텝) | 20 |
| 난류 모델 | DDF | DDF |
| IBM | `DEACTIVATE_IBM=1` (건물 제거) | 활성 |
| `part1` 배열 | 1,370 MB | 5,726 MB |
| GPU 메모리 (P1+SP1+P2+SP2+P3) | 4.4 GB | 18.3 GB |

---

## 5. CaseAS 결과 (6.9M 입자)

### 5.1 Wall-clock

`Computing time` 출력(50스텝 주기, `clock()` 기준) 차분:

| 스텝 구간 | 소요 | ms/step |
|---|---:|---:|
| 1–50 (초기화 포함) | 38.97 s | — |
| 51–100 | 36.66 s | 733.2 |
| 101–150 | 36.92 s | 738.4 |
| 151–200 | 48.11 s | 962.2 |
| 201–250 | 74.27 s | 1485.4 |
| 251–300 | 76.33 s | 1526.6 |

입력 파일(874 MB) 로딩에 별도로 약 100 s.

### 5.2 커널별 GPU 시간 — 스텝 250–310 구간 (nsys, 90 s / 58 스텝)

| 항목 | 호출/스텝 | 1회 평균 | ms/step | 비중 |
|---|---:|---:|---:|---:|
| `KERNEL_advection_force3D` | 1.00 | 377.5 ms | 377.5 | 25.5% |
| `KERNEL_PPE_prep3D` | 1.10 | 179.1 ms | 197.7 | 13.3% |
| `KERNEL_compute_vorticity3D` | 1.00 | 186.8 ms | 186.8 | 12.6% |
| `KERNEL_pressureforce3D` | 0.98 | 184.9 ms | 181.7 | 12.3% |
| `KERNEL_clc_prep3D` | 1.00 | 148.5 ms | 148.5 | 10.0% |
| `KERNEL_clc_correction_KGC_3D` | 1.00 | 98.9 ms | 98.9 | 6.7% |
| `KERNEL_clc_gradient_correction_3D` | 1.00 | 76.7 ms | 76.7 | 5.2% |
| `IBM_spreading_interpolation3D` | 1.00 | 51.5 ms | 51.5 | 3.5% |
| `KERNEL_time_update_projection` | 1.00 | 49.1 ms | 49.1 | 3.3% |
| `KERNEL_compute_div_un1` | 0.21 | 208.5 ms | 43.1 | 2.9% |
| `[memcpy DtoH]` | 1.00 | 28.9 ms | 28.9 | 1.9% |
| `KERNEL_Neumann_boundary3D` | 2.00 | 10.3 ms | 20.7 | 1.4% |
| `KERNEL_clc_projection` | 1.00 | 5.2 ms | 5.2 | 0.4% |
| `[memset]` | 8.41 | 0.4 ms | 3.0 | 0.2% |
| `[memcpy DtoD]` | 3.00 | 0.7 ms | 2.1 | 0.1% |
| `KERNEL_periodic_boundary3D` | 1.00 | 1.9 ms | 1.9 | 0.1% |
| **합계 (GPU busy)** | | | **1483** | **100%** |

- GPU busy 합계 1,483 ms/step vs wall-clock 1,527 ms/step → **GPU 가동률 97 %**. CPU 측 오버헤드나 커널 런치 갭은 3 % 수준이며, 병목은 전적으로 커널 자체.
- 이웃탐색 루프를 도는 커널들의 합 = **95.6 %**.

### 5.3 메모리 오퍼레이션 (같은 58 스텝 창)

| 오퍼레이션 | 총 시간 | 총 전송량 | 횟수 | ms/step | 실효 대역폭 |
|---|---:|---:|---:|---:|---:|
| `memcpy DtoH` (출력) | 1.674 s | 6,097 MB | 58 | 28.9 | **3.6 GB/s** |
| `memset` (P3 초기화 등) | 0.175 s | 77,964 MB | 488 | 3.0 | 446 GB/s |
| `memcpy DtoD` (SP1←P1, SP2←P2) | 0.119 s | 99,862 MB | 174 | 2.1 | 839 GB/s |

- DtoH 실효 대역폭이 **3.6 GB/s** 로 매우 낮습니다. `file_P1/P2/P3` 가 `malloc`(pageable) 이라 스테이징 버퍼를 거치기 때문입니다. `cudaHostAlloc` (pinned) 으로 바꾸면 이 구간만 3–5배 빨라집니다.
- 단일 최대 전송 395 ms (1,437 MB = `part1` 전체).

---
## 6. Kori 결과 (28.9M 입자)

### 6.1 Wall-clock

`time_end=2.0 s`, `dt=0.01` → 200 스텝 완주. CaseAS 와 달리 **스텝당 비용이 거의 완전히 일정**합니다 (표준편차 0.1 % 미만). dt 가 10배 크고 총 스텝이 200 뿐이라 유동이 아직 §2.3 의 램프 구간에 못 들어간 것으로 보입니다.

| 스텝 구간 | 소요 | ms/step |
|---|---:|---:|
| 1–50 (초기화 포함) | 164.05 s | — |
| 51–100 | 151.72 s | 3034 |
| 101–150 | 151.29 s | 3026 |
| 151–200 | 151.40 s | 3028 |

입력 파일(3.6 GB) 로딩 + 초기화에 별도로 약 70 s.

### 6.2 커널별 GPU 시간 (nsys, 전체 201 스텝)

| 항목 | 호출/스텝 | 1회 평균 | ms/step | 비중 |
|---|---:|---:|---:|---:|
| `KERNEL_advection_force3D` | 1.00 | 651.7 ms | 652 | 21.2% |
| `KERNEL_PPE_prep3D` | 1.10 | 352.5 ms | 389 | 12.6% |
| `KERNEL_pressureforce3D` | 1.00 | 366.7 ms | 367 | 11.9% |
| `KERNEL_compute_vorticity3D` | 1.00 | 361.4 ms | 361 | 11.7% |
| `KERNEL_clc_prep3D` | 1.00 | 294.9 ms | 295 | 9.6% |
| `KERNEL_clc_correction_KGC_3D` | 1.00 | 205.5 ms | 206 | 6.7% |
| `[CUDA memcpy DtoH]` | 1.08 | 169.8 ms | 184 | 6.0% |
| `KERNEL_clc_gradient_correction_3D` | 1.00 | 165.0 ms | 165 | 5.4% |
| `IBM_spreading_interpolation3D` | 1.00 | 115.4 ms | 115 | 3.7% |
| `KERNEL_time_update_projection` | 1.00 | 93.6 ms | 94 | 3.0% |
| `KERNEL_Neumann_boundary3D` | 2.00 | 46.5 ms | 93 | 3.0% |
| `KERNEL_compute_div_un1` | 0.21 | 410.1 ms | 86 | 2.8% |
| `IBM_force_interpolation3D` | 1.00 | 11.7 ms | 12 | 0.4% |
| `KERNEL_clc_projection` | 1.00 | 9.4 ms | 9 | 0.3% |
| `[CUDA memcpy DtoD]` | 3.94 | 2.2 ms | 8 | 0.3% |
| `KERNEL_boundary3D` | 1.00 | 5.6 ms | 6 | 0.2% |
| `IBM_corrector_fluid` | 1.00 | 5.3 ms | 5 | 0.2% |
| `KERNEL_set_ncell` | 1.00 | 5.2 ms | 5 | 0.2% |
| `KERNEL_convective_outflow3D` | 1.00 | 4.9 ms | 5 | 0.2% |
| `KERNEL_PPE_calc3D` | 1.00 | 4.6 ms | 5 | 0.1% |
| `KERNEL_sponge_ceiling` | 1.00 | 3.8 ms | 4 | 0.1% |
| `KERNEL_freeslip_boundary3D` | 1.00 | 3.7 ms | 4 | 0.1% |
| `KERNEL_periodic_boundary3D` | 1.00 | 3.6 ms | 4 | 0.1% |
| `[CUDA memset]` | 2.04 | 1.6 ms | 3 | 0.1% |
| **합계 (GPU busy)** | | | **3079** | **100%** |


- GPU busy 합계 3,079 ms/step vs wall-clock 3,028 ms/step → **가동률 100 % 수준**. (합계가 wall 보다 약간 큰 것은 스텝 0 의 NNPS 정렬·초기화가 201 개 인스턴스에 섞여 들어간 탓.)
- CaseAS 대비 입자 수 4.18배 → 스텝 시간 1.98배. 커널이 대역폭 바운드에 가까워 선형보다 잘 확장됩니다.

### 6.3 메모리 오퍼레이션

| 오퍼레이션 | 총 시간 | 총 전송량 | ms/step | 실효 대역폭 |
|---|---:|---:|---:|---:|
| `memcpy DtoH` (출력) | 37.02 s | 113,739 MB | **184** | **3.1 GB/s** |
| `memcpy DtoD` (SP1←P1, SP2←P2) | 1.71 s | 1,428,772 MB | 8.5 | 837 GB/s |
| `memset` (P3 초기화 등) | 0.65 s | 1,145,609 MB | 3.2 | 1766 GB/s |

**단일 DtoH 전송 최대 4.37 초** (6,005 MB = `part1` 전체를 pageable 버퍼로). CaseAS(1.9 %)와 달리 Kori 에서는 출력 DtoH 가 **6.0 %** 로 올라옵니다 — 입자 수에 비례해 커지므로 대용량일수록 더 아픕니다.

---

## 7. 병목 해석

### 7.1 `compute_vorticity3D` — CaseAS 12.6 % / Kori 11.7 %, 순수 낭비

`P3.omx/omy/omz` 는 **VTK 출력 함수에서만 읽힙니다** (`functions_OUTPUT.cuh` 외 소비처 없음 — 확인함). 그런데 커널은 매 스텝 실행되고, 실제 소비는 `save_vtk_bin_single_fluid_condition` (100스텝) / `save_vtk_bin_single_vorticity` (1000스텝) 뿐입니다.

→ `ISPH_Calc.cuh:619` 의 launch 를 출력 memcpy 와 같은 주기(`count%25==0`)로 묶으면 비용이 1/25 로 줍니다 (**약 12 % → 0.5 %**).

### 7.2 KGC 보정행렬 — CaseAS 11.9 % / Kori 12.1 %, 중복 + 시간불변

두 커널이 **완전히 동일한 대칭행렬**을 계산해 표현만 다르게 저장합니다:

| 커널 | 저장 위치 | CaseAS | Kori |
|---|---|---:|---:|
| `KERNEL_clc_correction_KGC_3D` (`functions_PREP.cuh:213`) | `P3.Cm[3][3]` | 98.9 ms | 205.5 ms |
| `KERNEL_clc_gradient_correction_3D` (`functions_KNL.cuh:118`) | `P3.inv_cm_*` | 76.7 ms | 165.0 ms |

대응 관계: `Cm[0][0]==inv_cm_xx`, `Cm[0][1]==Cm[1][0]==inv_cm_xy`, `Cm[1][2]==Cm[2][1]==inv_cm_yz`, `Cm[0][2]==Cm[2][0]==inv_cm_zx`, `Cm[1][1]==inv_cm_yy`, `Cm[2][2]==inv_cm_zz`.

게다가 이 행렬은 **위치·`h`·`m`·`rho` 에만 의존**하고, 이 네 가지 모두 런 내내 불변입니다:
- 위치: Eulerian (`eli` 가 입력에 없어 0 → `KERNEL_time_update_projection` 의 위치 갱신항이 `*(P1[i].eli)` 로 소거)
- `rho`: 입력 파싱과 초기화 이후 **어디서도 갱신되지 않음** (`grep` 확인)

→ 커널 하나로 통합 + `count==0` 에서만 실행하면 **약 12 % → 0 %**.

### 7.3 `PPE_DIAG` — 약 4.4 %, 진단 스위치

`Parameters.cuh:142` 의 `#define PPE_DIAG 1` 때문에 10스텝마다 `compute_div_un1` 2회(2.9 %) + 여분 `PPE_prep3D` (1.3 %) + `compute_ppe_residual` 이 돕니다. explicit/implicit 비교가 끝났다면 `0` 으로 두십시오.

### 7.4 출력 DtoH — CaseAS 1.9 % / **Kori 6.0 %**, 그중 3/4 는 버려짐

`ISPH_Calc.cuh:721` 이 `count%25==0` 마다 P1+P2+P3 를 전부 내려받지만, 실제 저장은 `%100`(fluid_condition) / `%1000`(vorticity) / `%10000`(fluid) 에서만 일어납니다. **4회 중 3회는 전송 후 그냥 버립니다.**

→ memcpy 를 각 저장 분기 안으로 이동 + pinned 버퍼. `file_P3` 는 vorticity 저장에만, `file_P2` 는 사실상 불필요.

**입자 수에 비례해 비중이 커집니다** (6.9M → 1.9 %, 28.9M → 6.0 %). 실효 대역폭이 3.1 GB/s 밖에 안 나오는 것은 호스트 버퍼가 `malloc`(pageable) 이기 때문이며, `cudaHostAlloc` 으로 바꾸면 이 부분만 3–5 배 개선됩니다.

### 7.5 FP64 오염 — 남은 이웃탐색 커널 전반

`Real=float` 인데 이웃 루프 안의 `2.0` / `0.5` / `0.01` / `1e-20` 같은 **접미사 없는 부동소수 리터럴이 double** 이라, `2.0*rhoi*rhoj/(rhoi+rhoj)` 같은 식이 통째로 FP64 로 컴파일됩니다. A100 은 FP64 FMA 가 FP32 의 1/2, **FP64 나눗셈은 소프트웨어 합성**이라 훨씬 비쌉니다 (`-use_fast_math` 는 float 에만 적용).

현재 바이너리의 정적 SASS 계수 (ncu 불가로 정적 분석):

| 커널 | 총 명령 | FP64 | 배정밀도 나눗셈 | REG |
|---|---:|---:|---:|---:|
| `advection_force3D` | 6,272 | 735 | 36 | 112 |
| `clc_prep3D` | 5,488 | 657 | 30 | 95 (+STACK 40, 스필) |
| `PPE_prep3D` | 3,936 | 542 | 31 | 90 |
| `pressureforce3D` | 3,968 | 543 | 33 | 64 |
| `compute_vorticity3D` | 3,712 | 506 | 28 | 78 |

**검증 실험** — 리터럴을 `Real` 에 묶고(`(Real)2.0` 또는 UDL `2.0_r`), `calc_kernel_wij/dwij/tmpA` 를 실사용 커널(Wendland6) 전용으로 특수화한 뒤 재컴파일:

| 커널 | 총 명령 | FP64 | dbl-div | REG |
|---|---:|---:|---:|---:|
| `advection_force3D` | 6,272 → **1,520** | 735 → 65 | 36 → 4 | 112 → **91** |
| `clc_prep3D` | 5,488 → **1,216** | 657 → 65 | 30 → 4 | 95 → **68** (스필 해소) |
| `PPE_prep3D` | 3,936 → **576** | 542 → 0 | 31 → 0 | 90 → **56** |
| `pressureforce3D` | 3,968 → **688** | 543 → 0 | 33 → 0 | 64 → **40** |

명령 수 감소분의 상당 부분은 실행되지 않는 죽은 분기(Gaussian/Quintic/Quartic 커널의 `pow(double,double)`)라 그대로 속도가 되지는 않습니다. 다만 **레지스터 90→56 은 점유율 31%→50% 로 직결**되고, 실제로 실행되던 FP64 나눗셈 제거분이 그 위에 얹힙니다.

> 주의: `2.0f` 로 하드코딩하지 마십시오. `Real` 을 `double` 로 바꿔 쓸 여지가 없어집니다. `(Real)2.0` 캐스트나 `constexpr Real operator"" _r(long double)` 사용자 정의 리터럴을 쓰면 두 빌드 모두에서 정확히 원하는 코드가 나옵니다 (SASS 로 확인함).
>
> 참고: 현재 코드에는 이미 `f` 접미사 리터럴 257개와 `fabsf`/`fmaxf`/`sqrtf` 계열 64개가 박혀 있어 **`Real=double` 스위치가 사실상 깨져 있습니다** (컴파일은 되지만 그 지점들에서 조용히 float 로 잘림).

---

## 8. 이전 코드리뷰 대비 정정

측정 결과 사전 코드리뷰에서 제시한 우선순위 일부를 바로잡아야 합니다.

| 항목 | 리뷰 당시 추정 (Kori 기준) | 실측 (Kori) | 판정 |
|---|---|---|---|
| 스텝당 구조체 복사 (`SP1←P1`, `SP2←P2`) + `P3` memset | "스텝당 ~6 ms, 유의미" | 11.7 ms/step, **전체의 0.4 %** | 절대값은 비슷하나 커널이 워낙 무거워 **비중은 무시할 수준**. 우선순위 하향. |
| 출력 DtoH | "스텝당 ~100 ms" | 184 ms/step, **6.0 %** (그중 3/4 낭비) | **오히려 과소평가**. Kori 규모에서는 실질 병목. 우선순위 유지. |
| `cudaDeviceSynchronize` 81 개 | "스텝당 ~0.3 ms" | GPU 가동률 ~100 % → CPU 갭 거의 없음 | 효과 미미. 우선순위 하향. |
| 이웃탐색 커널이 지배적 | "이웃 루프 12개" | 커널이 GPU 시간의 94 %, 그중 이웃 루프 90.5 % | 확인됨. |

**측정 전에는 크기를 알 수 없었던 항목**이 훨씬 큽니다. `compute_vorticity3D` (11.7 %) 와 KGC 중복 (12.1 %) 은 정적 코드 분석만으로는 비중을 추정할 수 없었고, 둘 다 물리를 전혀 건드리지 않고 제거됩니다.

**스케일 의존성**: 출력 DtoH 비중은 입자 수에 비례해 커집니다 (6.9M → 1.9 %, 28.9M → 6.0 %). 반대로 커널 시간은 입자 수 4.18배에 1.98배만 늘어 대역폭 바운드 특성을 보입니다. 따라서 **CaseAS 결과만으로 Kori 규모의 병목을 판단하면 안 됩니다.**

---

## 9. 최적화 후 재측정 체크리스트

1. §2.1 과 동일한 옵션으로 빌드했는가
2. 전용 GPU (`utilization.gpu == 0%`) 에서 돌렸는가
3. **동일 스텝 구간**(CaseAS 250–310)에서 nsys 를 떴는가
4. 물리 결과가 보존되는가 — `tu_max`, `tt_max/min`, `ftotal_max` 를 기존 로그와 대조. 특히 FP64→FP32 변경 후에는 필수
5. `BENCH_*` 출력 디렉토리를 썼는가 (실제 연구 결과 보호)


---

## 부록 A. 권장 작업 순서 (실측 근거 기준)

| # | 작업 | 예상 이득 (Kori) | 물리 영향 | 난이도 |
|---|---|---:|---|---|
| 1 | `compute_vorticity3D` 를 출력 주기로 gating | −11.2 % | 없음 (출력 전용 필드) | 5분 |
| 2 | KGC 커널 통합 + `count==0` 에서만 실행 | −12.1 % | 없음 (시간불변 값) | 1시간 |
| 3 | `PPE_DIAG 0` (비교 실험 종료 시) | −4.4 % | 없음 (진단) | 1분 |
| 4 | 출력 memcpy 를 저장 분기 안으로 + pinned 버퍼 | −4.5 % | 없음 | 30분 |
| 5 | FP64 제거 (`(Real)2.0` / `_r`) + 커널함수 특수화 | 미측정 (§7.5) | 미세 (float 중간연산) | 1–2시간 |
| 6 | `float4 pos_h[]` 분리 → 이웃 리스트 캐싱 | 미측정 | 없음 | 큼 |

1–4 만으로 **약 32 %** 단축이 예상되며 전부 물리 결과를 바꾸지 않습니다. 5번은 §2.3 프로토콜로 반드시 회귀 확인이 필요합니다.

## 부록 B. 산출물

| 파일 | 내용 |
|---|---|
| `BENCH_CaseAS/`, `BENCH_Kori/` | 벤치마크 실행이 남긴 VTK (기존 연구 결과와 분리) |
| nsys 리포트 | 재현하려면 §2.4 명령을 다시 실행 |

원시 로그와 nsys 리포트는 세션 스크래치패드에 있으며 영구 보관되지 않습니다. 장기 보관이 필요하면 §2.4 를 재실행해 저장하십시오.
