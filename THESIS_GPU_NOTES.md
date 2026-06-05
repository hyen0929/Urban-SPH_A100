# SOPHIA_gpu — GPU 가속 구현 기술 노트 (학위논문 2.3장 작성용)

> 본 문서는 학위논문 **2.3 GPU-Accelerated Implementation** 절(2.3.1 Parallelization strategy and data layout / 2.3.2 Neighbor-search algorithm and memory management)을 작성하기 위해 SOPHIA_gpu의 구현 세부를 정리한 자료이다. 두 절에 직접 사용할 수 있는 내용 외에, 절 구성 변경 시 활용할 수 있도록 IBM, PPE, 시간 적분, LDM 연성 등 GPU 관점에서 의미 있는 항목까지 포함한다.

---

## 0. 한눈에 보는 설계 요약

| 항목 | 채택안 | 비고 |
|---|---|---|
| 기본 정밀도 | `float` (`typedef float Real`) | A100/RTX30의 단정도 처리량 활용 |
| 컴파일 모델 | 단일 변환단위(single TU) | 모든 `.cuh`를 `SOPHIA_gpu.cu`에 직접 `#include` |
| 컴파일러/플래그 | `nvcc -arch=sm_80 -O3 -use_fast_math -expt-relaxed-constexpr` | A100/Ampere 기준; CUB 1.8.0 헤더 동봉 |
| 멀티 GPU | x축 도메인 분할 + pthread 1-스레드/GPU + P2P | `pthread_barrier_t`로 스텝 동기 |
| 입자 데이터 구조 | 3개 구조체로 분리(AoS) — `part1/part2/part3` | hot/cold/scratch 분리 |
| 이웃 탐색 | 셀 해시 + CUB Radix Sort + 셀 시작/끝(g_str/g_end) | 스텝마다 재정렬 |
| 정수/실수 솔버 파라미터 | 64/32 슬롯 배열 → `__constant__` 복사 | `vii[]`, `vif[]` |
| 커널 launch 표준 형 | `threads.x = 128`, `blocks.x = (N-1)/128 + 1` | warp 4개/블록, 점유율 우선 |
| 압력 해법 | Projection-based ISPH, 반복 가능한 행렬-자유(matrix-free) Jacobi 형태 | 한 스텝당 1 반복 기본 |
| 시간 적분 | Predictor–Corrector (projection) | Pre_Cor / Euler 선택 가능 |
| 출력 | VTK Binary POLYDATA, `freq_output` 스텝마다 D2H copy | ParaView 가시화 |

---

## 1. 코드 구조와 컴파일 모델

### 1.1 단일 변환단위 컴파일

[Makefile](Makefile)는 단 한 줄이다.

```bash
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr \
     SOPHIA_gpu.cu -o SOPHIA_gpu -I./cub-1.8.0/ -lpthread
```

[SOPHIA_gpu.cu](SOPHIA_gpu.cu#L110-L131)의 마지막 부분에서 모든 헤더가 한꺼번에 `#include`되어, **별도의 분리 컴파일 없이 nvcc 한 번으로 전체 솔버가 단일 객체로 빌드**된다. 이 방식의 효과:

- 인라인/링크 시간 최적화(LTO 없이도) — 핵심 inline `__device__` 함수(커널 함수, 거리 계산, gradient correction 등)가 호출처에서 자연스럽게 인라이닝
- 빌드 의존성 그래프가 단순 → CI/스크립트화가 용이
- 단점: 증분 빌드가 사실상 불가능, 헤더 한 줄 변경 시 전체 재컴파일

### 1.2 디렉터리/파일별 역할

[SOPHIA_gpu.cu](SOPHIA_gpu.cu) — `main()`, 전역 변수 선언, 헤더 통합 인클루드.
주요 헤더:

| 파일 | 라인수 | 역할 |
|---|---:|---|
| [Variable_Type.cuh](Variable_Type.cuh) | 5 | `typedef float Real; typedef int int_t;` — 정밀도 정책 |
| [Parameters.cuh](Parameters.cuh) | 163 | 컴파일 타임 상수, CASE_NAME, 입출력 경로, MOST/LDM 상수 |
| [class_Cuda_Particle_Array.cuh](class_Cuda_Particle_Array.cuh) | 169 | `part1/part2/part3`, `L_part1/L_part2`, `p2p_part3` 정의 |
| [function_init.cuh](function_init.cuh) | 1113 | `vii[]/vif[]` 매크로 alias, solv 파일 파서, 입자 입력 파싱 |
| [case_config.cuh](case_config.cuh) | 318 | CASE_ID별 도메인/지형 설정 |
| [functions_NNPS.cuh](functions_NNPS.cuh) | 777 | 셀 인덱싱·정렬·재배치, 멀티 GPU exchange 커널 |
| [functions_KNL.cuh](functions_KNL.cuh) | 290 | SPH 커널(W, ∇W) — Wendland2/4/6, Quintic, Quartic, Gaussian |
| [functions_PREP.cuh](functions_PREP.cuh) | 629 | 변형률, 난류 모델, KGC, 색깔함수 등 SPH 사전계산 |
| [functions_IBM.cuh](functions_IBM.cuh) | 549 | 침지경계법(IBM) 3단계: interpolation → spreading → corrector |
| [functions_PPE.cuh](functions_PPE.cuh) | 885 | Pressure Poisson Equation prep/계산/Neumann 경계 |
| [functions_TIME_ISPH.cuh](functions_TIME_ISPH.cuh) | 431 | Projection 기반 predictor/corrector, 위치/속도 갱신 |
| [functions_BC.cuh](functions_BC.cuh) | 1621 | no-slip, periodic, free-slip, open-boundary |
| [functions_MOST.cuh](functions_MOST.cuh) | 310 | Monin-Obukhov 표면층 경계조건 |
| [functions_LDM.cuh](functions_LDM.cuh) | 1149 | Lagrangian dispersion 입자 전송 |
| [functions_DD.cuh](functions_DD.cuh) | 191 | 건식 침적 / 재부유 |
| [functions_Diffusion.cuh](functions_Diffusion.cuh) | 1451 | Eulerian 스칼라 수송 (BiMacCormack 포함) |
| [functions_OUTPUT.cuh](functions_OUTPUT.cuh) | 3427 | VTK 바이너리 출력, 변수 선택 로직 |
| [ISPH_Calc.cuh](ISPH_Calc.cuh) | 836 | **시간 루프 본체** + 메모리 할당, GPU 셋업 |
| [ISPH.cuh](ISPH.cuh) | 187 | GPU device 조회, 입자 read, NI/NJ/NK 계산, pthread 분기 |

전체 약 **16.6k 라인** 규모. 약 8할이 CUDA `__global__`/`__device__` 코드.

---

## 2. 병렬화 전략 (논문 2.3.1용)

SOPHIA_gpu는 **2단계의 병렬화 계층**을 갖는다.

### 2.1 1차 병렬화 — 멀티 GPU 도메인 분할

[ISPH.cuh:115](ISPH.cuh#L115) — `calc_area = ceil(NI/ngpu);`

x축 방향으로만 단순 분할(slab decomposition):

```
GPU 0     GPU 1     GPU 2     ...     GPU ngpu-1
[--------][--------][--------]        [--------]
 calc_area cells     calc_area
```

- y, z 축은 분할하지 않는다. 도시 풍환경 시뮬레이션의 경우 **풍방향(x축)이 가장 길고, ABL 정상화에 따라 하류로 갈수록 통계량 수렴**하기 때문에 x축 슬랩이 자연스럽게 부하균형에 가깝다.
- 각 GPU는 호스트 측 [ISPH.cuh:166-169](ISPH.cuh#L166-L169)에서 **pthread 한 개**로 띄워진다.
  ```c
  for(i=0;i<ngpu;i++){
      tid[i]=i;
      pthread_create(&solve_thread[i], NULL, ISPH_Calc, (void*)&tid[i]);
  }
  ```
- 동기화는 단일 글로벌 `pthread_barrier_t barrier`로 처리. 한 스텝 내 데이터 교환 직전·직후에 `pthread_barrier_wait(&barrier)` 호출로 모든 GPU의 같은 페이즈 도달을 보장.

#### 2.1.1 입자 i_type 표지

[functions_NNPS.cuh:443-571](functions_NNPS.cuh#L443-L571) `initial_particle` / `inner_outer_particle` 커널은 각 입자의 위치를 보고 다음의 4 종류 라벨을 매긴다.

| `i_type` | 의미 |
|---|---|
| 0 | inner: 다른 GPU와 무관하게 자체 계산 가능 |
| 1 | inner-buffer: 옆 GPU에 고스트로 보낼 후보(자기 영역 안쪽 2 셀) |
| 2 | outer-buffer: 옆 GPU에서 받은 고스트(영역 밖, 2 셀 두께) |
| 3 | dummy: 메모리 슬롯 차지만 함 |

이 라벨은 모든 SPH 커널의 **early return** 분기로 사용된다 — 예: PPE prep의 [functions_PPE.cuh:5](functions_PPE.cuh#L5) `if(P1[i].i_type>i_type_crt) return;`. 즉 dummy/outer-buffer는 자기 자신의 RHS를 갱신하지 않고, 다른 입자의 stencil에만 기여한다.

#### 2.1.2 P2P 교환

[ISPH.cuh:144-150](ISPH.cuh#L144-L150)에서 모든 GPU 쌍에 `cudaDeviceEnablePeerAccess`를 설정한다. 이후 인접 GPU 간 입자 전송은 **호스트 메모리를 거치지 않고 GPU↔GPU 직통**(NVLink/PCIe P2P).

[functions_NNPS.cuh:326-410](functions_NNPS.cuh#L326-L410)의 `right_send_particle`, `left_send_particle`, `reorder_data_p2p`, `p2p_copyData` 커널이 다음 한 사이클을 구성한다.

1. `left/right_send_particle`: 보낼 입자(`xi >= maxb − 2·dcell` 또는 `xi < minb + 2·dcell`)를 골라 송신용 SoA(`send_P1`/`send_P3`)에 깃발(`p2p_af=1`) + 인덱스(`p2p_idx`) 기록.
2. CUB radix sort로 `p2p_af` 키 기준 정렬 → 보낼 입자가 앞으로 압축.
3. `reorder_data_p2p`: 압축된 슬롯에 입자 데이터 모음.
4. `p2p_copyData`: 수신 GPU의 `dev_SP1` 끝부분(slot `num_part2-i-1`)에 직접 기록.

전송 통신량 추정: 각 GPU 2·dcell 두께 슬랩 × YZ 단면적의 입자 수 × `sizeof(part1)+sizeof(p2p_part3)`. 이 때문에 [ISPH.cuh:119](ISPH.cuh#L119)에서 `num_p2p = num_part*4/NI*C_p2p` (`C_p2p=1.5`)로 송수신 버퍼 크기를 보수적으로 잡아둔다.

#### 2.1.3 스트림과 단계 중첩

[ISPH_Calc.cuh:416-417](ISPH_Calc.cuh#L416-L417)에서 GPU당 두 개의 스트림을 만든다(`str1`, `str2`). 설계 의도(주석)는 **연산 커널과 P2P memcpy의 동시 실행**(comm-compute overlap). 다만 현 구현은 매 커널 호출 직후 `cudaDeviceSynchronize()`를 호출하므로(예: [ISPH_Calc.cuh:27](ISPH_Calc.cuh#L27)) 사실상 직렬화에 가깝다 — 향후 최적화 여지로 명시 가능.

### 2.2 2차 병렬화 — GPU 내부 입자-병렬

모든 SPH 연산은 **입자 단위 병렬**이다. 즉 한 입자 i에 대한 모든 이웃 j 합산을 한 CUDA 스레드가 처리한다 — *gather 패턴*.

표준 launch 형:

```c
dim3 b, t;
t.x = 128;
b.x = (num_part2 - 1) / t.x + 1;
KERNEL_xxx<<<b, t>>>(...);
cudaDeviceSynchronize();
```

- 블록 크기 128 = 4 warp는 점유율과 register pressure 사이의 합리적 절충. SPH 커널은 입자당 register 사용량이 높아(보통 60+ regs) 256으로 키우면 점유율이 떨어지는 경향이 있다.
- block 수는 `num_part2`에 비례. A100(108 SM, 64 warp/SM 한도)에서 약 7k 블록 이상이면 SM이 충분히 포화.
- **gather** 패턴이라 입자 i가 수정하는 메모리는 자기 자신뿐이므로 atomic이 거의 필요 없다(예외: 일부 KERNEL_clc_concentration 등의 누적). data race 없는 단순 병렬 구조.

### 2.3 한 타임스텝의 커널 시퀀스

[ISPH_Calc.cuh](ISPH_Calc.cuh) `SOPHIA_single_ISPH()`가 한 스텝마다 다음 순서로 커널을 호출한다(라인 번호 표기):

| 단계 | 커널 | 대상 | 라인 |
|---|---|---|---:|
| 1 | `KERNEL_set_ncell` | 적응 셀 반경 | [26](ISPH_Calc.cuh#L26) |
| 2 | `KERNEL_index_particle_to_cell` | 셀 해시 계산 | [35](ISPH_Calc.cuh#L35) |
| 3 | `cub::DeviceRadixSort::SortPairs` | (g_idx, p_idx) 정렬 | [39](ISPH_Calc.cuh#L39) |
| 4 | `KERNEL_reorder` | 정렬에 따라 P→SP 재배치 + g_str/g_end | [44](ISPH_Calc.cuh#L44) |
| 5 | `gradient_correction` + `KERNEL_clc_gradient_correction_3D` | KGC 행렬 | [150–152](ISPH_Calc.cuh#L150-L152) |
| 6 | `KERNEL_clc_prep3D` | Shepard filter, ρ_ref, p_type 전환 | [156](ISPH_Calc.cuh#L156) |
| 7 | `KERNEL_boundary3D` | no-slip wall vel 보정 | [162](ISPH_Calc.cuh#L162) |
| 8 | `KERNEL_MOST_boundary3D` | Monin-Obukhov 마찰속도 갱신 | [168](ISPH_Calc.cuh#L168) |
| 9 | `KERNEL_advection_force3D` | 점성/난류/중력 등 압력 외 힘 | [182](ISPH_Calc.cuh#L182) |
| 10 | `KERNEL_clc_projection` | predictor (u* 계산) | [192](ISPH_Calc.cuh#L192) |
| 11 | `IBM_force_interpolation3D` | IBM 마커가 유체 속도 보간 | [207](ISPH_Calc.cuh#L207) |
| 12 | `IBM_spreading_interpolation3D` | 마커→유체 보정력 spread | [211](ISPH_Calc.cuh#L211) |
| 13 | `IBM_corrector_fluid` | 유체에 fb 적용 | [220](ISPH_Calc.cuh#L220) |
| 14 | `KERNEL_zero_interior_velocity` | 건물 내부 유체 0화 | [224](ISPH_Calc.cuh#L224) |
| 15 | `KERNEL_PPE_prep3D` | PPE 행렬-벡터 RHS 구성 | [235](ISPH_Calc.cuh#L235) |
| 16 | `KERNEL_PPE_calc3D` | 압력 갱신(Jacobi step) | [237](ISPH_Calc.cuh#L237) |
| 17 | `KERNEL_Neumann_boundary3D` | wall pres Neumann | [240](ISPH_Calc.cuh#L240) |
| 18 | `KERNEL_pressureforce3D` | ∇p 힘 | [249](ISPH_Calc.cuh#L249) |
| 19 | `KERNEL_time_update_projection` | corrector — 위치/속도 최종 갱신 | [256](ISPH_Calc.cuh#L256) |
| 20 | `KERNEL_donothing_buffer3D`, `KERNEL_periodic_boundary3D`, `KERNEL_freeslip_boundary3D` | 경계 후처리 | [274–280](ISPH_Calc.cuh#L274-L280) |
| 21 | `KERNEL_compute_vorticity3D` | 보조량 | [284](ISPH_Calc.cuh#L284) |
| 22 | LDM block (간헐) | 분산입자 advect/deposit | [291–320](ISPH_Calc.cuh#L291-L320) |
| 23 | Scalar transport(BiMC, 간헐) | Eulerian 농도장 | [334–347](ISPH_Calc.cuh#L334-L347) |
| 24 | `cudaMemcpy` D2H + VTK 작성(간헐) | 출력 | [358–378](ISPH_Calc.cuh#L358-L378) |

매 50스텝마다 [ISPH_Calc.cuh:728-758](ISPH_Calc.cuh#L728-L758)에서 CUB의 `DeviceReduce::Max/Min`으로 진단량(최대 속도, 밀도, 힘, 온도 등)을 환산해 표준출력에 보고한다.

### 2.4 LDM 입자의 부분 빈도 처리

LDM 입자(`L_part1`)는 `freq_LDM` 스텝마다만 NNPS와 advection을 수행한다([ISPH_Calc.cuh:72](ISPH_Calc.cuh#L72), [291](ISPH_Calc.cuh#L291)). 시간 단계는 `freq_LDM*dt`로 키운다. 이는

- LDM 통계 수렴에 SPH만큼 잦은 갱신이 불필요하고,
- LDM 농도 계산(`KERNEL_clc_concentration`)이 자체 NNPS 비용을 갖기 때문이다.

---

## 3. 데이터 레이아웃 (논문 2.3.1용)

### 3.1 3-구조체 분리 (hot/cold/scratch)

[class_Cuda_Particle_Array.cuh](class_Cuda_Particle_Array.cuh)의 핵심 설계는 입자 정보를 **세 구조체로 쪼개서** 매 스텝의 메모리 트래픽을 줄이는 것이다.

| 구조체 | 의미 | 매 스텝 | 사이즈 추정* |
|---|---|---|---:|
| **`part1`** | "hot" — 매 스텝 read+write되는 1차 변수 (x,y,z, ux,uy,uz, m, h, temp, pres, rho, k_turb, e_turb, vis_t, fbx/y/z, fmx/y/z, PPE1, PPE2, fp_x/y/z, fv_x/y/z, fe_x/y/z, ncell, source) | R/W 빈번 | ~ 200 B |
| **`part2`** | "cold" — 초기값(레퍼런스). x0/y0/z0, ux0/uy0/uz0, rho_ref, temp0, enthalpy0 | predictor/corrector에서만 R | ~ 36 B |
| **`part3`** | "scratch" — 매 스텝 0으로 클리어되는 중간량(drho, dconcn, ftotalx/y/z, vorticity, KGC 3×3 행렬 Cm[3][3], inv_cm_** 6성분, 색깔필드 grad 등) | 매 스텝 시작 시 `cudaMemset(0)` | ~ 120 B |

\* 정렬 패딩 무시한 추정. 실제로는 컴파일러 정렬에 따라 다소 증가.

이 분리의 효과:

- `cudaMemset(dev_P3, 0, sizeof(part3)*N)`로 **scratch 영역만** 초기화 → 매 스텝 hot/cold는 그대로 보존.
- predictor/corrector에서 P2의 초기위치/초기속도는 read-only, P1은 갱신 대상 → register/캐시 압력이 다른 두 흐름을 분리.
- L1/L2 캐시 지역성: 같은 스텝 내에서 part1만 집중적으로 액세스되는 커널이 많기 때문에, 큰 part 구조체 하나로 합칠 때보다 **유효 cache line 활용률이 높다**.

### 3.2 Array of Structures (AoS) 채택 — 그리고 그 trade-off

`part1*` 등은 typed AoS이다. SPH gather 커널의 각 thread는 자신의 입자 i를 한 번 읽고, 이웃 j의 모든 필드를 *연속 위치에서* 가져온다(같은 캐시라인 안에 j의 x,y,z,ux,uy,uz가 함께 있는 형태). 이는

- ✅ **이웃 j 액세스의 L1 hit율 양호** — 한 입자의 다수 필드를 연달아 사용
- ❌ 같은 필드(예: P1[].x)만 이웃들에 대해 broadcast로 쓰는 패턴에서는 **bandwidth 효율 ~30%만 활용** (필요 4B를 위해 cache line 128B를 끌고 옴)

본 코드는 한 커널 안에서 한 이웃의 다수 필드를 함께 사용하는 패턴이 많아 AoS가 합리적 선택이지만, **scalar transport / vorticity / 단일 양 reduction 커널**은 SoA가 더 빠를 수 있다 — 향후 mixed-layout 또는 KGC 6-성분 분리 등이 잠재적 최적화 후보.

### 3.3 두 번 사본: P1/SP1 더블 버퍼

[ISPH_Calc.cuh:512-513](ISPH_Calc.cuh#L512-L513) — 모든 입자 배열은 **`dev_P1`(이전 스텝 결과)**과 **`dev_SP1`(정렬된 작업본)** 두 벌이 있다.

스텝 시작 시:

- 0스텝([ISPH_Calc.cuh:35-51](ISPH_Calc.cuh#L35-L51)): cell index → sort → reorder를 **수행**한 뒤, `cudaMemcpy(dev_P1, dev_SP1, ..., D2D)`로 양쪽을 동기화. 입력파일 순서를 잃은 채로 정렬된 메모리 배열을 사용.
- 이후 스텝([ISPH_Calc.cuh:61](ISPH_Calc.cuh#L61)): 이미 정렬되어 있다고 가정하고 `dev_SP1 ← dev_P1` 단순 복사. (입자가 셀 사이를 거의 이동하지 않으므로 매 스텝 sort를 강요하지 않는 설계 — 단, 이는 large dt나 빠른 입자에서 NNPS 정확도가 떨어질 수 있는 trade-off다.)

이 더블버퍼 구조가 SPH 커널 작성을 단순하게 만든다: **읽기는 SP1에서, 쓰기는 P1으로** — RAW hazard 없는 깔끔한 데이터 흐름.

### 3.4 입자 종류와 분기

[Parameters.cuh](Parameters.cuh#L9-L13)에 정의된 `p_type`:

| `p_type` | 종류 | 처리 |
|---|---|---|
| 0 | BOUNDARY (지면, 고체벽) | velocity=0, density 고정 |
| 1 | FLUID (공기) | 정상 풀 SPH |
| 2 | AIR (2상에 사용; 본 케이스 거의 미사용) | — |
| 9 | MOVING (이동 벽) | predictor에서 강제 |
| ≥1000 | IBM 마커 (건물 표면) | 별도 IBM 3-step |

각 커널 진입부에서 `if(P1[i].p_type ...)`로 early return을 통해 분기를 처리. 분기 자체는 warp divergence를 일으키나, **공간 정렬 후에는 같은 종류 입자가 같은 셀에 모이는 경향**이 있어 동일 warp 안 분기 비율이 낮다 — sort-and-reorder가 일종의 분기 압축 효과를 함께 가져온다.

`buffer_type`은 추가적인 상태 비트:

- 0: active (정상)
- 1–5: open boundary buffer (inlet/outlet/left/right/ceiling)
- 99: 건물 내부에 갇힌 유체 입자 (BFS 외부 마킹 후 반대로 식별; PPE/시간갱신에서 skip)
- (LDM의 `L_part1.buffer_type`은 0=active / 1=deposited로 의미 다름)

### 3.5 솔버 파라미터: `vii[64]`/`vif[32]`

[function_init.cuh:1-83](function_init.cuh#L1-L83) — 정수/실수 솔버 파라미터를 두 개의 짧은 배열로 모아두고, **매크로 alias**로 의미 있는 이름을 부여한다(`#define dt vif[2]` 등). GPU 측에는 `__constant__ int_t k_vii[vii_size]`/`__constant__ Real k_vif[vif_size]`로 한 번에 복사([SOPHIA_gpu.cu:48-49](SOPHIA_gpu.cu#L48-L49), [ISPH_Calc.cuh:557-558](ISPH_Calc.cuh#L557-L558)).

이점:

- 모든 커널이 broadcast로 0-cycle latency에 가깝게 솔버 옵션 액세스 (constant cache).
- 신규 옵션 추가 시 한 슬롯 +`#define` 한 줄.
- 단점: 슬롯 번호 직접 관리(현재 vii[48]까지 사용) — 이름 충돌 위험.

물성 테이블(T, h, k, cp, vis)도 [SOPHIA_gpu.cu:52-58](SOPHIA_gpu.cu#L52-L58)에서 같은 방식으로 `__constant__`에 적재.

---

## 4. 이웃 탐색 (논문 2.3.2용)

### 4.1 균일 셀 격자

[ISPH.cuh:96-104](ISPH.cuh#L96-L104):

```c
Real h0 = HP1[0].h_ref;
dcell = cell_reduction_factor * kappa * h0 / float(ncell_init);   // 1.1 * κh / 2
NI = (int)((x_max-x_min)/dcell)+1;
NJ = (int)((y_max-y_min)/dcell)+1;
NK = (int)((z_max-z_min)/dcell)+1;
```

- 셀 변길이 `dcell`은 **반경 κh보다 작게** 잡혀 있다(`ncell_init=2`이므로 `dcell ≈ 0.55·κh`). 따라서 한 입자의 영향권을 덮으려면 ±2 셀(중심 포함 5×5×5=125 셀)을 훑어야 한다 — `Parameters.cuh:103`의 `ncell_init=2`.
- **이점**: dcell이 작을수록 셀당 입자 수가 적어 stencil 안에서 실제 거리 검사로 걸러지는 잉여 입자가 적다. 셀당 평균 입자 수 약 8~20 정도로 균형.
- 셀 수가 (Nx·Ny·Nz)+1로 1 슬롯 추가되는데(`clc_num_cells`, [functions_NNPS.cuh:294-306](functions_NNPS.cuh#L294-L306)) **마지막 슬롯은 dummy**(`i_type==3`) 전용. 즉 dummy 입자들은 모두 셀 인덱스 `num_cells`로 흘려보내져 정렬 후 배열의 끝에 모인다.

### 4.2 적응 stencil 반경 `ncell`

[functions_NNPS.cuh:765-777](functions_NNPS.cuh#L765-L777):

```c
__global__ void KERNEL_set_ncell(part1*P1){
    ...
    P1[i].ncell = ceil(float(ncell_init) / float(h0/h));
}
```

가변 smoothing length를 사용하는 입자(APR/multi-resolution)는 자신의 h가 h_ref와 다를 수 있다. 그 경우 **stencil 반경 ncell을 입자별로 다르게** 설정. 이후 모든 SPH 루프는 자기 자신의 `P1[i].ncell`을 그대로 사용한다(`for(z=-ncell; z<=ncell; z++)`).

### 4.3 NNPS 4단계: hash → sort → mark → reorder

#### 4.3.1 Hash — 셀 인덱스 계산

[functions_NNPS.cuh:671-699](functions_NNPS.cuh#L671-L699) `KERNEL_index_particle_to_cell`:

```c
icell = min(floor((P1[idx].x - k_x_min)/(k_x_max - k_x_min) * k_NI), k_NI-1);
g_idx[idx] = idx_cell(icell, jcell, kcell);
p_idx[idx] = idx;     // identity permutation
if(i_type==3){ g_idx[idx]=k_num_cells; }   // dummy → bucket 끝
```

`idx_cell(i,j,k) = i + NI·j + NI·NJ·k` ([functions_NNPS.cuh:317-324](functions_NNPS.cuh#L317-L324)) — z-축이 가장 느리게 변하는 row-major 형. 시뮬레이션 영역이 ABL인 경우 입자가 z 방향으로 거의 안 움직이므로 같은 z-슬랩에 있는 이웃들이 메모리에서 인접해진다 → ∇²p 등 **수직 stencil의 cache locality가 좋아진다**.

#### 4.3.2 Sort — CUB Radix Sort (key-value)

[ISPH_Calc.cuh:39](ISPH_Calc.cuh#L39):

```c
cub::DeviceRadixSort::SortPairs(
    dev_sort_storage, *sort_storage_bytes,
    g_idx_in, g_idx,             // keys (cell indices)
    p_idx_in, p_idx,             // values (original particle indices)
    num_part2);
```

- 키: 32-bit unsigned 셀 인덱스 (`int_t = int`, 양수만 들어감).
- 값: 정렬 전 슬롯 번호. 정렬 후 `p_idx[k]`는 "k번째 정렬 슬롯에 들어갈 입력 슬롯 번호"가 된다.
- 임시 storage는 [ISPH_Calc.cuh:659-665](ISPH_Calc.cuh#L659-L665)에서 dummy call로 크기를 받아낸 뒤 `cudaMalloc` (CUB의 표준 패턴).
- A100에서 ~수천만 키의 SortPairs는 ~수 ms 수준. 매 스텝 호출되는 가장 비싼 단일 연산 중 하나.

#### 4.3.3 Mark — 셀 시작/끝 인덱스 (CSR-like)

[functions_NNPS.cuh:573-604](functions_NNPS.cuh#L573-L604) `KERNEL_reorder`는 정렬 후 두 가지 일을 동시에 한다.

1. 같은 hash가 시작·끝나는 위치를 표시:
   ```c
   if(idx==0 || hash != sharedHash[threadIdx.x]){
       g_str[hash] = idx;
       if(idx>0) g_end[sharedHash[threadIdx.x]] = idx;
   }
   if(idx == k_num_part2-1) g_end[hash] = idx+1;
   ```
   - `g_str[c]` = 셀 c의 첫 입자 슬롯, `g_end[c]` = 마지막 입자 슬롯+1
   - 빈 셀은 `g_str[c]==cu_memset(0xFFFFFFFF)`로 표시되어 [functions_PREP.cuh:84](functions_PREP.cuh#L84) 등에서 `if(g_str[k]!=cu_memset)`로 skip.
2. 데이터 재배치: `SP1[idx] = P1[p_idx[idx]];` — *gather copy*. 이로써 같은 셀에 속한 입자들이 메모리에서 연속해진다.

`__shared__ int sharedHash[t.x+1]` 사용은 인접 thread가 자기 hash와 직전 thread의 hash를 비교할 때 *전역 메모리 재읽기*를 피하기 위함. shared memory size = `sizeof(int)*(t.x+1)`이 [ISPH_Calc.cuh:20](ISPH_Calc.cuh#L20)에서 launch 시 동적으로 잡힌다.

#### 4.3.4 Reorder — 정렬된 배열로 활용

위 단계 직후 `cudaMemcpy(dev_P1, dev_SP1, ..., D2D)`로 P1을 정렬본과 일치시킨다([ISPH_Calc.cuh:51](ISPH_Calc.cuh#L51)). 이후 스텝부터는 입자가 거의 이동 안 했다고 가정하고 sort를 생략 — 단, 0스텝 외에는 [ISPH_Calc.cuh:29-69](ISPH_Calc.cuh#L29-L69)에서 보듯 **단순 복사만 수행**. (`count==0` 분기에서만 본격 NNPS 재계산.)

> 노트: 현재 구현은 매 스텝 NNPS를 다시 호출하지 **않는다**. 이는 ABL 시뮬레이션에서 입자 변위가 cell 크기 대비 작다는 가정 위에 있다. 정확성을 우려한다면 N 스텝마다 강제 재정렬을 추가하는 것이 안전책이며, 향후 검증 항목으로 명시 가능.

### 4.4 Stencil 루프 (SPH 합산 패턴)

전형적 SPH 가우터 합산은 [functions_PPE.cuh:51-71](functions_PPE.cuh#L51-L71) 등에서 다음 형태로 등장:

```c
int_t ncell = P1[i].ncell;
for(int_t z=-ncell; z<=ncell; z++){
  for(int_t y=-ncell; y<=ncell; y++){
    for(int_t x=-ncell; x<=ncell; x++){
      int_t k = idx_cell(icell+x, jcell+y, kcell+z);
      if(out-of-domain) continue;
      if(g_str[k] != cu_memset){
        for(int_t j=g_str[k]; j<g_end[k]; j++){
          // 거리 검사 + W, ∇W 평가 + 누적
        }
      }
    }
  }
}
```

특징:

- 외부 3중 루프: 5×5×5=125 셀 (3D, ncell=2). 통상적으로 8할은 `g_str==cu_memset`(빈 셀)이거나 거리 검사에서 컷.
- 내부 j-루프: **해당 셀에 모인 입자만 순회**. CSR-like 표현 덕분에 zero-padding 없이 메모리 연속.
- 거리 검사: `if(tdist < search_range)` — `search_range = κ · max(hi, hj)`. multi-resolution을 가정한 fmax.
- 모든 변수는 register/shared cache에 머물고, 누적 결과는 마지막에 `P3[i].ftotalx = ...` 한 번만 global write.

이 형태의 메모리 액세스는

- ✅ 같은 셀의 입자 j가 메모리에서 연속 → **j-루프의 coalesced load**
- ✅ 인접 thread (인접 입자 i)가 **거의 같은 이웃 셀**을 본다 → **L2 cache reuse** 매우 높음
- ❌ warp 내 분기 발산: 인접 입자가 같은 ncell을 가질 때만 효율적. APR/h_change 사용 시 ncell이 발산하면 동기화 손실.

### 4.5 LDM의 별도 NNPS

[ISPH_Calc.cuh:78-97](ISPH_Calc.cuh#L78-L97) — LDM 입자는 **자체적인 g_str_LDM/g_end_LDM 한 쌍**을 가진다(셀 격자는 동일하게 공유). 이는 SPH 입자 수와 LDM 입자 수가 같지 않기 때문이며, 또 LDM 농도 계산에서 LDM-LDM 이웃이 필요하기 때문이다. SPH 입자 정보가 필요할 때(예: 풍속 보간 — `KERNEL_time_update_LDM`)는 SPH의 g_str/g_end를 그대로 쓰고 LDM 입자 위치만 입력으로 넣는다.

---

## 5. 메모리 관리 (논문 2.3.2용)

### 5.1 메모리 계층별 사용

| 영역 | 활용 | 코드 위치 |
|---|---|---|
| **Global** | 입자 배열(P1/SP1, P2/SP2, P3, LP1/SLP1, LP2/SLP2), NNPS 인덱스, sort temp storage, max-reduce temp | `cudaMalloc` 일괄 (ISPH_Calc.cuh) |
| **Constant (`__constant__`)** | `k_vii`, `k_vif`, `k_Tab_T/h/k/cp/vis`, table 인덱스 | 64+32 + 5×50 = ~620 슬롯 |
| **Shared (`__shared__`)** | `KERNEL_reorder`의 `sharedHash[]` 만 | dynamic, t.x+1 ints |
| **Texture/Surface** | 사용 안 함 | (잠재적 최적화 항목) |
| **Pinned host** | 멀티 GPU 교환용 `send_*`, `recv_*` 버퍼 | (멀티 GPU 모드에서만) |

> Shared memory를 거의 쓰지 않는 이유: SPH의 stencil이 데이터 의존적(이웃 입자 수가 동적)이고 cell 단위 협업이 깔끔하지 않다. 대신 L1/L2 cache에 의존 — A100에서 SM당 192KB SRAM이 L1+shared로 자동 분할되므로 **자동 캐시 활용이 충분**한 형태.

### 5.2 글로벌 메모리 풋프린트 추정

`num_part2`개 입자(보통 num_part의 1.0~1.2배 + 버퍼) 기준:

- `dev_P1` + `dev_SP1` ≈ 2 × 200 B/part × N
- `dev_P2` + `dev_SP2` ≈ 2 × 36 B × N
- `dev_SP3` ≈ 120 B × N
- LDM 4벌: 각 `L_part1` ≈ 180 B × N_LDM
- NNPS 인덱스(g_idx, p_idx, g_idx_in, p_idx_in) = 4·int·N
- g_str / g_end = 2·int·num_cells
- CUB sort temp = O(N) (수 ×N bytes)

도시 풍환경 케이스(N ≈ 5–20 M 입자, num_cells ≈ 수십 M)에서 **A100 80GB 1장에 무리 없이 적재**. 다만 num_cells = NI·NJ·NK가 영역 종횡비에 따라 매우 커질 수 있으므로(예: ABL 시뮬레이션에서 NI×NJ×NK가 수억) **num_cells 메모리가 입자 배열보다 커지는 케이스**도 발생 — 이때는 (i,j,k)→hash 매핑을 z-curve(Morton code) 등으로 압축하는 것이 향후 최적화 후보.

### 5.3 매 스텝 메모리 트래픽

매 스텝 transferred = (part1 R+W + part2 R + part3 R+W + g_idx/p_idx R+W + sort temp …) × N. SPH 커널이 ~10–20개 호출되며 각각 part1을 한번 이상 훑으므로, 한 스텝당 글로벌 전송이 대략 10·sizeof(part1)·N에 비례. A100의 HBM2e 1.55 TB/s를 활용하기 위해서는 **각 커널이 가능한 한 많은 일을 한 번의 read에서 끝내야** 한다. 현재 코드의 `KERNEL_clc_prep3D`처럼 여러 보조량을 한 커널에서 같이 계산하는 패턴은 이를 의도한 것이다.

### 5.4 정렬·축소(reduction) 보조 저장

[ISPH_Calc.cuh:651-665](ISPH_Calc.cuh#L651-L665) — CUB의 표준 2-pass 패턴:

```c
size_t bytes = 0;
cub::DeviceRadixSort::SortPairs(NULL, bytes, ...);   // pass 1: query
cudaMalloc(&dev_sort_storage, bytes);                // allocate once
// 매 스텝 같은 storage 재사용
cub::DeviceRadixSort::SortPairs(dev_sort_storage, bytes, ...);  // pass 2: do
```

매 스텝 alloc/free를 피해 메모리 단편화를 줄이고 launch overhead를 절약한다. `DeviceReduce::Max/Min`도 동일 패턴([ISPH_Calc.cuh:736-742](ISPH_Calc.cuh#L736-L742)).

### 5.5 명시적 0-초기화 정책

매 스텝 시작 시 [ISPH_Calc.cuh:48,67](ISPH_Calc.cuh#L48):

```c
cudaMemset(dev_P3, 0.0, sizeof(part3)*num_part2);
```

`part3`의 모든 누적 변수(ftotalx, drho, KGC 행렬 등)가 0이 됨이 보장. 누적 커널은 +=로 안전하게 작성.

`g_str`은 `cu_memset = 0xFFFFFFFF`로 채워져 빈 셀 표시([ISPH_Calc.cuh:31](ISPH_Calc.cuh#L31)).

---

## 6. 침지 경계법(IBM)의 GPU 구현

[functions_IBM.cuh](functions_IBM.cuh)는 본 코드의 도시 형상(건물) 처리 핵심이다. **마커 단일층** 표현 + **3-단계 알고리즘**을 GPU에서 다음과 같이 매핑.

| 단계 | 커널 | 병렬 단위 | 누적 방향 |
|---|---|---|---|
| ① 보간 | `IBM_force_interpolation3D` | IBM 마커 입자별(p_type ≥ 1000) | 자기 자신에만 write |
| ② 확산 | `IBM_spreading_interpolation3D` | 유체 입자별 | 자기 자신에만 write (gather) |
| ③ 보정 | `IBM_corrector_fluid` | 유체 입자별 | 속도 덮어쓰기 |

핵심은 **②가 push가 아니라 gather**라는 점이다. naive 구현은 IBM 마커가 주변 유체에 force를 *push*(atomicAdd)할 것이지만, 본 코드는 유체 입자가 IBM 마커를 *gather*하므로 atomic이 필요 없다. 수학적으로 같은 결과지만 GPU 효율은 압도적으로 낫다.

---

## 7. 압력 Poisson 방정식(PPE) 솔버

[functions_PPE.cuh](functions_PPE.cuh):

- `KERNEL_PPE_prep3D` — 우변(divergence of u*) + system 행렬 `Aij`, `AijPij` 합산을 **행렬을 명시적으로 저장하지 않고** 입자별로 즉석 계산.
- `KERNEL_PPE_calc3D` — Jacobi-style update: `pres_i ← bi/Aij` 또는 유사한 한 step.
- [ISPH_Calc.cuh:232](ISPH_Calc.cuh#L232) `for(int_t it=0; it<=0; it++)` — 현재는 1회 반복 (외부 반복 가능 구조만 유지). 더 큰 dt나 강한 압력 변동에서 it를 늘려야 할 수 있음.

행렬-자유(matrix-free) 형태 → **메모리 풋프린트 = O(N)**. CSR 행렬을 명시 저장하지 않아 NNPS 비용이 매 PPE 반복마다 발생하지만(이웃 다시 훑기), 메모리 측면에서 큰 이득.

---

## 8. 시간 적분

[functions_TIME_ISPH.cuh](functions_TIME_ISPH.cuh):

- `KERNEL_clc_projection` (predictor): `u* = u^n + Δt · (f_advection + f_visc + f_turb + f_geo)`
- IBM 단계에서 `u** = u* + Δt · f_IBM/ρ`
- PPE 풀고 ∇p 계산
- `KERNEL_time_update_projection` (corrector): `u^{n+1} = u** − Δt · ∇p / ρ`, 위치 갱신

특이점:

- ABL 정상상태를 위해 [functions_TIME_ISPH.cuh:49](functions_TIME_ISPH.cuh#L49)에서 **지오트로피컬 압력 구배** `f_geo = u*²/ABL_h`를 가산.
- predictor와 corrector 모두 입자별 독립 → 100% 병렬, 분기 발산 거의 없음.

---

## 9. 출력 — VTK Binary

[functions_OUTPUT.cuh](functions_OUTPUT.cuh) (3.4k 줄)는 ParaView가 직접 읽는 VTK POLYDATA 바이너리를 작성한다. 매 `freq_output` 스텝마다 `cudaMemcpy(file_P1, dev_P1, ..., D2H)`로 호스트로 가져와 파일로 출력 — D2H는 보통 sub-second 수준이지만 **빈도가 높을수록 PCIe 대역폭이 병목**이 된다.

CASE_ID 분기에 따라 fluid, building, vorticity, LDM 등 별도 파일로 분리 출력.

---

## 10. 수치적·구현적 trade-off 정리 (논문 토의용)

| 결정 | 이득 | 비용/한계 |
|---|---|---|
| 단정도 float | 처리량 ×2, 메모리 절반 | ABL 누적 통계의 round-off error 가능 |
| 단일 변환단위 컴파일 | 인라이닝, 단순 빌드 | 증분 빌드 불가 |
| AoS 3-구조체 분리 | hot 데이터 cache locality | 같은 필드 broadcast 시 효율 감소 |
| 매 스텝 NNPS 생략(0스텝 외) | 정렬 비용 절약 | 큰 dt에서 hash가 stale될 위험 |
| Cell 변길이 < κh | 거리 검사 컷오프율↑ | 5×5×5=125 셀 stencil 필요 |
| Matrix-free PPE | 메모리 O(N) | 매 PPE iter마다 NNPS 비용 |
| pthread+barrier 멀티 GPU | 단순 동기 모델 | barrier가 가장 느린 GPU에 묶임 |
| `cudaDeviceSynchronize` 매 커널 | 디버깅·정확성 보장 | 잠재 동시성(stream)이 직렬화됨 |
| Constant memory에 옵션 적재 | 0-cycle broadcast | 슬롯 수동 관리 (현재 vii[48]까지) |

---

## 11. 논문 작성 시 활용 팁 (메모)

1. **2.3.1**: 본 문서 §2 (병렬화 전략) + §3 (데이터 레이아웃) 그대로 본문화. Figure는 (a) 멀티 GPU 슬랩 분할 + i_type, (b) 3-구조체 hot/cold/scratch 다이어그램, (c) 한 스텝 커널 시퀀스 박스도가 효과적.
2. **2.3.2**: 본 문서 §4 (NNPS) + §5 (메모리). Figure는 (a) cell hash → CUB sort → reorder의 4단계 다이어그램, (b) g_str/g_end CSR-like 그림, (c) stencil 루프 의사코드.
3. **수식과 코드의 연결**: 코드 라인 인용은 부록에 두고 본문은 의사코드(pseudocode)로 추상화. 본 문서의 §4.4 stencil 루프 그대로 리스팅으로 사용 가능.
4. **벤치마크 절(향후)**: 만약 §2.3.3을 추가한다면, A100 1장 vs 4장 weak/strong scaling, NNPS/PPE/IBM의 단계별 wall-clock 비율, PCIe vs NVLink P2P 비교가 표준 지표.
5. **선행 연구와 차별점**: ① 단일 TU 빌드, ② 3-구조체 hot/cold 분리, ③ matrix-free PPE+gather IBM 조합은 학위논문 단계에서 강조할 만한 구현 결정.

---

## 부록 A. 주요 수치 상수 (Parameters.cuh 발췌)

- `kappa = vif[0]` — 커널 컷오프 배수 (대개 2.0; runtime)
- `ncell_init = 2` — stencil 반경 (셀 단위)
- `IBM_length = 0.5` — IBM 커널 반경 배수
- `Pressure_length = 1.0`
- `MAX_ITER = 200`, `TOLERANCE = 1e-5` — MOST 반복 솔버
- `Max_GPU = 10` — 동시 사용 가능 GPU 수 상한
- `vii_size = 64`, `vif_size = 32` — solver 옵션 슬롯 수
- `table_size = 50` — 물성 테이블 크기

## 부록 B. CASE_ID 약속 (case_config.cuh)

| 범위 | 의미 |
|---|---|
| 100–199 | Case A 계열 |
| 200–209 | Case C scaled |
| 210–219 | 안정도 비교 (unstable/neutral/stable) |
| 300–399 | Case C 정규 |
| 1000–1099 | Leipzig 외 ABL |
| 2001–2003 | Scalar transport 검증 |
| 3001 | MOST 검증 |

상세는 [case_config.cuh](case_config.cuh).

---

*본 문서는 학위논문 작성 시점(2026-04-26)의 코드 상태(branch: LDM-couple)를 기준으로 작성되었다. 코드 변경 시 본 문서도 함께 갱신할 것.*
