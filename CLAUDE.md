# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**SOPHIA_gpu** — GPU-accelerated Incompressible Smoothed Particle Hydrodynamics (ISPH) solver for atmospheric boundary layer flow over urban terrain. Developed at Seoul National University (ESLAB). Written in CUDA C++ with pthreads for multi-GPU support.

Primary applications: urban wind flow, pollutant/scalar dispersion (Eulerian + Lagrangian), building aerodynamics.

## Build & Run

```bash
# Build (single command, no dependencies beyond CUDA toolkit + bundled CUB 1.8.0)
make

# Run (argument = number of GPUs)
./SOPHIA_gpu 1

# Target architecture: sm_80 (A100 / RTX 30-series)
# To change GPU arch, edit -arch=sm_80 in Makefile
```

There are no tests, linters, or CI pipelines. Validation is done via Python analysis scripts against analytical solutions or experimental data.

## Architecture

### Compilation Model

All code is in `.cuh` headers `#include`'d into a single translation unit (`SOPHIA_gpu.cu`). There is no separate compilation — the entire solver compiles as one nvcc invocation. Include order matters (defined in SOPHIA_gpu.cu lines 110-130).

### Execution Flow

```
main() → read_solv_input() → read_table() → ISPH()
  └→ ISPH_Calc() [per-GPU pthread, runs the main time loop]
       Each timestep:
         1. NNPS (Neighbor search: cell index → RadixSort → reorder)
         2. Gradient correction (KGC matrix)
         3. Preparation (strain rate, turbulence, filter)
         4. MOST boundary (Monin-Obukhov surface layer)
         5. Predictor (explicit time advance)
         6. IBM pass (force interpolation → spreading → fluid correction)
         7. PPE (Pressure Poisson Equation — iterative solver)
         8. Pressure force
         9. Corrector (pressure correction → velocity update)
        10. Boundary conditions (Neumann, open boundary)
        11. Diffusion (scalar transport, energy)
        12. LDM (Lagrangian dispersion — sort, advect, deposit)
        13. Output (VTK binary)
```

### Key Data Structures

Three-struct particle design for GPU memory efficiency:
- **part1** — primary fields: position, velocity, mass, density, pressure, forces, turbulence (hot data, read/written every step)
- **part2** — reference/initial values: x0/y0/z0, ux0/uy0/uz0, rho_ref (used in predictor-corrector)
- **part3** — computed intermediates: correction matrix Cm[3][3], total forces, gradients (zeroed each step)

LDM has separate structs: `L_part1`, `L_part2`, `L_part3`.

Defined in `class_Cuda_Particle_Array.cuh`. Types: `Real` = float, `int_t` = int (see `Variable_Type.cuh`).

### Particle Types

- `p_type == 0`: Boundary (ground, walls)
- `p_type == 1`: Fluid (air)
- `p_type == 2`: Air (secondary phase, rarely used)
- `p_type >= 1000`: IBM markers (building surfaces, single-layer)

### Multi-GPU

Domain decomposition along x-axis. Each GPU runs `ISPH_Calc()` on a pthread. Synchronization via `pthread_barrier_wait()`. Inter-GPU particle exchange uses pinned host buffers with peer-access (`cudaDeviceEnablePeerAccess`).

### Configuration

**Compile-time** (`Parameters.cuh`): physical constants, IBM_length, MOST parameters, CASE_NAME, input/output paths. Changing cases requires editing `CASE_NAME` and recompiling.

**Runtime** (`input/<CASE_NAME>/<CASE_NAME>_solv.txt`): solver type, kernel, turbulence model, time step, boundary conditions, output variables.

**Case geometry** (`case_config.cuh`): domain extents (L1-L4, L_left/right/ceiling), CASE_ID-based `#if` blocks. Each CASE_ID maps to a domain configuration.

**Particle input** (`input/<CASE_NAME>/<CASE_NAME>.txt`): tab-separated, first row is column indices (1=x, 2=y, 3=z, 7=mass, 8=p_type, 9=h, 12=rho, etc.). Columns defined in `function_init.cuh` (inp_x through inp_m0).

### IBM (Immersed Boundary Method)

Building surfaces represented as single-layer particles (p_type=1000). Three-step process:
1. **Force interpolation** (`IBM_force_interpolation3D`): IBM particle interpolates fluid velocity via SPH kernel → computes velocity correction fb = u_desired - u_interpolated
2. **Spreading** (`IBM_spreading_interpolation3D`): fluid particles accumulate corrections from nearby IBM particles via `Σ fb_j * W(r, h_IBM) * m_j/ρ_j`
3. **Fluid correction** (`IBM_corrector_fluid`): applies accumulated fb to fluid velocity

IBM kernel size controlled by `IBM_length * h` (Parameters.cuh). Search range = `kappa * IBM_length * h`.

**Known issue**: IBM particle mass must be consistent with actual particle spacing. Sub-grid IBM particles (from STL mesh conversion) with mass set for coarser spacing cause over-counting in the spreading step.

**Known issue (PPE iteration count is fixed at 1)**: [ISPH_Calc.cuh:250](ISPH_Calc.cuh#L250) hardcodes `for(int_t it=0; it<=0; it++)` — only one Jacobi sweep per timestep. The solv-file fields `minimum-iteration` / `maximum-iteration` are parsed (`minIteration`, `maxIteration`) but **not** used in the PPE loop. Confirmed empirically that increasing `it<n` does **not** resolve wake/bluff-body velocity spikes, so don't bother suggesting "raise PPE iteration count" as a fix for those.

### Output

VTK binary (POLYDATA) files written to `./<CASE_NAME>/`. Readable by ParaView.
- `fluid_<step>stp.vtk` — fluid particles
- `building_<step>stp.vtk` — IBM/building particles
- `LDM_<step>stp.vtk` — Lagrangian dispersion particles

### Dry Deposition & Resuspension (LDM)

Implemented in `functions_DD.cuh`, called from `ISPH_Calc.cuh` once per `freq_LDM` steps after `KERNEL_time_update_LDM`.

**Particle state** tracked via `L_part1.buffer_type` (added field):
- `buffer_type == 0`: active — subject to Lagrangian advection + deposition check
- `buffer_type == 1`: deposited — immobile, subject to resuspension check only

Deposited particles are skipped in `KERNEL_clc_predictor_Dispersion` and `KERNEL_time_update_LDM` (early return).

**Deposition velocity** (Wesely 1989, Zhang et al. 2001, Slinn 1982):
```
Vd = 1/(Ra + Rb) + Vs

Ra = Phi_h(OB_length, z, z0) / (κ · u*)   [aerodynamic resistance]
Rb = 1 / (u* · ε₀ · (E_B + E_IM))         [viscous sublayer resistance]
   E_B  = Sc^(-2/3)                         Brownian diffusion efficiency
   E_IM = [St/(St+0.5)]²                    inertial impaction efficiency
Vs = 8.9×10⁻⁵ m/s                          gravitational settling (Cs-137, 1 μm, constant)
```

`u*` and `OB_length` are read from global `__device__` variables `us_update` and `OB_length_global` (both updated by `KERNEL_MOST_boundary3D` at the domain reference point).

**Stochastic deposition** (Stohl et al. 2005 / FLEXPART):
```
P_dep = 1 - exp(-Vd · Δt / dz)
```
Triggered when particle height above ground `dz < DD_h_dep` (default 2 m).

**Stochastic resuspension** (Biasi et al. 2001, Loosmore 2003):
```
λ_res = C_res · max(0, u*² - u*_t²)   [s⁻¹]
P_res = 1 - exp(-λ_res · Δt)
```

**Tunable `#define` constants in `functions_DD.cuh`**:

| Constant | Default | Description |
|----------|---------|-------------|
| `DD_h_dep` | `2.0` m | Height above ground to trigger deposition check |
| `DD_ustar_thresh` | `0.2` m/s | Threshold friction velocity for resuspension |
| `DD_C_res` | `1e-3` s⁻¹/(m/s)² | Resuspension rate coefficient |
| `DD_dp` | `1e-6` m | Particle diameter (Cs-137) |
| `DD_tau_p` | `9e-6` s | Particle relaxation time |
| `DD_Vs_const` | `8.9e-5` m/s | Gravitational settling velocity |
| `DD_eps0` | `3.0` | Viscous sublayer empirical constant |

**Resuspension rate examples** (u*_t = 0.2 m/s):
- u* = 0.3 m/s → λ ≈ 5×10⁻⁵ s⁻¹ (half-life ~4 h)
- u* = 0.5 m/s → λ ≈ 2×10⁻⁴ s⁻¹ (half-life ~55 min)
- u* = 1.0 m/s → λ ≈ 1×10⁻³ s⁻¹ (half-life ~12 min)

**Note**: RNG for deposition/resuspension uses the same `LP1[i].seed` chain as the Langevin step but with `curand_init` offset +17 for statistical independence.

## Key Files

| File | Purpose |
|------|---------|
| `SOPHIA_gpu.cu` | Entry point, global declarations, includes |
| `ISPH_Calc.cuh` | Main time-stepping loop, kernel launch sequence |
| `ISPH.cuh` | Memory allocation, multi-GPU setup, I/O orchestration |
| `function_init.cuh` | Input file parsing, solver parameter mapping |
| `case_config.cuh` | Per-case domain geometry definitions |
| `Parameters.cuh` | Compile-time constants, case name, file paths |
| `functions_IBM.cuh` | Immersed boundary method kernels |
| `functions_PPE.cuh` | Pressure Poisson equation solver |
| `functions_TIME_ISPH.cuh` | Predictor-corrector time integration |
| `functions_MOST.cuh` | Monin-Obukhov similarity theory (ABL) |
| `functions_LDM.cuh` | Lagrangian dispersion model |
| `functions_DD.cuh` | Dry deposition & resuspension (Cs-137, ground) |
| `functions_Diffusion.cuh` | Eulerian scalar/concentration transport |
| `functions_BC.cuh` | Boundary conditions (no-slip, open, Neumann) |
| `functions_NNPS.cuh` | Neighbor search, cell indexing, particle init |
| `functions_KNL.cuh` | SPH kernel functions (Wendland, Quintic, etc.) |
| `functions_PREP.cuh` | Strain rate, turbulence prep, correction matrix |
| `functions_OUTPUT.cuh` | VTK binary output, plot variable selection |

## Solver Parameters

Global solver state stored in two arrays accessed via macros:
```cpp
int_t vii[64];  // Integer params: solver_type=vii[0], dim=vii[1], ngpu=vii[29], etc.
Real vif[32];   // Float params: kappa=vif[0], dt=vif[2], time_end=vif[4], etc.
```
Copied to GPU as `__constant__` arrays (`k_vii`, `k_vif`). Macro aliases defined in `function_init.cuh` (e.g., `#define solver_type vii[0]`, `#define dt vif[2]`).

## CUDA Patterns

- All kernels use `thread.x=128`, `block.x=(N-1)/128+1`
- Particle sorting via CUB RadixSort each timestep (cell-based NNPS)
- Two device arrays: `dev_SP1` (sorted input), `dev_P1` (output). At step start: copy P1→SP1 (or SP1→P1 at step 0)
- `cudaDeviceSynchronize()` after every kernel launch
- `__constant__` memory for solver params and property tables
