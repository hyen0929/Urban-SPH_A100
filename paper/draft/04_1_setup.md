# 4.1 OSM-based geometry and simulation setup — 교체판

> Draft 2026-10-20. 원고 4.1 을 정사각 도메인(13 소스) 자료로 교체한 판.
> **도메인 여백·AIJ 권고 관련 단락은 의도적으로 비워 두었다** — 사용자가 직접 채운다.
> 수치 출처: `runs/2026-09-15_leipzig-dt01-sponge-uy/env/`, `case_config.cuh` CASE_ID 1002,
> `input/Leipzig_s10_dt01_full/`, `RESULTS.md` §1.

---

The claim that a meshless carrier flow removes the real-city geometry bottleneck is tested
by applying the framework, without modification, to the centre of Leipzig, Germany. Because
the SPH-IBM solver of Section 2.3 requires only an unstructured point distribution on the
solid surfaces rather than a body-fitted volumetric mesh, the preprocessing reduces to a
sequence of automated geometric operations: building footprints are retrieved from OSM [43]
as closed polygons carrying a height and, where available, descriptors of non-prismatic roof
forms; these are exported as STL surfaces and sampled into uniformly distributed points at
the fluid particle spacing, each carrying a position and outward normal, and supplied as
immersed-boundary input (Fig. 7). Matching the two spacings ensures that the kernel-based
distribution of the IBM forcing acts consistently on both sides of the solid–fluid
interface. Terrain across the densely built centre is treated as a single ground plane at
the domain-mean elevation.

The resulting surface representation comprises 553,904 immersed-boundary markers. Over the
built-up area the fabric consists of 34,255 building cells with a mean height of 18.9 m and
a median of 17.9 m, with a single isolated tower reaching 115 m; it is a typical
central-European perimeter-block morphology in which the blocks enclose interior courtyards.
Fluid particles that fall inside buildings are deactivated before the simulation starts.
They are identified from a column-maximum height field constructed from the surface markers,
a particle being interior if its height does not exceed the maximum marker height in its
horizontal column. This requires neither a watertight mesh nor surface normals, which
matters because the OSM-derived geometry is not watertight.

The computational domain measures 1,315 m in the streamwise direction, 1,315 m laterally and
380 m vertically, with the ground plane at z = 0 and the built-up area occupying
x = 200–1,118 m and y = −452–453 m. The fluid phase is discretised with a uniform particle
spacing of 2.5 m, yielding 43,046,441 SPH particles, and is integrated with a fixed time step
of 0.01 s for 300 s of physical time on a single NVIDIA A100 GPU, which takes 24.5 h of
wall-clock time.

> **[여백·AIJ 권고 단락 — 사용자가 채울 자리]**
> 참고 수치: 상류 200 m, 하류 198 m, 측면 205 m, 연직 380 m. H_max = 115 m 기준으로는
> 각각 1.74H / 1.72H / 1.78H / 3.30H 이고, 캐노피 평균고 18.9 m 기준이면 하류가 10.5 배다.
> 측면은 주기 경계이므로 AIJ 의 이격 기준이 적용되지 않는다(주기 1,297.5 m = 519 × 2.5 m).

The lateral boundaries are periodic, with the period set to an exact multiple of the particle
spacing so that the buffer planes map onto interior planes without interpolation. The inlet
and outlet are open, each sustained by a four-plane buffer band, and a Rayleigh damping layer
is applied near the domain top. The ground is a no-slip surface formed by boundary particles.
The Monin–Obukhov surface-layer parameterisation of Section 2.3, which is exercised in the
stability verification of Section 3.1, is not activated in the present runs; the stratification
is neutral throughout and the ground stress is carried by the resolved no-slip condition.

The inflow is the power-law profile

  u(z) = U_ref (z / z_ref)^α,

with U_ref = 10 m s⁻¹, z_ref = 371.25 m and α = 0.19, an exponent representative of a
moderately rough suburban surface [33] and consistent with the AIJ range [6]. At the release
height of 5 m this gives 4.41 m s⁻¹.

Thirteen continuous point sources are placed within the built-up area at a release height of
5 m, at locations spanning the principal morphological situations of a European city centre:
open intersections, dense street canyons, the immediate vicinity of isolated tall structures,
enclosed courtyards and open plazas adjacent to compact blocks. Each source emits 100,000
tracers continuously between t = 150 s and t = 300 s, giving 1,299,844 tracers for the
ensemble. Emission rate, emission schedule and random seed are identical across the thirteen
sources, so that the only difference between them is position and any difference in the
result is attributable to the urban fabric and to the flow it produces. The carrier flow is
taken to be statistically stationary from t = 150 s, which is when the dispersion model is
activated; time is measured from activation in what follows. The complete parameter set is
given in Table 3.

**Table 3.** Simulation conditions for the Leipzig city-centre case

| Item | Value |
|---|---|
| Domain | 1,315 m × 1,315 m × 380 m |
| Built-up extent | x = 200–1,118 m, y = −452–453 m |
| Building cells / IBM markers | 34,255 / 553,904 |
| Building height (mean / median / max) | 18.9 m / 17.9 m / 115 m |
| Particle spacing | 2.5 m |
| SPH particles | 43,046,441 |
| Time step | 0.01 s |
| Physical simulation time | 300 s |
| LDM activation | t = 150 s |
| Inflow | power law; U_ref = 10 m s⁻¹ at z_ref = 371.25 m; α = 0.19 |
| Lateral boundaries | periodic (period 1,297.5 m) |
| Inlet / outlet | open, four-plane buffer each |
| Top boundary | Rayleigh damping layer |
| Ground | no-slip boundary particles; MOST not activated |
| Stratification | neutral |
| Sources × tracers each | 13 × 100,000 (1,299,844 total) |
| Release height | 5 m |
| Release window | t = 150–300 s, continuous |
| Hardware / wall clock | single NVIDIA A100 GPU / 24.5 h |

---

## 원고 대비 바뀐 수치 (교체 확인용)

| 항목 | 원고 (7소스) | 교체판 (13소스) |
|---|---|---|
| 도메인 | 2,250 × 1,750 × 400 m | 1,315 × 1,315 × 380 m |
| SPH 입자 | 63,598,991 | 43,046,441 |
| 물리 시간 | 500 s | 300 s |
| 소스 × 트레이서 | 7 × 100,000 | 13 × 100,000 |
| z_ref | 200 m | 371.25 m |
| 측면 경계 | (명시 없음) | 주기 |
| 지표 | MOST wall model | no-slip, MOST 미적용 |
