<!--
1편 확산 검증 절 — 방법 부분만 (2026-10-28). 원고(BnE_manuscript_v2.docx)에는 아직 안 넣었다.
위치는 3장 끝(3.3)을 가정했다. 결과가 나오면 3.3.3 을 채우고 원고에 옮긴다.
[[ ]] = 확인·작성 필요. 런: runs/2026-10-28_cedval-c2-puff
-->

### 3.3 Validation of dispersion against a wind-tunnel puff experiment

Sections 3.1–3.2 test the carrier flow and the deposition step. The dispersion step is tested separately here, against the CEDVAL-LES Complexity 2 experiment of the Environmental Wind Tunnel Laboratory (EWTL), University of Hamburg [[ref: CEDVAL-LES data set / Hertwig et al.]]. The experiment is chosen for three reasons. Its geometry is a building array rather than an isolated obstacle, so plume transport is controlled by street channelling and canopy exchange, as in Section 4. The source is at ground level and most receptors are at pedestrian height. And the data are repeated puff releases rather than a single continuous release, so they give both the mean dispersion and its realisation-to-realisation variability. Only a time-resolved model can reproduce the second.

#### 3.3.1 Experiment

The experiment was carried out in the WOTAN boundary-layer wind tunnel at a scale of 1:225. Values are given at full scale below. The approach flow is a very rough boundary layer (VDI class BL3) with a power-law exponent of 0.27 and a roughness length of 1.53 m. The array consists of 218 buildings of 18 m × 18 m footprint, arranged in 15 staggered rows with 18 m wide streets in both directions, giving a 522 m × 522 m model area and a plan-area density of λp = 0.26. The buildings have four roof types: 109 flat roofs at 18 m, 37 flat roofs at 27 m, and 72 pitched roofs with 18 m eaves and a 27 m ridge, half with the ridge parallel and half perpendicular to the wind. The mixture of heights makes the canopy top irregular, unlike the uniform arrays of Section 3.1.

Ethane was released at ground level on the array centreline, 17.6 m upwind of the first row, in puffs of 0.49 s (model scale). Concentration time series were recorded at 44 receptors with 151–339 puffs per receptor. Twenty-five receptors lie at z = 2 m along four cross-streets, at 45, 117, 225 and 477 m downwind of the source, spanning the plume laterally. The other nineteen form vertical profiles up to z = 90 m in the same cross-streets.

#### 3.3.2 Simulation set-up

The array was reproduced at full scale with the same numerical settings as Section 4 (Table 3), so the comparison tests the configuration actually used for the city. The particle spacing is Δx = 2.5 m, which resolves a building width with 7.2 particles. The domain is 815 m × 735 m × 216 m, with 9.2 × 10⁶ particles. The buildings, including the pitched roofs and gable ends, are represented by 4.0 × 10⁵ surface markers at about 1 m spacing. The inflow is a power-law profile with the experimental exponent of 0.27; the reference speed at the receptor reference height of 13.5 m is 4.7 m/s. The lateral boundaries are periodic with a period of 717.5 m, leaving a 196 m gap between the array and its periodic image. A sponge layer damps the flow above z = 150 m, which is above the highest receptor. The ground uses the MOST wall function of Section 2.3 with z₀ = 0.1 m. Ethane does not deposit, so the deposition step of Section 2.4.4 is switched off and particles reflect at the ground and at building surfaces.

The flow was developed for 200 s. Sixty puffs were then released at 5 s intervals, each consisting of 2 × 10⁴ Lagrangian particles emitted over 0.5 s from z = 1 m. The interval is about five times the canopy time scale H/U_ref ≈ 4–6 s, so successive puffs meet different realisations of the canopy turbulence. Each particle carries the index of its puff, so overlapping puffs are counted separately. The simulation ran to t = 1000 s, long enough for the last puff to leave the farthest receptor line. Particle positions were stored every 1 s. The concentration at a receptor was obtained by counting the particles of each puff in a sampling volume of 5 m × 5 m × 4 m centred on the receptor (from the ground to 4 m for the pedestrian-height receptors).

The comparison uses three quantities.

1. **Ensemble-mean dosage.** For each puff the dosage is the time integral of concentration at a receptor. In stationary flow, the mean dosage of an instantaneous release equals, up to the source strength, the mean concentration of a continuous release at the same point. The dosage pattern therefore tests the mean dispersion, the same quantity that a steady-state or mean-field model predicts.
2. **Puff-to-puff variability.** The distribution over puffs of the dosage, the peak concentration and the peak arrival time at each receptor. These are what a time-averaged flow cannot supply, and they connect directly to the transient-versus-mean comparison of Sections 3.2 and 4.4.
3. **Statistical performance.** The fraction of predictions within a factor of two of the observations (FAC2), the fractional bias (FB) and the normalised mean square error (NMSE), evaluated on the mean dosage at all receptors. The acceptance criteria for urban dispersion of Hanna and Chang [[ref: Hanna & Chang 2012, Bound.-Layer Meteorol.]] are used (FAC2 ≥ 0.30, |FB| ≤ 0.67, NMSE ≤ 6).

Concentrations and times are compared in the dimensionless form of the data set, C* and T* [[to confirm: the documentation that comes with the data does not state the reference length and velocity that define T* and C*; check with EWTL before converting]]. Because the conversion is not yet fixed, the shapes of the lateral and vertical dosage profiles and the ratios of arrival times between receptor lines are also compared. These do not depend on the conversion.

#### 3.3.3 Results

[[Results to be written after the run.]]

<!--
그림·표 계획 (번호는 원고에 넣을 때 맞춤)

Fig. A  실험 배치 (3.3.1 에서 인용) — 결과 없이 지금 그릴 수 있음
  (a) 평면도: 건물 218동을 지붕 종류별 색/무늬로, 방출점 ★, 측정점 + (지면 2 m 는 채운 점, 연직 분포 위치는 빈 원),
      바람 화살표, 측정선 거리 45/117/225/477 m 표기
  (b) 지붕 네 종류 측면 모식 (18 평, 27 평, 경사 2종)
  자료: input/CEDVAL_C2/buildings.csv, validation/cedval/setup/roof_pattern.txt

Fig. B  평균 dosage 횡분포 (핵심 그림)
  네 측정선(45/117/225/477 m) 각 한 패널, 가로축 y, 세로축 평균 dosage(각 선의 합 또는 최댓값으로 정규화)
  실험: 점 + 퍼프 간 25–75 % 막대 / 모델: 선 + 25–75 % 음영
  → 플룸 폭·중심 치우침·길 채널링을 한눈에

Fig. C  평균 dosage 연직분포
  측정점 연직 분포가 있는 지점(117, 225, 477 m 중심선 등) 패널, 가로축 dosage, 세로축 z/H(H=18 m), 캐노피 높이 18·27 m 점선

Fig. D  퍼프 간 변동 (전이 계산의 고유 결과 — 3.2·4.4 와 연결)
  대표 측정점 3 곳(가까움·중간·멂, 모두 중심선 z=2 m)
  (a) 도착 시간(최대농도 시각) 분포  (b) 최대 농도 분포  (c) dosage 분포 — 실험 히스토그램 vs 모델
  보조: 모든 측정점의 dosage 변동계수(CV) 실험 vs 모델 산점도

Fig. E  산점도 — 모든 측정점 평균 dosage, 모델 vs 실험, 1:1 선과 factor-2 띠, 지면/연직 기호 구분

Table X  FAC2, FB, NMSE (지면 측정점, 전체) + 합격 기준

우선순위: 지면 공간이 좁으면 Fig. B + E + D 세 장으로 충분, A 는 setup 그림과 합칠 수 있음, C 는 부록으로.
-->
Fig. B: lateral profiles of mean dosage on the four cross-streets, model vs. experiment, with puff-to-puff interquartile ranges. Fig. C: vertical profiles of mean dosage. Fig. D: distributions of peak arrival time, peak concentration and dosage at three centreline receptors. Fig. E: scatter of modelled vs. observed mean dosage at all receptors with factor-of-two bounds. Table X: FAC2, FB, NMSE against the urban acceptance criteria.

[[Expected limitation to discuss: the inflow is a power-law profile without resolved turbulent fluctuations (Section 2.3), whereas the wind-tunnel approach flow is fully turbulent. Because the source is only 17.6 m upwind of the array, the plume's initial spread is set mostly by the array rather than by the approach flow, but the first receptor line (45 m) is where any difference will be largest.]]
