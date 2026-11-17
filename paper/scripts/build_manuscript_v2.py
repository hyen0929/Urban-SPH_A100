#!/usr/bin/env python3
"""BnE 원고 v2 — 기존 원고(archive/BnE_manuscript(1-4chapter).docx)를 이어 쓴다.

새로 씀 : Abstract, 4.3, 4.4, 5(5.1–5.3), 6
교체    : 4.1 영역·설정 문단·Table 3·Fig 7(c), 4.2 전체(Fig 8–10), 4.5 결과 문단, 4.5 수식 내 수치
고침    : 서론 마지막 문단과 기여 3·4, Table 1 채울 수 있는 칸
표시    : 노란 형광 = 저자 확인·작성이 필요한 자리

문장 안 표기:  ^{..} 위첨자,  _{..} 아래첨자,  [[HL:..]] 노란 형광
"""
import copy, re, sys
import docx
from docx.shared import Cm
from docx.enum.text import WD_COLOR_INDEX, WD_ALIGN_PARAGRAPH
from docx.text.paragraph import Paragraph

SRC, OUT = sys.argv[1], sys.argv[2]
FIG = '/home/kdh/Urban-SPH/paper/figs/paper1'
d = docx.Document(SRC)
P = list(d.paragraphs)
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'

# ── 도우미 ───────────────────────────────────────────────────────────
TOK = re.compile(r'(\^\{[^}]*\}|_\{[^}]*\}|\[\[HL:.*?\]\])', re.S)


def runs_from(p, text, bold=False):
    for t in TOK.split(text):
        if not t:
            continue
        if t.startswith('^{'):
            r = p.add_run(t[2:-1]); r.font.superscript = True
        elif t.startswith('_{'):
            r = p.add_run(t[2:-1]); r.font.subscript = True
        elif t.startswith('[[HL:'):
            r = p.add_run(t[5:-2]); r.font.highlight_color = WD_COLOR_INDEX.YELLOW
        else:
            r = p.add_run(t)
        if bold:
            r.font.bold = True


def clear(p):
    for el in list(p._p):
        if el.tag != W + 'pPr':
            p._p.remove(el)


def set_text(p, text):
    clear(p); runs_from(p, text)


def after(anchor, template, text=None):
    """template 문단의 서식(pPr)을 복사한 새 문단을 anchor 뒤에 넣는다."""
    el = copy.deepcopy(template._p)
    for c in list(el):
        if c.tag != W + 'pPr':
            el.remove(c)
    anchor._p.addnext(el)
    p = Paragraph(el, anchor._parent)
    if text is not None:
        runs_from(p, text)
    return p


def picture(anchor, path, width_cm):
    p = after(anchor, IMGTPL)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(path, width=Cm(width_cm))
    return p


def caption(anchor, num, text):
    p = after(anchor, CAPTPL)
    runs_from(p, f'[Fig. {num}]', bold=True)
    runs_from(p, ' ' + text)
    return p


def drop(p):
    p._p.getparent().remove(p._p)


def math_replace(p, old, new):
    n = 0
    for el in p._p.iter(M + 't'):
        if el.text == old:
            el.text = new; n += 1
    assert n >= 1, (old, p.text[:60])


def tail_replace(p, marker, new):
    """일반 텍스트 run 중 marker 가 시작되는 지점부터 문단 끝까지를 new 로 바꾼다(수식은 보존)."""
    ts = [el for el in p._p.iter(W + 't')]
    for k, el in enumerate(ts):
        if el.text and marker in el.text:
            el.text = el.text[:el.text.index(marker)]
            for e2 in ts[k+1:]:
                e2.text = ''
            break
    else:
        raise ValueError(marker)
    runs_from(p, new)


def find(prefix):
    for p in d.paragraphs:
        if p.text.strip().startswith(prefix):
            return p
    raise ValueError(prefix)


BODY, H2, CAPTPL, IMGTPL = P[178], P[177], P[199], P[180]

# ══════════════════════════════════════════════════════════════════════
# Abstract
# ══════════════════════════════════════════════════════════════════════
set_text(P[5], (
    "Building-resolved prediction of urban dispersion and dry deposition is limited by the cost of meshing "
    "complex geometry and by the practice of advecting pollutants on time-averaged flows. A mesh-free framework is presented in which an explicit incompressible "
    "smoothed particle hydrodynamics large-eddy simulation with immersed-boundary buildings drives a "
    "Lagrangian stochastic dispersion model through the solver's own kernel, and dry deposition is resolved "
    "as individual tracer–surface events on the ground, walls and roofs. The framework is verified against "
    "Architectural Institute of Japan wind-tunnel benchmarks and an analytical solution and applied to a 1.3 km × 1.3 km district of Leipzig, Germany, built directly from OpenStreetMap and "
    "resolved with 43 million particles on a single GPU. Thirteen releases differing only in position yield "
    "exposed areas, peak pedestrian-level concentrations and deposited fractions differing by factors of "
    "14, 20 and 126, ordered by surrounding building density. Dense surroundings do not increase "
    "pedestrian-level exposure but concentrate it into smaller areas, and suppress ground deposition "
    "through the low friction velocity of sheltered ground while shifting deposition towards building "
    "envelopes. [[HL:In an idealised cube array, advecting the releases on the time-averaged flow "
    "conserves the total deposited mass but redistributes it between sources by more than a factor of 40;]] "
    "in the real district it also reduces the deposited mass by 21% and contracts deposition footprints by "
    "up to 62%, because averaging removes wake-driven transport and the effect of near-surface fluctuations "
    "on the nonlinear deposition velocity. Near-source exposure and surface contamination in dense urban "
    "districts require time-resolved, building-resolved simulation."))

# ══════════════════════════════════════════════════════════════════════
# 1 Introduction — 결과에 맞춰 마지막 문단과 기여 3·4 수정
# ══════════════════════════════════════════════════════════════════════
set_text(find('This construction also lets the framework'), (
    "This construction also lets the framework address the second open question directly. Driving the "
    "instantaneous and the time-averaged flow field separately from the same sources shows that the spatial "
    "distribution of deposition cannot be recovered from the averaged field. In an idealised array the total "
    "deposited mass is conserved while individual sources are redistributed by more than an order of "
    "magnitude; in real urban geometry, where the deposition velocity responds to the local near-surface "
    "flow, the averaged field additionally underestimates the total. The failure of real-city deposition "
    "prediction is thus more than a parameterisation problem: the pattern information is lost the moment the "
    "field is averaged in time."))
set_text(find('Quantification, by coupling'), (
    "Quantification, by coupling the instantaneous and the time-averaged flow field to the same sources, of "
    "the deposition error introduced by temporal averaging: a redistribution between sources by more than a "
    "factor of 40 at a total conserved to within 2% in an idealised array, and in real geometry a 21% "
    "underestimate of the total, reaching 65% for individual sources, with deposition footprints contracted "
    "by up to 62%."))
set_text(find('Application of the full pipeline'), (
    "Application of the full pipeline to a district of Leipzig, Germany, extracted from OSM, with thirteen "
    "controlled releases, showing that the surrounding building density organises near-source exposure and "
    "deposition: it concentrates rather than increases pedestrian-level exposure, and it partitions "
    "deposition between the ground and building envelopes through the local surface friction velocity."))

# ══════════════════════════════════════════════════════════════════════
# Table 1 — 본문에 있는 값으로 채우고, 참조값은 형광 빈칸
# ══════════════════════════════════════════════════════════════════════
T1 = d.tables[0]
rows = [
    ("Carrier flow", "AIJ Case A (isolated high-rise); Case C1–C3 (3 × 3 cube array) [6,38]",
     "time-mean streamwise velocity", "correlation coefficient (CC)",
     "A: 0.840 (vertical) / 0.899 (horizontal); C1: 0.657; [[HL:C2: ]]; C3: 0.815", "[[HL:reference CFD]]"),
    ("Wall / surface", "analytic rough-wall logarithmic law; imposed surface heat flux",
     "mean velocity profile; elevated-to-near-ground tracer ratio", "profile agreement; ratio",
     "log law recovered for z_{0} = 0.1, 1.0, 2.0 m; ratio 0.07 / 0.64 / 1.55 (stable / neutral / unstable)",
     "log law; monotonic stability dependence"),
    ("Dispersion", "Green's function, continuous point source in uniform flow; Eulerian SPH scalar transport",
     "centroid, variance growth, concentration PDF", "deviation from analytical solution",
     "all three reproduced; no spurious lateral drift", "analytical; Eulerian SPH"),
]
for r, vals in zip(T1.rows[1:], rows):
    for c, v in zip(r.cells, vals):
        p = c.paragraphs[0]; clear(p); runs_from(p, v)

# ══════════════════════════════════════════════════════════════════════
# 3.2 — 재계산 대기 표시
# ══════════════════════════════════════════════════════════════════════
after(find('3.2 Transient versus time-averaged carrier flow'), BODY, (
    "[[HL:(작성 메모 — 삭제할 것) 3.2절의 수치(Table 2, Figs. 3–6)는 2025-11 CaseCS 런 기준이며 재계산 후 갱신 "
    "예정. 재계산은 MOST_boundary = 0(균일 u*)을 유지해야 4.4절의 대비 논거가 성립함. 본문의 'time step 0.01 s, "
    "averaged over 50 s'는 복원된 설정(dt 0.005 s, 평균 30 s)과 다르므로 확인 필요.]]"))

# ══════════════════════════════════════════════════════════════════════
# 4.1
# ══════════════════════════════════════════════════════════════════════
img_c = P[181]
picture(img_c, f'{FIG}/Fig07_domain_sources.png', 10.5); drop(img_c)
set_text(P[182], "")
runs_from(P[182], "[Fig. 7]", bold=True)
runs_from(P[182], " OSM-to-SPH geometry pipeline for the Leipzig city centre: (a) building footprints retrieved "
          "from OSM; (b) surface particle distribution used as immersed-boundary input; (c) built-up area with "
          "the thirteen release points, grouped by the dispersion regime identified in Section 4.2.")
drop(P[183]); drop(P[184])

set_text(P[186], (
    "The computational domain measures 1,315 m in the streamwise and lateral directions and 380 m vertically, "
    "with the built-up area occupying x = 200–1,118 m and y = −452–453 m (Fig. 7c). Over this area the "
    "surface is represented by 553,904 immersed-boundary markers. On a 2.5 m raster the buildings occupy "
    "34,255 cells with a mean height H̄ = 18.9 m and a median of 17.9 m, and a single tower reaches 115 m; the "
    "fabric is a typical central-European perimeter-block morphology in which most blocks enclose courtyards. "
    "Fluid particles falling inside buildings are deactivated before the run on the basis of a column-maximum "
    "height field constructed from the surface markers, which requires neither a watertight surface nor "
    "surface normals. [[HL:(저자 작성: 도메인 여백·AIJ 권고 관련 서술. 참고 — 상류 200 m, 하류 198 m, 측면 205 m"
    "(주기 경계), 연직 380 m; H_max = 115 m 기준 1.74H / 1.72H / 1.78H / 3.30H, 평균고 기준 하류 10.5H̄.)]]"))
set_text(P[187], (
    "The fluid phase is discretised with a uniform particle spacing of 2.5 m, yielding 43,046,441 SPH "
    "particles, and integrated with a fixed time step of 0.01 s for 300 s of physical time on a single NVIDIA "
    "A100 GPU, which required 51 h of wall-clock time. The lateral boundaries are periodic, the inlet and "
    "outlet are open with four-plane buffer layers, a Rayleigh damping layer is applied below the top "
    "boundary, and the ground carries the MOST wall model of Section 2.3 with z_{0} = 0.1 m. The inflow is "
    "the power-law profile"))
math_replace(P[191], '200', '371.25')
tail_replace(P[191], 'The stratification is neutral.', (
    "The stratification is neutral. Thirteen continuous point sources are placed at a release height of 5 m, "
    "at locations spanning the principal morphological situations of a European city centre: open squares "
    "and intersections, street canyons of different orientation, the vicinity of tall structures, and "
    "enclosed courtyards (Fig. 7c). Each source emits 99,988 tracers at a constant rate between t = 150 s "
    "and t = 300 s. Emission rate, schedule and random seed are identical across sources, so that any "
    "difference in the outcome is attributable to the urban fabric and to the flow it produces. The complete "
    "parameter set is given in Table 3."))

T3 = d.tables[2]
t3 = [("Item", "Value"),
      ("Domain", "1,315 m × 1,315 m × 380 m"),
      ("Built-up area", "x = 200–1,118 m, y = −452–453 m"),
      ("Building height (mean / median / max)", "18.9 m / 17.9 m / 115 m"),
      ("Immersed-boundary markers", "553,904"),
      ("Particle spacing", "2.5 m"),
      ("SPH particles", "43,046,441"),
      ("Time step", "0.01 s"),
      ("Physical simulation time", "300 s"),
      ("Inflow", "power law; U_{ref} = 10 m s^{−1} at z_{ref} = 371.25 m; α = 0.19"),
      ("Boundaries", "lateral periodic; open inlet and outlet; Rayleigh damping layer at top; "
                     "ground no-slip with MOST wall model (z_{0} = 0.1 m)"),
      ("Stratification", "neutral"),
      ("Sources × tracers each", "13 × 99,988; release height 5 m; t = 150–300 s"),
      ("Hardware / wall clock", "single NVIDIA A100 GPU / 51 h")]
for r, vals in zip(T3.rows, t3):
    for c, v in zip(r.cells, vals):
        p = c.paragraphs[0]; clear(p); runs_from(p, v)

# ══════════════════════════════════════════════════════════════════════
# 4.2 — 전면 교체
# ══════════════════════════════════════════════════════════════════════
h42 = P[195]
set_text(h42, "4.2 Carrier flow and source-dependent dispersion")
for i in range(196, 212):
    drop(P[i])

a = after(h42, BODY, (
    "The resolved carrier flow is summarised in Fig. 8. At mid-canopy height the instantaneous vorticity "
    "(Fig. 8a) concentrates in the shear layers shed from windward building edges and along canyon walls, "
    "while open squares and broad streets carry comparatively weak vorticity; the approach flow is nearly "
    "irrotational upstream of the first buildings and acquires its turbulence within the first block rows. "
    "The time-averaged wind at the lowest computational level, z = 1.25 m (Fig. 8b), is fast only along "
    "streets aligned with the approach flow and in the gaps between detached buildings, and it is nearly "
    "stagnant in courtyards and in narrow canyons normal to the wind; its median over the open ground of the "
    "built-up area is 0.11 U_{ref}. The horizontally averaged mean wind (Fig. 8c) is reduced below roof "
    "level, to 40% of the undisturbed inflow profile at the lowest level, and rejoins the inflow profile "
    "above about three canopy heights."))
a = after(a, BODY, (
    "Above the roughness sublayer the averaged profile follows a logarithmic law with a canopy-scale friction "
    "velocity of 0.5–0.7 m s^{−1}, the value depending on the fitting range because the upper profile remains "
    "tied to the imposed inflow. The friction velocity returned by the wall model at the ground between the "
    "buildings is markedly smaller, 0.20 m s^{−1} as a domain average, so the ground carries roughly a tenth "
    "of the total momentum flux and the remainder is absorbed as form drag by the buildings, as expected for "
    "a dense urban canopy [45,46]. It is this ground value, not the canopy-scale one, that enters the "
    "deposition velocity of Section 2.4.4. The resolved kinetic energy is constant to within 1% from "
    "t ≈ 25 s, and the subgrid kinetic energy grows by about 5% per 50 s over the analysis window; the "
    "carrier flow is therefore taken as statistically stationary from t = 150 s, when the releases begin."))
a = picture(a, f'{FIG}/Fig08_carrier_flow.png', 15.9)
a = caption(a, 8, (
    "Carrier flow over the Leipzig district. (a) Instantaneous vorticity magnitude at mid-canopy height "
    "(z = 11.25 m) at t = 300 s; (b) time-averaged horizontal wind speed at the lowest computational level "
    "(z = 1.25 m), normalised by U_{ref}; (c) horizontally averaged time-mean wind speed over the fluid region "
    "compared with the imposed inflow profile, H̄ denoting the mean building height."))
a = after(a, BODY, (
    "The thirteen sources release simultaneously, and Fig. 9 shows the resulting time-integrated "
    "concentration in the pedestrian layer (z ≤ 2 m) over the release period, per released tracer. The "
    "joint field is a superposition of thirteen individually distinct plumes rather than statistically "
    "similar copies of a single canonical plume. Their directions are set by the street network as much as "
    "by the approach flow: the exposure centroids of the four plumes released in open space depart from the "
    "mean wind direction by 13–88°, and for six of the nine sources released within blocks the exposure "
    "centroid lies upwind of the release point, the signature of recirculating canyon and courtyard flow."))
a = picture(a, f'{FIG}/Fig09_pedestrian_exposure_all.png', 13.5)
a = caption(a, 9, (
    "Time-integrated concentration in the pedestrian layer (z ≤ 2 m) over the release period "
    "t = 150–300 s, per released tracer and summed over the thirteen sources. Open circles mark the release "
    "points."))
a = after(a, BODY, (
    "Examined individually, the plumes fall into three regimes set by the configuration immediately around "
    "the source (Figs. 7c and 10). Plumes released in open space (sources 8, 9, 11 and 13) are advected: "
    "they leave the release area within tens of seconds, spread laterally over the next blocks and reach the "
    "downstream edge of the built-up area (Fig. 10a). Plumes released in streets (sources 1, 3, 4, 5, 7 and "
    "12) are channelled: they follow the street axis, turn at junctions and leak into side streets, so that "
    "their footprint traces the street network rather than the wind direction; the plume from source 4 runs "
    "along an east–west canyon over more than 200 m while remaining one street wide (Fig. 10b). Plumes "
    "released in enclosed courtyards or yards (sources 2, 6 and 10) are trapped and remain within the "
    "enclosure throughout the release (Fig. 10c)."))
a = picture(a, f'{FIG}/Fig10_regimes.png', 15.9)
a = caption(a, 10, (
    "Representative dispersion regimes, shown as the time-integrated pedestrian-layer concentration of a "
    "single source on a common scale and in windows of equal size: (a) advected (source 8, open square); "
    "(b) channelled (source 4, east–west street canyon); (c) trapped (source 10, enclosed courtyard)."))
a = after(a, BODY, (
    "The regimes are not discrete classes but ranges of a continuum ordered by the building plan-area "
    "density λ_{p} within 75 m of the release point [45] (Fig. 11; Spearman rank correlation ρ_{s} with a "
    "permutation p-value, n = 13). With increasing λ_{p} the exposed area, defined as the area over which the "
    "time-integrated concentration exceeds 1% of its maximum, decreases from about 13,000 m^{2} to 900 m^{2} "
    "(ρ_{s} = −0.77, p = 0.004), the exposure-weighted mean distance from the source decreases from 76 m to "
    "7 m (ρ_{s} = −0.81, p = 0.002), and the peak concentration increases twenty-fold (ρ_{s} = +0.76, "
    "p = 0.004). The pedestrian-layer residence time per released tracer, in contrast, shows no dependence "
    "on λ_{p} (ρ_{s} = +0.01, p = 0.98): a tracer released in a courtyard spends about as long near the ground "
    "as one released in a square. Dense surroundings therefore do not raise the near-ground exposure "
    "produced by a release; they compress the same exposure into a smaller area, with correspondingly higher "
    "local concentrations. This distinction underlies the deposition and exposure results that follow."))
a = picture(a, f'{FIG}/Fig11_metrics_vs_lambda.png', 15.9)
a = caption(a, 11, (
    "Source metrics against the building plan-area density λ_{p} within 75 m of the release point (n = 13): "
    "(a) residence time in the pedestrian layer per released tracer; (b) exposed area, where the "
    "time-integrated pedestrian-layer concentration exceeds 1% of its maximum; (c) peak time-integrated "
    "concentration; (d) deposited fraction at t = 300 s; (e) share of the deposited tracers on walls and "
    "roofs; (f) exposure-weighted mean distance from the source. Symbols denote the dispersion regimes of "
    "Fig. 7c; ρ_{s} is the Spearman rank correlation with λ_{p}."))

# ══════════════════════════════════════════════════════════════════════
# 4.3
# ══════════════════════════════════════════════════════════════════════
h43 = find('4.3 Surface deposition and urban morphology')
nxt = Paragraph(h43._p.getnext(), h43._parent)
if not nxt.text.strip() and not nxt._p.xpath('.//w:drawing'):
    drop(nxt)
a = after(h43, BODY, (
    "By the end of the release 6.4% of the released tracers have deposited, and the deposited count is "
    "still increasing at that time; deposition is therefore compared between sources and between flow "
    "treatments rather than interpreted as a final inventory. Of the deposited tracers, 98.8% are on the "
    "ground and 1.2% on building envelopes, the latter divided between roofs (65%) and walls (35%) with "
    "median heights of 15 m and 13 m, close to the mean building height."))
a = after(a, BODY, (
    "The deposited fraction differs by more than two orders of magnitude between sources, from 0.2% for "
    "the release in an enclosed courtyard (source 10) to 26.6% for the release in an open square "
    "(source 13), and it falls steeply with λ_{p} (Fig. 11d; ρ_{s} = −0.86, p < 0.001). The share deposited "
    "on building envelopes follows the opposite trend (Fig. 11e; ρ_{s} = +0.83, p < 0.001): it is below "
    "0.5% for the four advected releases, reaches 13–15% for the releases within blocks (sources 3 and 6), "
    "and 45% for the courtyard release, although the last value rests on only 212 deposited tracers. Fig. 12 "
    "contrasts the two situations. The plume from the open square deposits on the ground in a broad patch "
    "that accumulates where it impinges on the windward face of the next block, along which the few envelope "
    "deposits are found; the plume trapped in the courtyard of source 6 deposits sparsely on the courtyard "
    "floor and preferentially along the surrounding walls and roofs."))
a = picture(a, f'{FIG}/Fig12_deposition_pathway.png', 13.5)
a = caption(a, 12, (
    "Deposited tracers at t = 300 s for (a) source 8, released in an open square, and (b) source 6, released "
    "in an enclosed courtyard. Blue: deposited on the ground; orange: deposited on walls and roofs. The "
    "envelope share of deposition is 0.4% in (a) and 12.9% in (b). The window in (a) is shifted downwind to "
    "include the deposition footprint."))
a = after(a, BODY, (
    "Because the residence time in the pedestrian layer does not depend on λ_{p} (Fig. 11a), the difference "
    "in deposition cannot be attributed to how long the plumes stay near the ground. It follows instead from "
    "where they stay. The deposition velocity of Section 2.4.4 is evaluated with the friction velocity "
    "returned by the wall model from the flow immediately above each ground element, and for particles of "
    "1 μm it depends steeply on that value through the Stokes-number dependence of the impaction efficiency "
    "[11–13]: V_{d} decreases 25-fold, from 3.1 × 10^{−2} to 1.2 × 10^{−3} m s^{−1}, when u_{*} decreases "
    "from 0.5 to 0.2 m s^{−1}. Weighted by the near-ground exposure of each plume, the local ground friction "
    "velocity is 0.12–0.32 m s^{−1} for the advected releases and 0.002–0.16 m s^{−1} for the trapped ones "
    "(ρ_{s} = −0.86 with λ_{p}), and the deposited fraction predicted from the exposure-weighted local "
    "deposition velocity reproduces the ranking of the simulated deposited fractions (ρ_{s} = 0.93). "
    "Sheltered ground in dense blocks thus receives little deposition although the plumes linger above it, "
    "because the near-surface wind that would carry particles through the quasi-laminar sublayer is weak."))
a = after(a, BODY, (
    "Envelope deposition is not subject to this sheltering in the same way. At most elevated wall and roof "
    "points no surface-layer sample of the ground model lies within the averaging stencil, and the envelope "
    "deposition velocity is evaluated with the domain-averaged friction velocity; envelope deposition is "
    "therefore governed by the exposure of the envelope to the plume rather than by local sheltering. As "
    "λ_{p} increases, ground deposition is suppressed while the plume is held among building faces, and the "
    "deposited share shifts towards the envelope. For decontamination planning after an accidental release "
    "this means that the surfaces requiring attention depend on where the release occurs: almost exclusively "
    "the ground for releases in open space, but walls and roofs to a substantial degree for releases within "
    "blocks. Given the simplicity of the envelope model, these shares are interpreted as relative rather "
    "than absolute values (Section 5.3)."))

# ══════════════════════════════════════════════════════════════════════
# 4.4
# ══════════════════════════════════════════════════════════════════════
h44 = find('4.4 Transient versus time-averaged deposition in real geometry')
nxt = Paragraph(h44._p.getnext(), h44._parent)
if not nxt.text.strip() and not nxt._p.xpath('.//w:drawing'):
    drop(nxt)
a = after(h44, BODY, (
    "The comparison of Section 3.2 is repeated in the real geometry. The carrier field of Section 4.1 is "
    "time-averaged over t = 150–299 s, the thirteen releases are advected on the averaged field with the "
    "same release schedule, random seeds and deposition model, and the outcome is compared with that on the "
    "transient field (Fig. 13). The averaged calculation takes 5.4 min on a single GPU."))
a = after(a, BODY, (
    "Three differences appear. First, the averaged field deposits less: 65,952 tracers against 83,103 on "
    "the transient field, 21% fewer in total, with a ratio between 0.35 and 1.13 across sources and a "
    "decrease of 65% for source 11 (Fig. 13c). Second, the deposition footprints contract: the area "
    "containing 90% of the deposited tracers is smaller on the averaged field for ten of the thirteen "
    "sources, by up to 62% (source 11), and branches of the footprint that are fed by unsteady wake "
    "transport disappear entirely, as does the south-eastern branch of the source 11 deposit around the "
    "downstream block (Fig. 13a,b). Third, the pedestrian-layer exposure contracts in the same way: the area "
    "receiving at least one tracer-second per cell decreases by 30%, while the domain-integrated residence "
    "increases by 11% because fewer tracers are removed by deposition."))
a = picture(a, f'{FIG}/Fig13_transient_vs_mean.png', 15.9)
a = caption(a, 13, (
    "Deposition on the transient and the time-averaged carrier flow. (a, b) Deposited tracers per unit area "
    "at t = 300 s for source 11 on (a) the transient and (b) the time-averaged field, on a common scale. "
    "(c) Ratio of the time-averaged to the transient result for each source: deposited tracer count "
    "(circles) and area containing 90% of the deposited tracers (squares)."))
a = after(a, BODY, (
    "The contraction of the footprints is the real-geometry counterpart of the result of Section 3.2: "
    "averaging removes the lateral and vertical exchange produced by unsteady wakes, and the plume is "
    "confined to its mean path. The underestimate of the total, which did not occur in the idealised array, "
    "has an additional origin that is specific to a deposition velocity responding to the local flow. The "
    "friction velocity entering V_{d} is evaluated from the near-surface wind; on the averaged field this "
    "wind is the time mean, whereas on the transient field it fluctuates about that mean. Because V_{d} is a "
    "strongly convex function of u_{*}, the fluctuations raise its mean value even where they leave the mean "
    "friction velocity almost unchanged. Evaluated over the ground of the built-up area, the near-surface "
    "friction velocity differs by only 7% between the instantaneous field at t = 300 s and the averaged field "
    "(0.205 against 0.191 m s^{−1}), but the mean deposition velocity is 44% higher on the instantaneous "
    "field (6.6 against 4.6 mm s^{−1}). In the idealised configuration of Section 3.2 the friction velocity "
    "entering V_{d} was uniform, so this effect was absent and the total was conserved. A time-averaged "
    "carrier flow therefore biases deposition in real geometry in two ways: it misplaces deposition, as in "
    "the idealised case, and it underestimates its magnitude wherever the deposition velocity depends "
    "nonlinearly on the local near-surface flow."))

# ══════════════════════════════════════════════════════════════════════
# 4.5 — 수식 내 수치, 결과 문단
# ══════════════════════════════════════════════════════════════════════
p217 = find('Application to near-field emergency response requires')
for el in p217._p.iter(W + 't'):
    if el.text and '25 m grid' in el.text:
        el.text = el.text.replace('25 m grid', '5 m grid')
p221 = find('where')  # 'where n is the local tracer number density ...'
for p in d.paragraphs:
    if p.text.strip().startswith('where') and 'local tracer number density' in p.text:
        p221 = p; break
math_replace(p221, '5', '2')
p225 = [p for p in d.paragraphs if p.text.strip().startswith('with') and 'LDM output interval' in p.text][0]
tail_replace(p225, 'The quantity is reported in arbitrary units', (
    "Here Δt_{out} = 1 s and T = 150 s. In Figs. 9–11, E_{ped} is divided by z_{ped} and by the number of "
    "released tracers, giving the layer-averaged time-integrated concentration per released tracer in "
    "s m^{−3}; it is a spatial indicator of cumulative near-surface exposure rather than a dose estimate, "
    "since conversion to a dose for a specific substance would require the corresponding source term and "
    "exposure pathway, which are beyond the present scope."))
p226 = find('The resulting field (Fig. 11)')
set_text(p226, (
    "The resulting field (Fig. 9) is highly anisotropic at the block scale and organised by the street "
    "canyons and open spaces of the actual city, in contrast to the symmetric footprint that a Gaussian "
    "plume would yield for the same releases. Two properties established above bear directly on exposure "
    "assessment. First, the exposure produced by a release is not larger in dense surroundings but is "
    "concentrated: compared with a release in open space, a release within a block exposes an area an order "
    "of magnitude smaller at a peak concentration several times higher (Fig. 11b,c). Conventional zoning on "
    "concentric perimeters around the source is poorly suited to either case, since it overstates the "
    "exposed area for trapped releases and misplaces it for advected ones, whose footprints are deflected by "
    "the street network. Second, the near-source exposure is carried by the time-resolved flow: on the "
    "averaged field the exposed area shrinks by 30% (Section 4.4), so that a mean-flow assessment would "
    "understate the area requiring protective action. The spatial information supplied by the framework is "
    "therefore of direct use to exposure assessment and zoning, at a cost compatible with a single GPU at "
    "the scale of a city district."))
# 옛 Fig. 11 (노출 지도) 은 Fig. 9 로 대체 — 그림·캡션 제거
for p in list(d.paragraphs):
    if p.text.strip().startswith('[Fig. 11] Cumulative pedestrian-level exposure'):
        prev = Paragraph(p._p.getprevious(), p._parent)
        if prev._p.xpath('.//w:drawing'):
            drop(prev)
        drop(p); break

# ══════════════════════════════════════════════════════════════════════
# 5 Discussion
# ══════════════════════════════════════════════════════════════════════
h5 = find('5. Discussion')
nxt = Paragraph(h5._p.getnext(), h5._parent)
if not nxt.text.strip():
    drop(nxt)
a = after(h5, H2, "5.1 What time- and building-resolved simulation adds")
a = after(a, BODY, (
    "The results of Sections 3 and 4 point to three properties of near-source dispersion and deposition in "
    "dense urban districts that are inaccessible to flat-terrain and mean-flow methods. First, the outcome of "
    "a release depends on its position within the urban fabric to a degree that is large in its own right: "
    "thirteen releases of identical strength and height in one district differ by more than two orders of "
    "magnitude in deposited fraction and by more than an order of magnitude in exposed area and peak "
    "concentration. Second, this dependence is organised. The surrounding building density orders the "
    "outcomes, so that a simple morphological indicator can rank release locations, but it acts through "
    "mechanisms, the confinement of the plume and the sheltering of the ground, that set magnitudes no "
    "indicator captures. Third, the spatial pattern of deposition and exposure is carried by the "
    "time-resolved flow. Averaging the carrier field removes the transport by unsteady wakes and, in real "
    "geometry, the effect of near-surface fluctuations on the deposition velocity; a mean-flow calculation "
    "can therefore be correct in total and wrong in pattern, as in the idealised array, or wrong in both, as "
    "in the city."))
a = after(a, BODY, (
    "These properties matter most where the consequence of a release is decided near the source, as in "
    "accidental or deliberate releases of hazardous or radioactive material within a city, where an "
    "assessment must identify the affected streets and building surfaces rather than the amount deposited "
    "over a district. For such assessments the present results suggest that protective zones should be "
    "derived from the street-scale footprint rather than from concentric perimeters, that decontamination "
    "planning should consider building envelopes explicitly for releases within blocks, and that the carrier "
    "flow should be resolved in time rather than averaged before the dispersion is computed."))
a = after(a, H2, "5.2 Computational cost")
a = after(a, BODY, (
    "The real-city calculation of Section 4, with 43 million SPH particles, 300 s of physical time and "
    "1.3 million tracers with deposition on ground and building surfaces, ran on a single NVIDIA A100 GPU in "
    "51 h of wall-clock time. The preprocessing from OSM to immersed-boundary input is automated and involves "
    "no volume mesh, so a new district or a modified geometry requires no manual meshing effort. Because the "
    "tracers are passive, releases do not interact, and any number of release points can be evaluated in a "
    "single calculation. Once a carrier field has been computed, scenarios that do not alter the flow, such "
    "as other release locations, schedules or deposition parameters, can be evaluated on the stored "
    "time-averaged field in minutes, as in Section 4.4, or, where the time-resolved flow is required, by "
    "repeating the transient calculation."))
a = after(a, H2, "5.3 Limitations")
a = after(a, BODY, (
    "Several limitations qualify these results. First, the release and tracking period is finite; "
    "deposition was still increasing at the end of the 150 s release, so deposited fractions are compared "
    "between sources and between flow treatments rather than interpreted as final inventories. Second, a "
    "single district, a single approach-flow direction and speed, and neutral stratification were "
    "considered. The source-to-source variability reported here is a property of this configuration; its "
    "dependence on wind direction and stability is the subject of ongoing work, and stratified conditions "
    "require the advection scheme for momentum to be extended. Third, the ground friction velocity is "
    "obtained from the wall model at the first computational level, 1.25 m above the ground and within the "
    "urban canopy, where the logarithmic profile underlying the wall function is not strictly valid. "
    "Because the deposition velocity of fine particles depends steeply on this value, absolute deposited "
    "fractions are sensitive to it, whereas their ranking between sources, which is set by the contrast "
    "between sheltered and exposed ground, is robust."))
a = after(a, BODY, (
    "Fourth, envelope deposition uses a 2 m capture layer, the domain-averaged friction velocity at most "
    "elevated points and no resuspension from walls, so the envelope shares are relative indicators. "
    "Resuspension from the ground is included in the model but is negligible here: at the highest local "
    "friction velocities encountered, about 0.3 m s^{−1}, the resuspension rate corresponds to less than 1% "
    "of the deposited tracers over the release period. Fifth, the dispersion layer has been verified against "
    "an analytical solution but not yet validated against measured concentrations within building arrays. "
    "[[HL:(저자 작성: 확산 검증 계획(CEDVAL-LES Complexity 2 등) — 데이터 확보 후 3.1절로 이동하거나 여기서 서술.)]] "
    "[[HL:(저자 작성: 도메인 여백 관련 한계.)]]"))

# ══════════════════════════════════════════════════════════════════════
# 6 Conclusions
# ══════════════════════════════════════════════════════════════════════
h6 = find('6. Conclusions')
nxt = Paragraph(h6._p.getnext(), h6._parent)
if not nxt.text.strip():
    drop(nxt)
a = after(h6, BODY, (
    "A mesh-free framework coupling an EISPH large-eddy simulation with immersed-boundary buildings to a "
    "Lagrangian stochastic dispersion and deposition model was presented, verified, and applied to a city "
    "district built directly from OpenStreetMap on a single GPU. The main conclusions are as follows."))
for t in (
    "(1) Building-resolved dispersion and deposition in real urban geometry are obtained without a volume "
    "mesh. The carrier flow reproduces the momentum partition of a dense urban canopy, with the ground "
    "carrying about a tenth of the total momentum flux and the buildings the remainder.",
    "(2) Thirteen identical releases in one district differ by factors of 14, 20 and 126 in exposed area, "
    "peak pedestrian-level concentration and deposited fraction. The surrounding building density orders "
    "these outcomes but does not increase the near-ground exposure produced by a release; it concentrates "
    "the same exposure into a smaller area.",
    "(3) Ground deposition is suppressed in dense surroundings by the low friction velocity of sheltered "
    "ground rather than by a shorter near-ground residence of the plume, and the deposited share shifts "
    "from the ground towards building envelopes, from below 0.5% for releases in open space to 13–45% for "
    "releases within blocks.",
    "(4) [[HL:Advecting the releases on the time-averaged flow conserves the total deposited mass in an "
    "idealised array but redistributes it between sources by more than a factor of 40.]] In real geometry it "
    "additionally underestimates the total by 21% and contracts the deposition footprints by up to 62%, "
    "because the deposition velocity depends nonlinearly on the fluctuating near-surface flow."):
    a = after(a, BODY, t)
a = after(a, BODY, (
    "Near-source exposure and surface contamination in dense urban districts therefore require simulation "
    "that resolves both the buildings and the time variation of the flow."))

# ══════════════════════════════════════════════════════════════════════
# 참고문헌 추가
# ══════════════════════════════════════════════════════════════════════
last = [p for p in d.paragraphs if p.text.strip().startswith('[44]')][0]
a = after(last, last, (
    "[45] C.S.B. Grimmond, T.R. Oke, Aerodynamic properties of urban areas derived from analysis of surface "
    "form, J. Appl. Meteorol. 38 (1999) 1262-1292. https://doi.org/10.1175/1520-0450(1999)038<1262:APOUAD>2.0.CO;2."))
a = after(a, last, (
    "[46] R.W. Macdonald, R.F. Griffiths, D.J. Hall, An improved method for the estimation of surface "
    "roughness of obstacle arrays, Atmos. Environ. 32 (1998) 1857-1864. "
    "https://doi.org/10.1016/S1352-2310(97)00403-2."))

# 원본에서 미완성 절 제목은 빨간색(EE0000)으로 표시돼 있었다 — 채웠으므로 해제
for p in d.paragraphs:
    if p.style.name.startswith('Heading'):
        for c in list(p._p.iter(W + 'color')):
            if c.get(W + 'val') == 'EE0000':
                c.getparent().remove(c)

d.save(OUT)
print('저장:', OUT)
