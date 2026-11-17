"""
Generate temperature equations documentation (DOCX) for SOPHIA_gpu.
Reflects current state of the code (post unit-bug fix, 2026-04-25).
"""
from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


def add_heading(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    return h


def add_para(doc, text, bold=False, italic=False, mono=False, size=11):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if mono:
        run.font.name = 'Consolas'
        rPr = run._element.get_or_add_rPr()
        rFonts = OxmlElement('w:rFonts')
        rFonts.set(qn('w:ascii'), 'Consolas')
        rFonts.set(qn('w:hAnsi'), 'Consolas')
        rFonts.set(qn('w:cs'), 'Consolas')
        rPr.append(rFonts)
    return p


def add_code(doc, code):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Cm(0.5)
    run = p.add_run(code)
    run.font.size = Pt(9)
    run.font.name = 'Consolas'
    rPr = run._element.get_or_add_rPr()
    rFonts = OxmlElement('w:rFonts')
    rFonts.set(qn('w:ascii'), 'Consolas')
    rFonts.set(qn('w:hAnsi'), 'Consolas')
    rFonts.set(qn('w:cs'), 'Consolas')
    rPr.append(rFonts)
    shading = OxmlElement('w:shd')
    shading.set(qn('w:val'), 'clear')
    shading.set(qn('w:color'), 'auto')
    shading.set(qn('w:fill'), 'F2F2F2')
    p._p.get_or_add_pPr().append(shading)
    return p


def add_eq(doc, eq):
    """Mathematical equation rendered as bold serif italic text."""
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Cm(1.0)
    run = p.add_run(eq)
    run.font.size = Pt(12)
    run.font.name = 'Cambria Math'
    rPr = run._element.get_or_add_rPr()
    rFonts = OxmlElement('w:rFonts')
    rFonts.set(qn('w:ascii'), 'Cambria Math')
    rFonts.set(qn('w:hAnsi'), 'Cambria Math')
    rPr.append(rFonts)
    return p


# =========================================================================
doc = Document()
style = doc.styles['Normal']
style.font.name = 'Calibri'
style.font.size = Pt(11)

# Title
title = doc.add_heading('SOPHIA_gpu Temperature Equations (Current State)', level=0)
add_para(doc, '2026-04-25 — post unit-bug fix in SGS / artificial thermal diffusion',
         italic=True, size=10)

# -------------------------------------------------------------------------
add_heading(doc, '0. Overview', level=1)
add_para(doc,
    'Temperature evolution per fluid particle is computed as a sum of contributions '
    'accumulated into P3[i].dtemp, then time-integrated by predictor-corrector. The '
    'four contributing kernels (in launch order) are:')
add_para(doc, '  1. KERNEL_MOST_boundary3D — surface sensible heat flux (ground/MOST_buffer particles)')
add_para(doc, '  2. KERNEL_advection_force3D — IBM wall conduction + fluid–fluid Fourier conduction + SGS thermal diffusion + artificial thermal diffusion + advection')
add_para(doc, '  3. KERNEL_time_update_projection — explicit Euler time integration: temp = temp0 + dtemp · dt')
add_para(doc, '  4. KERNEL_zero_interior_velocity — resets interior-fluid (buffer_type==99) particles to T_ini each step')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '1. Time integration (predictor-corrector)', level=1)
add_para(doc, 'Mathematical form:')
add_eq(doc, 'T_i^{n+1} = T_i^{0} + (dT/dt)_i · Δt · 𝟙(buffer_type = 0)')
add_para(doc, 'Code (functions_TIME_ISPH.cuh):')
add_code(doc,
    'Real temp_new = P2[i].temp0 + P3[i].dtemp * t_dt * (buffer_type==0);\n'
    'TP1[i].temp = temp_new;')
add_para(doc,
    'Note: the (buffer_type==0) gate freezes inlet/outlet/MOST particles. The temp0 '
    'reference is the value at the start of the timestep stored in P2.')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '2. P3[i].dtemp accumulation', level=1)
add_para(doc, 'Final assembly at the end of KERNEL_advection_force3D (functions_PPE.cuh:549):')
add_code(doc, 'P3[i].dtemp += tmp_Rc/cpi + tmp_Rc_diff − eulert;')
add_para(doc, 'where:')
add_para(doc, '  • tmp_Rc — IBM wall + fluid Fourier conduction (W/kg form, divided by cp_i to get K/s)')
add_para(doc, '  • tmp_Rc_diff — SGS + artificial diffusion (already K/s, no cp division)')
add_para(doc, '  • eulert — Eulerian advection of temperature (K/s, subtracted)')
add_para(doc,
    'Additionally, MOST surface flux contributes to P3[i].dtemp directly in '
    'KERNEL_MOST_boundary3D BEFORE this kernel runs, so the kernel above accumulates '
    'with += rather than =.')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '3. Surface sensible heat flux (MOST)', level=1)
add_para(doc, 'Mathematical form (Monin–Obukhov similarity):')
add_eq(doc, 'H = −ρ · c_p · u* · θ*    [W/m²]')
add_eq(doc, 'u* = κ · |U_h| / Φ_m(ζ)')
add_eq(doc, 'θ* = κ · (T_i − θ_s) / Φ_h(ζ)')
add_eq(doc, '(dT/dt)_MOST = H / (Δz · ρ · c_p)    [K/s]')
add_para(doc, 'Code (functions_MOST.cuh:256–262):')
add_code(doc,
    'Real sensible_heat_flux = -rhoi*Air_cp*u_star*theta_star;       // W/m^2\n'
    'if(apply_heat) P3[i].dtemp += sensible_heat_flux/(space*rhoi*Air_cp);')
add_para(doc,
    'Δz = "space" (particle layer thickness). Applied only to MOST_buffer particles '
    '(apply_heat=true). Sets us_update / OB_length_global at the domain reference '
    'point for use by other kernels (e.g., dry deposition).')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '4. IBM wall conduction (Cleary–Monaghan harmonic, fluid ⇄ IBM)', level=1)
add_para(doc, 'Mathematical form:')
add_eq(doc,
    '(dE/dt)_i,IBM = Σ_{j∈IBM} m_j · 4·k_air·k_i / [(k_air + k_i) · ρ_i · ρ_j] · '
    '(T_i − T_j) · r_ij · ∇W_ij / r_ij²    [W/kg]')
add_para(doc,
    'IBM particles (p_type ≥ 1000) carry T = T_con. Both sides use Air_k '
    '(thin air-layer model — the wall transfers heat through a film of air, not '
    'through bulk solid). Search range uses IBM_length·h.')
add_para(doc, 'Code (functions_PPE.cuh:321–335):')
add_code(doc,
    'if(k_con_solve && P1[j].p_type>=1000){\n'
    '    Real hj_ibm = IBM_length * P1[j].h;\n'
    '    Real sum_ibm = 4.0f * mj_ibm * Air_k * kci * (tempi - tempj_ibm) * tdwij_ibm;\n'
    '    sum_ibm /= (tdist_ibm + 1e-10f) * rhoi * rhoj_ibm * (kci + Air_k);\n'
    '    tmp_Rc += sum_ibm;\n'
    '}')
add_para(doc, 'Sign: T_i > T_j and ∇W < 0 ⇒ sum_ibm < 0 ⇒ T_i cools toward wall ✓')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '5. Fluid–fluid Fourier conduction (Cleary–Monaghan)', level=1)
add_para(doc, 'Mathematical form:')
add_eq(doc,
    '(dE/dt)_i,fluid = Σ_{j∈fluid} m̃_j · 4·k_i·k_j / [(k_i+k_j) · ρ_i · ρ̃_j] · '
    '(T_i − T_j) · ∇W_ij / r_ij    [W/kg]')
add_para(doc,
    'Mass and density adjusted (m̃_j, ρ̃_j) so that solid–fluid mixed pairs use the '
    'solid-side ρ to maintain enthalpy continuity. IBM particles excluded here '
    '(handled separately in §4). Interior fluid (buffer_type==99) excluded.')
add_para(doc, 'Code (functions_PPE.cuh:448–459):')
add_code(doc,
    'if(k_con_solve && P1[j].buffer_type!=99){\n'
    '    kcj = Solid_k*(ptypej==0) + Air_k*(ptypej==1||ptypej==-1);\n'
    '    Real rhoj_temp, mj_temp;\n'
    '    if(ptypej==1||ptypej==-1) rhoj_temp = rhoj;\n'
    '    else rhoj_temp = Solid_rho;\n'
    '    mj_temp = mj * rhoj_temp/rhoj;\n'
    '    sum_con_H = 4.0*mj_temp*kcj*kci*(tempi-tempj)*tdwij;\n'
    '    sum_con_H /= (tdist+1e-10)*rhoi*rhoj_temp*(kci+kcj);\n'
    '    tmp_Rc += sum_con_H;\n'
    '}')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '6. SGS thermal diffusion (Deardorff turbulence model)', level=1)
add_para(doc, 'Mathematical form (Brookshaw Laplacian, K/s units):')
add_eq(doc,
    '(dT/dt)_i,SGS = Σ_{j∈fluid} m_j/ρ_j · (ν_t,i + ν_t,j)/Pr_t · '
    '(T_i − T_j) · ∇W_ij / r_ij    [K/s]')
add_para(doc,
    'Pr_t = 0.71 (turbulent Prandtl). ν_t = subgrid eddy viscosity from Deardorff '
    'one-equation model. Active only when k_turbulence_model == DDF. Interior fluid '
    '(buffer_type==99) excluded.')
add_para(doc, 'Code (functions_PPE.cuh:461–464):')
add_code(doc,
    'if(k_turbulence_model==DDF && P1[j].buffer_type!=99){\n'
    '    tmp_Rc_diff += mj/rhoj * ((P1[i].vis_t + P1[j].vis_t)/Pr_t)\n'
    '                 * (tempi - tempj) * tdwij/(tdist+1e-20);\n'
    '}')
add_para(doc,
    'Note: This term was previously written as m_j/(ρ_i·ρ_j)·… and accumulated into '
    'tmp_Rc, then divided by cp at the end. That gave units m³·K/(kg·s) and an '
    'extra /cp ⇒ effective diffusion rate was ~830× too small. Fixed 2026-04-25.',
    italic=True)
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '7. Artificial thermal diffusion (corner stabilization)', level=1)
add_para(doc, 'Mathematical form (Brookshaw Laplacian with single-particle α):')
add_eq(doc,
    '(dT/dt)_i,art = Σ_j 2·α_art · m_j/ρ_j · (T_i − T_j) · ∇W_ij / r_ij    [K/s]')
add_eq(doc,
    'α_art = 0.5 · h_ij · u_sig    [m²/s]')
add_eq(doc,
    'u_sig = |u_ij · r_ij| / r_ij,   h_ij = (h_i + h_j)/2')
add_para(doc, 'Activation gate:')
add_eq(doc, 'apply if (u_ij · r_ij < 0) ∨ (flt_s_i < 0.7)')
add_para(doc,
    'Converging flow always activates the term (Monaghan-style AV). For separation '
    'regions where u_ij·r_ij ≥ 0 but flt_s reflects a low-neighbor configuration '
    '(near IBM walls), it also activates to prevent eulert runaway.')
add_para(doc, 'Code (functions_PPE.cuh:465–478):')
add_code(doc,
    'if(k_con_solve){\n'
    '    Real uij_xij_t = (uxi-uxj)*(xi-xj) + (uyi-uyj)*(yi-yj) + (uzi-uzj)*(zi-zj);\n'
    '    if(uij_xij_t < 0 || P1[i].flt_s < 0.7f){\n'
    '        Real h_ij = (hi+hj)*0.5f;\n'
    '        Real u_sig = fabsf(uij_xij_t)/(tdist+1e-20f);\n'
    '        Real alpha_art = 0.5f*h_ij*u_sig;       // m^2/s\n'
    '        tmp_Rc_diff += 2.0f*mj/rhoj * alpha_art\n'
    '                     * (tempi-tempj) * tdwij/(tdist+1e-20f);\n'
    '    }\n'
    '}')
add_para(doc,
    'Factor 2 comes from the standard Brookshaw derivation: ∇²T ≈ 2·Σ m_j/ρ_j · '
    '(T_i−T_j) · (r·∇W)/r².')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '8. Eulerian advection of temperature', level=1)
add_para(doc, 'Mathematical form (subtractive — pair-symmetric advection of T):')
add_eq(doc,
    'eulert_i = Σ_{j∈fluid} 2 · m_j/(ρ_i+ρ_j) · (T_j − T_i) · '
    '[u_ij · r_ij] · ∇W_ij / r_ij · (1 − eli_i)    [K/s]')
add_para(doc,
    'eli is the ALE blending factor (1 = pure Lagrangian, 0 = pure Eulerian). For '
    'Eulerian-leaning particles the advective transport term is included. '
    'Interior fluid (buffer_type==99) excluded as both i and j.')
add_para(doc, 'Code (functions_PPE.cuh:439–442):')
add_code(doc,
    'if(k_con_solve && ptypej<1000 && P1[j].buffer_type!=99){\n'
    '    Real raw_dot_t = ((uxi-uxj)*(xi-xj) + (uyi-uyj)*(yi-yj) + (uzi-uzj)*(zi-zj))\n'
    '                   * tdwij/(tdist+1e-20f);\n'
    '    eulert -= 2.0*mj/(rhoi+rhoj)*(tempj-tempi)*raw_dot_t*(1-P1[i].eli);\n'
    '}')
add_para(doc, 'Final dtemp uses (− eulert) — see §2.')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '9. Interior fluid reset (buffer_type==99)', level=1)
add_para(doc, 'After projection step, particles flagged as inside-building have temperature pinned to T_ini:')
add_eq(doc, 'T_i ← T_ini    if buffer_type_i = 99')
add_para(doc, 'Code (functions_IBM.cuh:504–515):')
add_code(doc,
    '__global__ void KERNEL_zero_interior_velocity(part1* P1) {\n'
    '    if (P1[i].buffer_type != 99) return;\n'
    '    P1[i].ux = P1[i].uy = P1[i].uz = 0.0;\n'
    '    P1[i].pres = 0.0;\n'
    '    P1[i].temp = T_ini;   // prevent temperature accumulation inside buildings\n'
    '}')
add_para(doc,
    'This is a safety net — these particles are also excluded from all dtemp '
    'contributions above, and the time-integration gate (buffer_type==0) prevents '
    'them from updating. The reset is defensive against any leakage path.')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '10. Call site / launch order', level=1)
add_para(doc, 'Per timestep, kernels touching dtemp / temp run in this order (ISPH_Calc.cuh):')
add_code(doc,
    'cudaMemset(dev_P3, 0, sizeof(part3)*num_part2);   // dtemp = 0\n'
    '...\n'
    'KERNEL_MOST_boundary3D       (line ~168)  → P3.dtemp += MOST flux\n'
    'KERNEL_advection_force3D     (line ~182)  → P3.dtemp += conduction + diffusion − advection\n'
    '...\n'
    'KERNEL_time_update_projection (line ~256) → P1.temp = P2.temp0 + P3.dtemp·dt\n'
    'KERNEL_zero_interior_velocity (line ~224) → P1.temp = T_ini for buffer_type==99')
doc.add_paragraph()

# -------------------------------------------------------------------------
add_heading(doc, '11. Symbol reference', level=1)
add_para(doc, 'm_j        — particle mass [kg]', mono=True)
add_para(doc, 'ρ_i, ρ_j   — particle density [kg/m³]', mono=True)
add_para(doc, 'k_i, k_j   — thermal conductivity [W/(m·K)] (Air_k=0.025, Solid_k=1.0)', mono=True)
add_para(doc, 'cp_i       — specific heat [J/(kg·K)] (Air_cp=1005, Solid_cp=900)', mono=True)
add_para(doc, 'ν_t        — SGS eddy viscosity [m²/s]', mono=True)
add_para(doc, 'Pr_t       — turbulent Prandtl number (=0.71)', mono=True)
add_para(doc, 'r_ij       — |x_i − x_j| [m]', mono=True)
add_para(doc, '∇W_ij      — kernel gradient at i due to j [1/m^4]', mono=True)
add_para(doc, 'h_i, h_j   — smoothing length [m]', mono=True)
add_para(doc, 'eli        — ALE blending factor (0=Eulerian, 1=Lagrangian)', mono=True)
add_para(doc, 'flt_s      — Shepard filter (≈ Σ m_j/ρ_j · W_ij)', mono=True)
add_para(doc, 'u*, θ*     — friction velocity [m/s], temperature scale [K]', mono=True)
add_para(doc, 'Φ_m, Φ_h   — MOST stability functions', mono=True)
add_para(doc, 'T_ini      — initial reference temperature [K]', mono=True)
add_para(doc, 'T_con      — IBM wall (building surface) temperature [K]', mono=True)
add_para(doc, 'space      — particle spacing / layer thickness [m]', mono=True)

# -------------------------------------------------------------------------
out = '/home/kdh/Urban-SPH/temperature_equations.docx'
doc.save(out)
print(f'Saved: {out}')
