from docx import Document
from docx.shared import Pt, RGBColor, Cm, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import copy

doc = Document()

# ── 기본 여백 설정 ──────────────────────────────────────────────
for section in doc.sections:
    section.top_margin    = Cm(2.5)
    section.bottom_margin = Cm(2.5)
    section.left_margin   = Cm(3.0)
    section.right_margin  = Cm(2.5)

# ── 스타일 헬퍼 ────────────────────────────────────────────────
def heading(doc, text, level=1):
    p = doc.add_heading(text, level=level)
    run = p.runs[0] if p.runs else p.add_run(text)
    run.font.color.rgb = RGBColor(0x1F, 0x49, 0x7D)
    return p

def body(doc, text):
    p = doc.add_paragraph(text)
    p.paragraph_format.space_after = Pt(4)
    for run in p.runs:
        run.font.size = Pt(11)
    return p

def italic_body(doc, text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.italic = True
    run.font.size = Pt(11)
    p.paragraph_format.space_after = Pt(4)
    return p

def math_block(doc, text):
    """수식 전용 단락 (들여쓰기 + 고정폭)"""
    p = doc.add_paragraph()
    p.paragraph_format.left_indent  = Cm(1.5)
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after  = Pt(2)
    run = p.add_run(text)
    run.font.name = 'Courier New'
    run.font.size = Pt(10.5)
    return p

def code_block(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent  = Cm(1.0)
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after  = Pt(2)
    run = p.add_run(text)
    run.font.name = 'Courier New'
    run.font.size = Pt(9.5)
    run.font.color.rgb = RGBColor(0x20, 0x5F, 0x20)
    return p

def label(doc, text):
    """수식 번호·라벨 오른쪽 정렬"""
    p = doc.add_paragraph(text)
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    p.paragraph_format.space_after = Pt(6)
    for run in p.runs:
        run.font.size = Pt(10)
    return p

def hr(doc):
    p = doc.add_paragraph()
    pPr = p._p.get_or_add_pPr()
    pBdr = OxmlElement('w:pBdr')
    bottom = OxmlElement('w:bottom')
    bottom.set(qn('w:val'), 'single')
    bottom.set(qn('w:sz'), '6')
    bottom.set(qn('w:space'), '1')
    bottom.set(qn('w:color'), 'AAAAAA')
    pBdr.append(bottom)
    pPr.append(pBdr)
    return p

# ══════════════════════════════════════════════════════════════
# 표지
# ══════════════════════════════════════════════════════════════
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
title.paragraph_format.space_before = Pt(40)
r = title.add_run('Artificial Thermal Diffusion in SOPHIA_gpu')
r.bold = True
r.font.size = Pt(20)
r.font.color.rgb = RGBColor(0x1F, 0x49, 0x7D)

sub = doc.add_paragraph()
sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
sub.paragraph_format.space_before = Pt(6)
r2 = sub.add_run('수식 정리 및 코드 적용 위치')
r2.font.size = Pt(13)
r2.font.color.rgb = RGBColor(0x44, 0x44, 0x44)

meta = doc.add_paragraph()
meta.alignment = WD_ALIGN_PARAGRAPH.CENTER
meta.paragraph_format.space_before = Pt(10)
r3 = meta.add_run('SOPHIA_gpu  ·  functions_PPE.cuh  ·  2026-04-13')
r3.font.size = Pt(10)
r3.font.color.rgb = RGBColor(0x88, 0x88, 0x88)

doc.add_page_break()

# ══════════════════════════════════════════════════════════════
# 1. 배경
# ══════════════════════════════════════════════════════════════
heading(doc, '1. 배경 및 문제 정의', level=1)

body(doc,
    'SOPHIA_gpu는 Lagrangian Eulerian 혼합 SPH 시간적분을 사용한다. '
    '온도 이류는 Eulerian SPH 형태(eli = 0)로 계산되며, '
    '건물(IBM) 벽면 코너와 같이 흐름 분리점(separation point)이 형성되는 '
    '영역에서 다음과 같은 양의 피드백(positive feedback) 메커니즘이 발생한다.')

body(doc,
    '① 코너 외부 유체 입자 i가 약간 가열됨 (Tᵢ > Tⱼ)\n'
    '② 분리 흐름: (uᵢ − uⱼ)·(xᵢ − xⱼ) > 0 (발산 흐름)\n'
    '③ dW/dr < 0  →  raw_dot = [(u_ij · x_ij) / r] · dW/dr < 0\n'
    '④ eulert = −2 Σ mⱼ/(ρᵢ+ρⱼ) · (Tⱼ−Tᵢ) · raw_dot > 0  →  입자 i 추가 가열\n'
    '⑤ 양의 피드백으로 온도가 수천 K까지 폭주')

body(doc,
    '또한 IBM 벽면 근방 입자는 kernel support가 불완전하여 '
    'SPH 합산(Σ mⱼ/ρⱼ · W_ij = flt_s)이 1보다 현저히 작다. '
    '이 조건에서는 단순 수렴/발산 구분만으로는 불충분하므로 '
    'flt_s < 0.7 조건을 추가 트리거로 적용한다.')

hr(doc)

# ══════════════════════════════════════════════════════════════
# 2. 온도 방정식 전체 구조
# ══════════════════════════════════════════════════════════════
heading(doc, '2. SPH 온도 방정식', level=1)

body(doc, '온도 업데이트 (corrector step, KERNEL_time_update_projection):')
math_block(doc, 'T_i^(n+1) = T_i^(n) + (dT/dt)_i · Δt')

body(doc, '여기서 온도 변화율 (dT/dt)_i 는 다음 세 항의 합으로 구성된다:')
math_block(doc,
    '(dT/dt)_i  =  (1/c_p,i) · R_i  −  A_i^euler')

doc.add_paragraph()

# ── 2.1 R_i (열전도 + SGS + artificial diffusion)
heading(doc, '2.1  열 플럭스 항  R_i', level=2)

body(doc,
    'R_i는 세 가지 열 플럭스의 합이다:')
math_block(doc,
    'R_i = R_cond + R_SGS + R_art')

# (a) 분자 열전도
body(doc, '(a) 분자 열전도 (Fourier):')
math_block(doc,
    'R_cond = Σⱼ  4 · (mⱼ* · kᵢ · kⱼ) / ((kᵢ + kⱼ) · ρᵢ · ρⱼ*) · (Tᵢ − Tⱼ) · (dW_ij / r_ij)')
math_block(doc,
    '  kᵢ = Air_k  (유체),  kᵢ = Solid_k  (벽면 p_type=0)\n'
    '  mⱼ*, ρⱼ* : IBM 벽 이웃의 경우 Solid_rho 기준으로 보정 적용')
body(doc, '※ IBM 입자(p_type ≥ 1000) 및 내부 유체(buffer_type = 99)는 이웃에서 제외')

doc.add_paragraph()

# (b) SGS 열확산
body(doc, '(b) SGS 난류 열확산 (Deardorff DDF 모델):')
math_block(doc,
    'R_SGS = Σⱼ  mⱼ/(ρᵢ·ρⱼ) · (ν_t,i + ν_t,j)/Pr_t · (Tᵢ − Tⱼ) · (dW_ij / r_ij)')
math_block(doc,
    '  ν_t = C_k · Δ · √k_SGS     (DDF 1-equation model)\n'
    '  Pr_t = 0.71                 (turbulent Prandtl number)\n'
    '  C_k  = 0.1,   Δ = max(σ_ij, hᵢ)')
body(doc, '※ buffer_type = 99 이웃 제외')

doc.add_paragraph()

# (c) artificial
body(doc, '(c) Artificial thermal diffusion (이번 적용 항):')
math_block(doc,
    'R_art = Σⱼ  mⱼ/(ρᵢ·ρⱼ) · α_art · (Tᵢ − Tⱼ) · (dW_ij / r_ij)')
math_block(doc,
    '  α_art = 0.5 · h_ij · |u_ij · x_ij| / r_ij    [m²/s]')
math_block(doc,
    '  h_ij  = (hᵢ + hⱼ) / 2')

body(doc, '적용 조건 (OR):')
math_block(doc,
    '  조건 1 : (u_ij · x_ij) < 0          (수렴 흐름, 압축성 진동 감쇠)\n'
    '  조건 2 : flt_s,i < 0.7              (kernel 불완전 – IBM 벽 근방)')
body(doc,
    '부호 검증: dW_ij/dr < 0이므로 Tᵢ > Tⱼ이면 R_art < 0 → 입자 i 냉각 (올바른 방향)')

hr(doc)

# ── 2.2 Eulerian advection
heading(doc, '2.2  Eulerian 이류 항  A_i^euler', level=2)

body(doc, '온도 이류 (Eulerian SPH, eli = 0 조건에서만 활성):')
math_block(doc,
    'A_i^euler = −Σⱼ  2mⱼ/(ρᵢ+ρⱼ) · (Tⱼ−Tᵢ) · raw_dot · (1 − eli_i)')
math_block(doc,
    '  raw_dot = [(uᵢ−uⱼ)·(xᵢ−xⱼ) / r_ij] · dW_ij/dr')
body(doc,
    '※ KGC 보정 미적용(raw gradient) — 코너에서 KGC 행렬 증폭 방지\n'
    '※ IBM 이웃(p_type ≥ 1000), 내부 유체(buffer_type = 99) 제외')

hr(doc)

# ══════════════════════════════════════════════════════════════
# 3. flt_s 정의
# ══════════════════════════════════════════════════════════════
heading(doc, '3. kernel 충전도 지표  flt_s', level=1)

body(doc,
    'flt_s는 입자 i 주변의 SPH kernel 합산으로 정의된다:')
math_block(doc,
    'flt_s,i = Σⱼ  (mⱼ/ρⱼ) · W(|xᵢ−xⱼ|, hᵢ)')
body(doc,
    '완전히 둘러싸인 입자: flt_s ≈ 1\n'
    'IBM 벽면 코너 입자: 이웃이 한쪽에만 존재 → flt_s ≪ 1\n'
    '임계값 0.7은 벽 근방 한 층 정도를 커버하도록 경험적으로 설정')

hr(doc)

# ══════════════════════════════════════════════════════════════
# 4. 코드 위치
# ══════════════════════════════════════════════════════════════
heading(doc, '4. 코드 적용 위치', level=1)

body(doc, '파일: functions_PPE.cuh  →  KERNEL_advection_force3D')
body(doc, '관련 항목과 위치:')

tbl = doc.add_table(rows=1, cols=3)
tbl.style = 'Table Grid'
hdr = tbl.rows[0].cells
hdr[0].text = '항목'
hdr[1].text = '코드 위치 (approximate line)'
hdr[2].text = '설명'
for cell in hdr:
    for run in cell.paragraphs[0].runs:
        run.bold = True
        run.font.size = Pt(10)

rows_data = [
    ('R_cond', '~L431–441', '분자 열전도\nbuffer_type!=99 및 p_type<1000 이웃만'),
    ('R_SGS', '~L443–446', 'DDF 난류 열확산\nbuffer_type!=99 이웃만'),
    ('R_art', '~L448–460', 'Artificial thermal diffusion\n조건: uij·xij<0  OR  flt_s<0.7'),
    ('A_euler', '~L419–424', '온도 이류\nraw gradient 사용, IBM·interior 제외'),
    ('dtemp 누적', '~L528–530', 'P3[i].dtemp += tmp_Rc/cp_i − eulert'),
    ('온도 업데이트', 'functions_TIME_ISPH.cuh ~L226',
     'T_new = T0 + dtemp · Δt  (buffer_type==0만)'),
]
for name, loc, desc in rows_data:
    row = tbl.add_row().cells
    row[0].text = name
    row[1].text = loc
    row[2].text = desc
    for c in row:
        for run in c.paragraphs[0].runs:
            run.font.size = Pt(9.5)

doc.add_paragraph()

hr(doc)

# ══════════════════════════════════════════════════════════════
# 5. 적용 영역 요약 (서술)
# ══════════════════════════════════════════════════════════════
heading(doc, '5. 적용 영역 및 물리적 의미', level=1)

body(doc,
    '아래 그림은 적용 영역을 개념적으로 나타낸다. '
    '건물 windward 코너를 중심으로 두 가지 조건에 의해 R_art가 활성화된다.')

body(doc, '[조건 1] 수렴 흐름  (uᵢⱼ · xᵢⱼ < 0)')
body(doc,
    '수렴 흐름 구역에서 SPH는 홀짝 진동(odd-even oscillation)이 발생하기 쉽다. '
    'Monaghan(1997) 인공 점성의 온도 아날로그로서 적용하여 이를 감쇠시킨다. '
    '일반 유체 영역에서도 활성화되나 α_art는 국소 속도 차에 비례하여 '
    '물리적 확산에 비해 작게 유지된다.')

body(doc, '[조건 2] kernel 불완전 영역  (flt_s < 0.7)')
body(doc,
    '건물 표면 바로 외부 입자는 이웃의 절반 이상이 IBM 고체 측에 가려져 '
    'kernel 적분이 불완전하다(flt_s ≈ 0.4–0.7). '
    '이 영역에서는 발산 흐름이라도 eulert의 부호 오류로 인한 온도 폭주가 '
    '발생할 수 있으므로 방향과 무관하게 확산을 강제한다. '
    '단, u_sig = 0이면 α_art = 0이어서 추가 확산 없음.')

body(doc,
    '결과적으로 벌크 흐름 영역(flt_s ≥ 0.7, 발산)에는 R_art가 적용되지 않아 '
    '일반 온도 이류 정확도를 유지한다.')

hr(doc)

# ══════════════════════════════════════════════════════════════
# 6. 파라미터 참조표
# ══════════════════════════════════════════════════════════════
heading(doc, '6. 관련 파라미터 (Parameters.cuh)', level=1)

tbl2 = doc.add_table(rows=1, cols=3)
tbl2.style = 'Table Grid'
hdr2 = tbl2.rows[0].cells
hdr2[0].text = '기호'
hdr2[1].text = '값'
hdr2[2].text = '설명'
for c in hdr2:
    for run in c.paragraphs[0].runs:
        run.bold = True
        run.font.size = Pt(10)

params = [
    ('Pr_t',       '0.71',    'Turbulent Prandtl number'),
    ('Air_k',      '0.025 W/mK', 'Air molecular thermal conductivity'),
    ('Air_cp',     '1005 J/kgK', 'Air specific heat at constant pressure'),
    ('C_k (DDF)',  '0.1',     'SGS TKE diffusion coefficient'),
    ('flt_s 임계', '0.7',     'Kernel 충전도 임계값 (경험적)'),
    ('α_art 계수', '0.5',     'Artificial diffusivity 스케일링 인자'),
]
for sym, val, desc in params:
    row = tbl2.add_row().cells
    row[0].text = sym
    row[1].text = val
    row[2].text = desc
    for c in row:
        for run in c.paragraphs[0].runs:
            run.font.size = Pt(9.5)

doc.add_paragraph()

hr(doc)

# ══════════════════════════════════════════════════════════════
# 저장
# ══════════════════════════════════════════════════════════════
out_path = '/home/kdh/Urban-SPH/artificial_thermal_diffusion.docx'
doc.save(out_path)
print(f'Saved: {out_path}')
