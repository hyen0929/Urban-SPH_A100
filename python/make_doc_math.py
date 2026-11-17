"""
OMML 수식을 직접 삽입하는 docx 생성 스크립트
"""
from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from lxml import etree
import copy

# ── OMML 네임스페이스 ──────────────────────────────────────────
OMML_NS = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
W_NS    = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'

def omml(tag, **attrs):
    return OxmlElement(f'm:{tag}')

def _e(tag, text=None, **kw):
    el = OxmlElement(tag)
    if text:
        el.text = text
    return el

# ── 수식 빌더 헬퍼 ────────────────────────────────────────────
def m(tag):           return OxmlElement(f'm:{tag}')
def mPr(parent):
    pr = m('rPr'); parent.append(pr); return pr

def mRun(text, italic=True, size_pt=None):
    """m:r 요소 생성"""
    r = m('r')
    rpr = m('rPr')
    if not italic:
        sty = m('sty'); sty.set(qn('m:val'), 'p'); rpr.append(sty)
    r.append(rpr)
    t = m('t'); t.text = text
    if text.startswith(' ') or text.endswith(' '):
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    r.append(t)
    return r

def mText(text):
    """일반 텍스트 run"""
    return mRun(text, italic=False)

def mSym(text):
    """이탤릭 수학 기호 run"""
    return mRun(text, italic=True)

def mFrac(num_els, den_els):
    """m:f (분수)"""
    f = m('f')
    num = m('num')
    for e in num_els: num.append(e)
    den = m('den')
    for e in den_els: den.append(e)
    f.append(num); f.append(den)
    return f

def mSup(base_els, sup_els):
    sSup = m('sSup')
    e = m('e')
    for x in base_els: e.append(x)
    sup = m('sup')
    for x in sup_els: sup.append(x)
    sSup.append(e); sSup.append(sup)
    return sSup

def mSub(base_els, sub_els):
    sSub = m('sSub')
    e = m('e')
    for x in base_els: e.append(x)
    sub = m('sub')
    for x in sub_els: sub.append(x)
    sSub.append(e); sSub.append(sub)
    return sSub

def mSubSup(base_els, sub_els, sup_els):
    sSubSup = m('sSubSup')
    e = m('e'); [e.append(x) for x in base_els]
    sub = m('sub'); [sub.append(x) for x in sub_els]
    sup = m('sup'); [sup.append(x) for x in sup_els]
    sSubSup.append(e); sSubSup.append(sub); sSubSup.append(sup)
    return sSubSup

def mNary(char, sub_els, sup_els, body_els, limLoc='undOvr'):
    """∑ 등 n-ary 연산자"""
    nary = m('nary')
    naryPr = m('naryPr')
    chrEl = m('chr'); chrEl.set(qn('m:val'), char)
    limLocEl = m('limLoc'); limLocEl.set(qn('m:val'), limLoc)
    naryPr.append(chrEl); naryPr.append(limLocEl)
    nary.append(naryPr)
    sub = m('sub'); [sub.append(x) for x in sub_els]
    sup = m('sup'); [sup.append(x) for x in sup_els]
    e = m('e'); [e.append(x) for x in body_els]
    nary.append(sub); nary.append(sup); nary.append(e)
    return nary

def mAbs(body_els):
    """|x| — m:d with |"""
    d = m('d')
    dPr = m('dPr')
    begChr = m('begChr'); begChr.set(qn('m:val'), '|')
    endChr = m('endChr'); endChr.set(qn('m:val'), '|')
    dPr.append(begChr); dPr.append(endChr)
    e = m('e'); [e.append(x) for x in body_els]
    d.append(dPr); d.append(e)
    return d

def mParen(body_els):
    d = m('d')
    e = m('e'); [e.append(x) for x in body_els]
    d.append(e)
    return d

def mSqrt(body_els):
    rad = m('rad')
    radPr = m('radPr')
    deg = m('deg')   # empty = sqrt
    e = m('e'); [e.append(x) for x in body_els]
    rad.append(radPr); rad.append(deg); rad.append(e)
    return rad

def oMath(display=True, *children):
    """oMath 블록 (display or inline)"""
    om = OxmlElement('m:oMath')
    if display:
        omPr = OxmlElement('m:oMathPara')
        omParaPr = OxmlElement('m:oMathParaPr')
        jc = OxmlElement('m:jc'); jc.set(qn('m:val'), 'center')
        omParaPr.append(jc)
        omPr.append(omParaPr)
        for c in children: om.append(c)
        omPr.append(om)
        return omPr
    else:
        for c in children: om.append(c)
        return om

def insertMath(para, display=True, *children):
    """단락에 수식 삽입"""
    om = OxmlElement('m:oMath')
    for c in children: om.append(c)
    if display:
        omPara = OxmlElement('m:oMathPara')
        omParaPr = OxmlElement('m:oMathParaPr')
        jc = OxmlElement('m:jc'); jc.set(qn('m:val'), 'center')
        omParaPr.append(jc)
        omPara.append(omParaPr)
        omPara.append(om)
        para._p.append(omPara)
    else:
        para._p.append(om)

# ── 문서 시작 ─────────────────────────────────────────────────
doc = Document()
for s in doc.sections:
    s.top_margin=Cm(2.5); s.bottom_margin=Cm(2.5)
    s.left_margin=Cm(3.0); s.right_margin=Cm(2.5)

def h1(text):
    p = doc.add_heading(text, level=1)
    p.runs[0].font.color.rgb = RGBColor(0x1F,0x49,0x7D)

def h2(text):
    p = doc.add_heading(text, level=2)
    if p.runs: p.runs[0].font.color.rgb = RGBColor(0x2E,0x6B,0x9E)

def body(text):
    p = doc.add_paragraph(text)
    p.paragraph_format.space_after = Pt(4)
    return p

def mathpara():
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after  = Pt(4)
    return p

# ════════════════════════════════════════════════════════
# 표지
# ════════════════════════════════════════════════════════
tp = doc.add_paragraph()
tp.paragraph_format.space_before = Pt(40)
from docx.enum.text import WD_ALIGN_PARAGRAPH
tp.alignment = WD_ALIGN_PARAGRAPH.CENTER
r=tp.add_run('Artificial Thermal Diffusion in SOPHIA_gpu')
r.bold=True; r.font.size=Pt(20); r.font.color.rgb=RGBColor(0x1F,0x49,0x7D)

sp=doc.add_paragraph(); sp.alignment=WD_ALIGN_PARAGRAPH.CENTER
r2=sp.add_run('수식 정리 및 코드 적용 위치  |  2026-04-13')
r2.font.size=Pt(12); r2.font.color.rgb=RGBColor(0x55,0x55,0x55)
doc.add_page_break()

# ════════════════════════════════════════════════════════
# 1. 온도 방정식 전체 구조
# ════════════════════════════════════════════════════════
h1('1. SPH 온도 방정식')
body('corrector step (KERNEL_time_update_projection) 에서의 온도 업데이트:')

# T_i^{n+1} = T_i^{n,0} + (dT/dt)_i · Δt
p = mathpara()
insertMath(p, True,
    mSubSup([mSym('T')],[mSym('i')],[mText('n+1')]),
    mText(' = '),
    mSubSup([mSym('T')],[mSym('i')],[mText('n,0')]),
    mText(' + '),
    mParen([
        mFrac([mSym('dT')],[mSym('dt')])
    ]),
    mSub([mText('')],[mSym('i')]),
    mText(' · Δt')
)

body('온도 변화율의 구성:')

# (dT/dt)_i = R_i / c_{p,i} − A_i^euler
p = mathpara()
insertMath(p, True,
    mParen([
        mFrac([mSym('dT')],[mSym('dt')])
    ]),
    mSub([mText('')],[mSym('i')]),
    mText(' = '),
    mFrac([mSub([mSym('R')],[mSym('i')])],
          [mSub([mSym('c')],[mText('p,i')])]),
    mText(' − '),
    mSup([mSub([mSym('A')],[mSym('i')])],[mText('euler')])
)

body('R_i 는 세 열 플럭스의 합:')

# R_i = R_cond + R_SGS + R_art
p = mathpara()
insertMath(p, True,
    mSub([mSym('R')],[mSym('i')]),
    mText(' = '),
    mSub([mSym('R')],[mText('cond')]),
    mText(' + '),
    mSub([mSym('R')],[mText('SGS')]),
    mText(' + '),
    mSub([mSym('R')],[mText('art')])
)

# ════════════════════════════════════════════════════════
# 2. 분자 열전도
# ════════════════════════════════════════════════════════
h2('1.1  분자 열전도  R_cond')

# R_cond = Σ_j  4 m_j* k_i k_j / ((k_i+k_j) ρ_i ρ_j*) (T_i−T_j) dW_ij/r_ij
p = mathpara()
insertMath(p, True,
    mSub([mSym('R')],[mText('cond')]),
    mText(' = '),
    mNary('∑',[mSym('j')],[],[
        mFrac(
            [mText('4'), mSup([mSym('m')],[mText('*')]),
             mSub([mSym('k')],[mSym('i')]),
             mSub([mSym('k')],[mSym('j')])],
            [mParen([
                mSub([mSym('k')],[mSym('i')]),
                mText('+'),
                mSub([mSym('k')],[mSym('j')])
             ]),
             mSub([mSym('ρ')],[mSym('i')]),
             mSup([mSub([mSym('ρ')],[mSym('j')])],[mText('*')])]
        ),
        mText(' ('),
        mSub([mSym('T')],[mSym('i')]),
        mText('−'),
        mSub([mSym('T')],[mSym('j')]),
        mText(') '),
        mFrac(
            [mSub([mText('dW')],[mSym('ij')])],
            [mSub([mSym('r')],[mSym('ij')])]
        )
    ])
)
body('※ IBM (p_type ≥ 1000) 및 내부 유체 (buffer_type = 99) 이웃 제외')
body('※ m*, ρ*: IBM 벽 이웃은 Solid_rho 기준으로 보정')

# ════════════════════════════════════════════════════════
# 3. SGS 열확산
# ════════════════════════════════════════════════════════
h2('1.2  SGS 난류 열확산  R_SGS  (Deardorff DDF)')

# R_SGS = Σ_j  m_j/(ρ_i ρ_j) · (ν_{t,i}+ν_{t,j})/Pr_t · (T_i−T_j) · dW_ij/r_ij
p = mathpara()
insertMath(p, True,
    mSub([mSym('R')],[mText('SGS')]),
    mText(' = '),
    mNary('∑',[mSym('j')],[],[
        mFrac([mSub([mSym('m')],[mSym('j')])],
              [mSub([mSym('ρ')],[mSym('i')]),mSub([mSym('ρ')],[mSym('j')])]),
        mText(' · '),
        mFrac(
            [mSub([mSym('ν')],[mText('t,i')]),mText('+'),mSub([mSym('ν')],[mText('t,j')])],
            [mSub([mText('Pr')],[mSym('t')])]
        ),
        mText(' · ('),
        mSub([mSym('T')],[mSym('i')]),
        mText('−'),
        mSub([mSym('T')],[mSym('j')]),
        mText(') · '),
        mFrac([mSub([mText('dW')],[mSym('ij')])],[mSub([mSym('r')],[mSym('ij')])])
    ])
)

body('여기서 SGS 와점도는 Deardorff 1-equation 모델:')

# ν_t = C_k · Δ · √k_SGS
p = mathpara()
insertMath(p, True,
    mSub([mSym('ν')],[mSym('t')]),
    mText(' = '),
    mSub([mSym('C')],[mSym('k')]),
    mText(' · Δ · '),
    mSqrt([mSub([mSym('k')],[mText('SGS')])])
)
body('C_k = 0.1,    Pr_t = 0.71,    Δ = max(σ_ij, h_i)')
body('※ buffer_type = 99 이웃 제외')

# ════════════════════════════════════════════════════════
# 4. Artificial thermal diffusion  ← 핵심
# ════════════════════════════════════════════════════════
h2('1.3  Artificial Thermal Diffusion  R_art  ★')

body('건물 코너 / kernel 불완전 영역의 온도 폭주를 방지하기 위해 추가된 항:')

# R_art = Σ_j  m_j/(ρ_i ρ_j) · α_art · (T_i − T_j) · dW_ij/r_ij
p = mathpara()
insertMath(p, True,
    mSub([mSym('R')],[mText('art')]),
    mText(' = '),
    mNary('∑',[mSym('j')],[],[
        mFrac([mSub([mSym('m')],[mSym('j')])],
              [mSub([mSym('ρ')],[mSym('i')]),mSub([mSym('ρ')],[mSym('j')])]),
        mText(' · '),
        mSub([mSym('α')],[mText('art')]),
        mText(' · ('),
        mSub([mSym('T')],[mSym('i')]),
        mText('−'),
        mSub([mSym('T')],[mSym('j')]),
        mText(') · '),
        mFrac([mSub([mText('dW')],[mSym('ij')])],[mSub([mSym('r')],[mSym('ij')])])
    ])
)

body('인공 열확산계수:')

# α_art = 0.5 · h_ij · |u_ij · x_ij| / r_ij
p = mathpara()
insertMath(p, True,
    mSub([mSym('α')],[mText('art')]),
    mText(' = 0.5 · '),
    mSub([mSym('h')],[mSym('ij')]),
    mText(' · '),
    mFrac(
        [mAbs([
            mParen([
                mSub([mSym('u')],[mSym('i')]),
                mText('−'),
                mSub([mSym('u')],[mSym('j')])])
            ,
            mText(' · '),
            mParen([
                mSub([mSym('x')],[mSym('i')]),
                mText('−'),
                mSub([mSym('x')],[mSym('j')])])
         ])],
        [mSub([mSym('r')],[mSym('ij')])]
    ),
    mText('    [m²/s]')
)

# h_ij = (h_i + h_j)/2
p = mathpara()
insertMath(p, True,
    mSub([mSym('h')],[mSym('ij')]),
    mText(' = '),
    mFrac(
        [mSub([mSym('h')],[mSym('i')]),mText('+'),mSub([mSym('h')],[mSym('j')])],
        [mText('2')]
    )
)

body('적용 조건 (다음 중 하나 이상):')

# 조건 1
p = mathpara()
insertMath(p, True,
    mText('조건 1 :  '),
    mParen([
        mSub([mSym('u')],[mSym('i')]),mText('−'),mSub([mSym('u')],[mSym('j')])
    ]),
    mText(' · '),
    mParen([
        mSub([mSym('x')],[mSym('i')]),mText('−'),mSub([mSym('x')],[mSym('j')])
    ]),
    mText(' < 0     (수렴 흐름: 압축성 진동 감쇠)')
)

# 조건 2
p = mathpara()
insertMath(p, True,
    mText('조건 2 :  '),
    mSub([mText('flt_s')],[mSym('i')]),
    mText(' < 0.7     (kernel 불완전 – IBM 벽 근방)')
)

body('부호 검증: dW/dr < 0 이므로, T_i > T_j 일 때 R_art < 0 → 입자 i 냉각 (올바른 방향)')

# ════════════════════════════════════════════════════════
# 5. Eulerian 이류
# ════════════════════════════════════════════════════════
h2('1.4  Eulerian 온도 이류  A_i^euler')

body('eli = 0 (Eulerian SPH) 입자에서만 활성:')

# A_i^euler = − Σ_j 2m_j/(ρ_i+ρ_j) (T_j−T_i) · raw_dot · (1−eli_i)
p = mathpara()
insertMath(p, True,
    mSup([mSub([mSym('A')],[mSym('i')])],[mText('euler')]),
    mText(' = −'),
    mNary('∑',[mSym('j')],[],[
        mFrac([mText('2'),mSub([mSym('m')],[mSym('j')])],
              [mSub([mSym('ρ')],[mSym('i')]),mText('+'),mSub([mSym('ρ')],[mSym('j')])]),
        mParen([
            mSub([mSym('T')],[mSym('j')]),mText('−'),mSub([mSym('T')],[mSym('i')])
        ]),
        mText(' · '),
        mSub([mText('raw_dot')],[mSym('ij')]),
        mText(' · (1−'),
        mSub([mText('eli')],[mSym('i')]),
        mText(')')
    ])
)

body('raw_dot 정의:')

# raw_dot = [(u_i−u_j)·(x_i−x_j) / r_ij] · dW/dr
p = mathpara()
insertMath(p, True,
    mSub([mText('raw_dot')],[mSym('ij')]),
    mText(' = '),
    mFrac(
        [mParen([
            mSub([mSym('u')],[mSym('i')]),mText('−'),mSub([mSym('u')],[mSym('j')])
         ]),
         mText(' · '),
         mParen([
            mSub([mSym('x')],[mSym('i')]),mText('−'),mSub([mSym('x')],[mSym('j')])
         ])],
        [mSub([mSym('r')],[mSym('ij')])]
    ),
    mText(' · '),
    mFrac([mSub([mText('dW')],[mSym('ij')])],[mSym('dr')])
)
body('※ KGC 보정 미적용 (코너에서 행렬 증폭 방지)  |  IBM, buffer_type=99 이웃 제외')

# ════════════════════════════════════════════════════════
# 6. kernel 충전도 flt_s
# ════════════════════════════════════════════════════════
h1('2. kernel 충전도 지표  flt_s')

# flt_s,i = Σ_j (m_j/ρ_j) W(|x_i−x_j|, h_i)
p = mathpara()
insertMath(p, True,
    mSub([mText('flt_s')],[mSym('i')]),
    mText(' = '),
    mNary('∑',[mSym('j')],[],[
        mFrac([mSub([mSym('m')],[mSym('j')])],[mSub([mSym('ρ')],[mSym('j')])]),
        mText(' · W('),
        mAbs([mSub([mSym('x')],[mSym('i')]),mText('−'),mSub([mSym('x')],[mSym('j')])]),
        mText(' , '),
        mSub([mSym('h')],[mSym('i')]),
        mText(')')
    ])
)
body('완전 지지(fully supported) 입자: flt_s ≈ 1')
body('IBM 벽 코너 입자: 이웃 절반 이상이 고체 측에 차단 → flt_s ≈ 0.4–0.7')
body('임계값 0.7: 벽면 1층 입자를 커버하는 경험적 값')

# ════════════════════════════════════════════════════════
# 7. 코드 위치 표
# ════════════════════════════════════════════════════════
h1('3. 코드 위치 (functions_PPE.cuh)')

tbl = doc.add_table(rows=1, cols=3)
tbl.style = 'Table Grid'
hdr = tbl.rows[0].cells
for c, t in zip(hdr, ['항목','위치 (approx. line)','설명']):
    c.text = t
    c.paragraphs[0].runs[0].bold = True

rows_data = [
    ('R_cond', '~L431–441', '분자 열전도, Fourier'),
    ('R_SGS',  '~L443–446', 'DDF SGS 열확산'),
    ('R_art',  '~L448–460', 'Artificial thermal diffusion\n조건: uij·xij<0  OR  flt_s<0.7'),
    ('A_euler','~L419–424', '온도 이류 (raw gradient)'),
    ('dtemp 누적','~L528–530','P3[i].dtemp += tmp_Rc/cp − eulert'),
    ('온도 업데이트','functions_TIME_ISPH.cuh ~L226','T_new = T0 + dtemp·Δt'),
]
for name, loc, desc in rows_data:
    row = tbl.add_row().cells
    row[0].text=name; row[1].text=loc; row[2].text=desc
    for c in row:
        c.paragraphs[0].runs[0].font.size = Pt(9.5)

doc.add_paragraph()

# ════════════════════════════════════════════════════════
# 8. 파라미터표
# ════════════════════════════════════════════════════════
h1('4. 관련 파라미터 (Parameters.cuh)')

tbl2 = doc.add_table(rows=1, cols=3)
tbl2.style='Table Grid'
hdr2=tbl2.rows[0].cells
for c,t in zip(hdr2,['기호','값','설명']):
    c.text=t; c.paragraphs[0].runs[0].bold=True

params=[
    ('Pr_t','0.71','Turbulent Prandtl number'),
    ('Air_k','0.025 W/mK','공기 분자 열전도율'),
    ('Air_cp','1005 J/kgK','공기 정압 비열'),
    ('C_k','0.1','DDF SGS TKE 확산 계수'),
    ('flt_s 임계값','0.7','kernel 불완전 판별 기준 (경험적)'),
    ('α_art 계수','0.5','인공 확산계수 스케일링 인자'),
]
for sym,val,desc in params:
    row=tbl2.add_row().cells
    row[0].text=sym; row[1].text=val; row[2].text=desc
    for c in row:
        c.paragraphs[0].runs[0].font.size=Pt(9.5)

out = '/home/kdh/Urban-SPH/artificial_thermal_diffusion.docx'
doc.save(out)
print(f'Saved: {out}')
