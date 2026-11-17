from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import copy

# ── 색상 팔레트 ──────────────────────────────────────────
C_BG       = RGBColor(0x1E, 0x1E, 0x2E)   # 진한 네이비
C_ACCENT   = RGBColor(0x89, 0xB4, 0xFA)   # 연한 블루
C_ACCENT2  = RGBColor(0xA6, 0xE3, 0xA1)   # 연한 그린
C_ACCENT3  = RGBColor(0xF3, 0x8B, 0xA8)   # 연한 핑크
C_TEXT     = RGBColor(0xCD, 0xD6, 0xF4)   # 밝은 회색
C_SUBTEXT  = RGBColor(0x93, 0x9A, 0xB0)   # 어두운 회색
C_BOX      = RGBColor(0x31, 0x32, 0x44)   # 박스 배경
C_WHITE    = RGBColor(0xFF, 0xFF, 0xFF)
C_YELLOW   = RGBColor(0xF9, 0xE2, 0xAF)
C_ORANGE   = RGBColor(0xFA, 0xB3, 0x87)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)

BLANK = prs.slide_layouts[6]   # 완전 빈 레이아웃

# ══════════════════════════════════════════════════════════
# 헬퍼 함수들
# ══════════════════════════════════════════════════════════
def add_bg(slide, color=C_BG):
    """슬라이드 배경 사각형"""
    w, h = prs.slide_width, prs.slide_height
    sh = slide.shapes.add_shape(1, 0, 0, w, h)
    sh.fill.solid()
    sh.fill.fore_color.rgb = color
    sh.line.fill.background()
    return sh

def add_rect(slide, l, t, w, h, fill=C_BOX, line=None, line_color=C_ACCENT, radius=False):
    from pptx.util import Emu
    sh = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    if line:
        sh.line.color.rgb = line_color
        sh.line.width = Pt(line)
    else:
        sh.line.fill.background()
    return sh

def add_text(slide, text, l, t, w, h,
             size=18, bold=False, color=C_TEXT,
             align=PP_ALIGN.LEFT, italic=False, wrap=True):
    txb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    txb.word_wrap = wrap
    tf = txb.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txb

def add_multiline(slide, lines, l, t, w, h, size=16, color=C_TEXT, bold_first=False, line_spacing=None):
    """lines: list of (text, bold, color) tuples"""
    txb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    txb.word_wrap = True
    tf = txb.text_frame
    tf.word_wrap = True
    for i, item in enumerate(lines):
        if isinstance(item, str):
            text, bold, col = item, False, color
        else:
            text, bold, col = item
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        if line_spacing:
            from pptx.util import Pt as UPt
            from pptx.oxml.ns import qn
            from lxml import etree
            pPr = p._pPr if p._pPr is not None else p._p.get_or_add_pPr()
        run = p.add_run()
        run.text = text
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = col
    return txb

def add_accent_bar(slide, t=0.85, h=0.07, color=C_ACCENT):
    sh = slide.shapes.add_shape(1, 0, Inches(t), prs.slide_width, Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = color
    sh.line.fill.background()

def slide_number(slide, n, total):
    add_text(slide, f"{n} / {total}", 12.3, 7.15, 0.8, 0.3,
             size=10, color=C_SUBTEXT, align=PP_ALIGN.RIGHT)

# ══════════════════════════════════════════════════════════
# Slide 1 — Title
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)

# 상단 장식 바
for i, (col, frac) in enumerate([(C_ACCENT, 0.45), (C_ACCENT2, 0.35), (C_ACCENT3, 0.2)]):
    w_bar = prs.slide_width * frac
    sh = sl.shapes.add_shape(1, 0, 0, int(w_bar), Inches(0.12))
    sh.fill.solid(); sh.fill.fore_color.rgb = col; sh.line.fill.background()

# 중앙 제목 박스
add_rect(sl, 1.2, 1.8, 10.9, 3.5, fill=C_BOX, line=1.5, line_color=C_ACCENT)

add_text(sl, "방사성 물질 건식 침착 모델", 1.4, 2.0, 10.5, 1.2,
         size=36, bold=True, color=C_ACCENT, align=PP_ALIGN.CENTER)
add_text(sl, "Dry Deposition Model for Radioactive Particles", 1.4, 2.95, 10.5, 0.7,
         size=20, bold=False, color=C_TEXT, align=PP_ALIGN.CENTER, italic=True)
add_text(sl, "Ra + Rb Resistance Framework  ·  Lagrangian Stochastic Approach", 1.4, 3.65, 10.5, 0.6,
         size=16, color=C_SUBTEXT, align=PP_ALIGN.CENTER)

# 하단 정보
add_rect(sl, 0, 6.8, 13.33, 0.7, fill=RGBColor(0x18,0x18,0x28))
add_text(sl, "ESLAB  ·  Seoul National University  ·  SOPHIA_gpu LDM Module", 0.3, 6.85, 12.7, 0.4,
         size=13, color=C_SUBTEXT, align=PP_ALIGN.CENTER)

slide_number(sl, 1, 6)

# ══════════════════════════════════════════════════════════
# Slide 2 — 개요: 왜 침착 모델이 필요한가
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)
add_accent_bar(sl, 0.0, 0.12, C_ACCENT)

add_text(sl, "연구 배경 및 목적", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT)
add_text(sl, "Background & Motivation", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 왼쪽 박스: 현재 LDM
add_rect(sl, 0.3, 1.4, 5.8, 4.5, fill=C_BOX, line=1, line_color=C_ACCENT)
add_text(sl, "현재 SOPHIA LDM", 0.5, 1.5, 5.4, 0.5,
         size=17, bold=True, color=C_ACCENT)
items_l = [
    "• 라그랑지안 입자 추적",
    "• 난류 확산 (SGS 모델)",
    "• 입자 방출 (시간 지연 포함)",
    "• 농도 계산 (커널 보간)",
    "",
    "✗  지면 침착 없음",
    "✗  건물 표면 침착 없음",
    "✗  방사성 세슘 물리 없음",
]
for i, txt in enumerate(items_l):
    col = C_ACCENT3 if txt.startswith("✗") else C_TEXT
    add_text(sl, txt, 0.5, 2.05 + i*0.35, 5.4, 0.38, size=14, color=col)

# 오른쪽 박스: 필요한 것
add_rect(sl, 6.5, 1.4, 6.5, 4.5, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "추가 필요 모델", 6.7, 1.5, 6.1, 0.5,
         size=17, bold=True, color=C_ACCENT2)
items_r = [
    "• 저항 모델 (Ra + Rb)",
    "  → 공기역학적·점성 저항",
    "• 입자 물성 (세슘, d≈1 μm)",
    "  → Brownian 확산, 관성 충돌",
    "• 중력 침강속도 Vs",
    "",
    "✓  확률적 침착 판정",
    "✓  침착량 공간 분포 누적",
]
for i, txt in enumerate(items_r):
    col = C_ACCENT2 if txt.startswith("✓") else C_TEXT
    add_text(sl, txt, 6.7, 2.05 + i*0.35, 6.1, 0.38, size=14, color=col)

# 화살표 텍스트
add_text(sl, "→", 6.0, 3.4, 0.6, 0.5, size=28, bold=True, color=C_YELLOW, align=PP_ALIGN.CENTER)

slide_number(sl, 2, 6)

# ══════════════════════════════════════════════════════════
# Slide 3 — 저항 모델 전체 구조
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)
add_accent_bar(sl, 0.0, 0.12, C_ACCENT2)

add_text(sl, "저항 기반 침착 속도 모델", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT2)
add_text(sl, "Resistance Analogy  —  Wesely (1989)", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 핵심 식 박스
add_rect(sl, 0.5, 1.35, 12.3, 1.0, fill=RGBColor(0x28,0x30,0x48), line=2, line_color=C_YELLOW)
add_text(sl, "Vd  =  1 / (Ra + Rb)  +  Vs", 0.7, 1.45, 11.9, 0.75,
         size=30, bold=True, color=C_YELLOW, align=PP_ALIGN.CENTER)

# 세 항목 박스
boxes = [
    ("Ra — 공기역학적 저항", "Aerodynamic Resistance",
     ["ln(z/z₀) - ψₕ(z/L)", "Ra  =  ─────────────────", "         κ · u*", "",
      "• MOST 커널의 u*, L 재사용", "• 안정도 보정 ψₕ 포함"],
     C_ACCENT, 0.4),
    ("Rb — 점성 경계층 저항", "Viscous Sublayer Resistance",
     ["       1", "Rb  =  ─────────────────────", "     u* · ε₀ · (E_B + E_IM)", "",
      "ε₀ = 3.0  (경험 상수)", "E_B = Sc⁻²/³  (Brownian)",
      "E_IM = [St/(St+0.5)]²  (관성)"],
     C_ACCENT2, 4.83),
    ("Vs — 중력 침강속도", "Gravitational Settling",
     ["    ρₚ d²ₚ g Cс", "Vs = ───────────", "        18μ", "",
      "Cс ≈ 1.17  (Cunningham, 1μm)", "Vs ≈ 8.9×10⁻⁵ m/s  (Cs, 1μm)"],
     C_ACCENT3, 9.26),
]

for title, sub, lines, col, left in boxes:
    add_rect(sl, left, 2.55, 3.9, 4.5, fill=C_BOX, line=1.2, line_color=col)
    add_text(sl, title, left+0.15, 2.62, 3.6, 0.45, size=14, bold=True, color=col)
    add_text(sl, sub,   left+0.15, 3.02, 3.6, 0.35, size=11, color=C_SUBTEXT, italic=True)
    # 구분선
    ln = sl.shapes.add_shape(1, Inches(left+0.1), Inches(3.38), Inches(3.7), Inches(0.03))
    ln.fill.solid(); ln.fill.fore_color.rgb = col; ln.line.fill.background()
    for i, txt in enumerate(lines):
        mono = txt.startswith("R") or txt.startswith("V") or txt.startswith("E") or txt.startswith("ε") or "/" in txt or "=" in txt
        add_text(sl, txt, left+0.15, 3.45+i*0.4, 3.65, 0.42,
                 size=13, color=C_TEXT if not mono else C_WHITE,
                 bold=mono and ("=" in txt))

slide_number(sl, 3, 6)

# ══════════════════════════════════════════════════════════
# Slide 4 — Rb 상세: 세슘 입자 물성
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)
add_accent_bar(sl, 0.0, 0.12, C_ACCENT3)

add_text(sl, "Rb 상세 — 세슘 입자 물성", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT3)
add_text(sl, "Zhang et al. (2001)  ·  Slinn (1982)  ·  Seinfeld & Pandis (2006)", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 왼쪽: 세슘 입자 물성표
add_rect(sl, 0.3, 1.3, 5.5, 5.5, fill=C_BOX, line=1, line_color=C_ACCENT3)
add_text(sl, "Cs-137 입자 기본 물성", 0.5, 1.38, 5.1, 0.45,
         size=16, bold=True, color=C_ACCENT3)

props = [
    ("입자 직경 dp",    "~1 μm  (accumulation mode)"),
    ("입자 밀도 ρp",    "2,500 kg/m³"),
    ("Cunningham Cc",   "1.17  (dp = 1 μm)"),
    ("Brownian 확산 DB","~2.8 × 10⁻¹¹ m²/s"),
    ("Schmidt 수 Sc",   "ν / DB  ≈  5.4 × 10⁵"),
    ("침강속도 Vs",     "8.9 × 10⁻⁵ m/s"),
    ("입자 이완시간 τp","~9.0 × 10⁻⁶ s"),
]
for i, (k, v) in enumerate(props):
    y = 1.95 + i*0.55
    add_rect(sl, 0.35, y, 2.65, 0.48, fill=RGBColor(0x28,0x28,0x40))
    add_text(sl, k, 0.42, y+0.05, 2.5, 0.38, size=12, color=C_SUBTEXT)
    add_rect(sl, 3.02, y, 2.65, 0.48, fill=RGBColor(0x22,0x22,0x38))
    add_text(sl, v, 3.08, y+0.05, 2.55, 0.38, size=13, bold=True, color=C_TEXT)

# 오른쪽: 효율 메커니즘 설명
add_rect(sl, 6.2, 1.3, 6.8, 5.5, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "수집 효율 메커니즘", 6.4, 1.38, 6.4, 0.45,
         size=16, bold=True, color=C_ACCENT2)

mechs = [
    ("E_B  — Brownian 확산", C_ACCENT,
     ["E_B = Sc⁻²/³",
      "입자가 작을수록 우세 (< 0.1 μm)",
      "1 μm Cs : E_B ≈ 1.4 × 10⁻³",
      "→ Rb에서 지배적 항"]),
    ("E_IM  — 관성 충돌 (Impaction)", C_ACCENT2,
     ["E_IM = [St / (St + 0.5)]²",
      "St = τp u*² / ν  (Stokes 수)",
      "u*=0.5 m/s → St ≈ 0.15",
      "→ E_IM ≈ 0.053"]),
    ("E_IN  — 차단 (Interception)", C_SUBTEXT,
     ["E_IN = 0.5 (dp/A)²",
      "A: 표면 특성 길이",
      "세슘 1 μm에서 무시 가능"]),
]
y0 = 1.95
for title, col, lines in mechs:
    add_text(sl, title, 6.4, y0, 6.4, 0.38, size=14, bold=True, color=col)
    for j, ln in enumerate(lines):
        add_text(sl, ln, 6.55, y0+0.38+j*0.32, 6.3, 0.35, size=12, color=C_TEXT)
    y0 += 0.38 + len(lines)*0.32 + 0.25

slide_number(sl, 4, 6)

# ══════════════════════════════════════════════════════════
# Slide 5 — 확률적 침착 판정 알고리즘
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)
add_accent_bar(sl, 0.0, 0.12, C_YELLOW)

add_text(sl, "확률적 침착 판정", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_YELLOW)
add_text(sl, "Stochastic Deposition  —  Stohl et al. (2005)  FLEXPART", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 핵심 식
add_rect(sl, 0.5, 1.32, 12.3, 0.95, fill=RGBColor(0x2A,0x2A,0x20), line=2, line_color=C_YELLOW)
add_text(sl, "P_dep  =  1 − exp( −Vd · Δt / Δz_particle )", 0.7, 1.42, 11.9, 0.75,
         size=28, bold=True, color=C_YELLOW, align=PP_ALIGN.CENTER)

# 변수 설명
vars_desc = [
    ("Vd", "= 1/(Ra+Rb) + Vs",    "전체 침착 속도 [m/s]"),
    ("Δt", "= LDM 타임스텝",      "freq_LDM × dt [s]"),
    ("Δz", "= 입자↔표면 거리",    "지면: z_particle,  건물: IBM 입자까지 거리 [m]"),
]
add_rect(sl, 0.5, 2.45, 12.3, 1.5, fill=C_BOX)
for i, (var, eq, desc) in enumerate(vars_desc):
    y = 2.52 + i*0.45
    add_text(sl, var, 0.65, y, 0.6, 0.42, size=15, bold=True, color=C_ACCENT)
    add_text(sl, eq,  1.3,  y, 3.2, 0.42, size=14, color=C_TEXT)
    add_text(sl, "→  " + desc, 4.5, y, 8.1, 0.42, size=13, color=C_SUBTEXT)

# 알고리즘 플로우
add_rect(sl, 0.4, 4.1, 12.5, 3.05, fill=C_BOX, line=1, line_color=C_ACCENT)
add_text(sl, "KERNEL_time_update_LDM 내 판정 로직", 0.6, 4.15, 12.0, 0.42,
         size=15, bold=True, color=C_ACCENT)

code_lines = [
    ("if  z_particle < h_dep  (지면)  or  dist_IBM < h_IBM  (건물):",  C_ACCENT2, False),
    ("    Vd = calc_Vd(u_star, OB_length, dp, rho_p)        // Ra+Rb+Vs",  C_TEXT, False),
    ("    P_dep = 1 - exp(-Vd * dt / dz_particle)",                         C_YELLOW, True),
    ("    r = curand_uniform(&state)     // 균일 난수 생성",                 C_TEXT, False),
    ("    if  r < P_dep  →  LP1[i].i_type = 0    // 침착 (비활성화)",       C_ACCENT3, False),
    ("                    deposition_map[cell] += mass    // 침착량 누적",   C_ACCENT3, False),
]
for i, (line, col, bold) in enumerate(code_lines):
    add_text(sl, line, 0.6, 4.65+i*0.37, 12.3, 0.38, size=12.5, color=col, bold=bold)

slide_number(sl, 5, 6)

# ══════════════════════════════════════════════════════════
# Slide 6 — 참고문헌 & 구현 요약
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)
add_accent_bar(sl, 0.0, 0.12, C_ACCENT)

add_text(sl, "참고문헌 & 구현 계획", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT)
add_text(sl, "References & Implementation Plan", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 참고문헌
add_rect(sl, 0.3, 1.3, 7.5, 5.8, fill=C_BOX, line=1, line_color=C_ACCENT)
add_text(sl, "참고문헌", 0.5, 1.38, 7.1, 0.42, size=16, bold=True, color=C_ACCENT)

refs = [
    ("[1]  Wesely (1989)",
     "Parameterization of surface resistances to",
     "gaseous dry deposition.  Atmos. Environ. 23, 1293."),
    ("[2]  Zhang et al. (2001)",
     "A size-segregated particle dry deposition scheme",
     "for an atmos. aerosol module.  Atmos. Environ. 35, 549."),
    ("[3]  Slinn (1982)",
     "Predictions for particle deposition to vegetative",
     "canopies.  Atmos. Environ. 16, 1785."),
    ("[4]  Stohl et al. (2005)",
     "Technical note: The Lagrangian particle dispersion",
     "model FLEXPART v6.2.  Atmos. Chem. Phys. 5, 2461."),
    ("[5]  Seinfeld & Pandis (2006)",
     "Atmospheric Chemistry and Physics,",
     "2nd ed. Wiley.  Chapter 9."),
]
for i, (tag, l1, l2) in enumerate(refs):
    y = 1.92 + i*1.0
    add_text(sl, tag, 0.5,  y,      7.1, 0.36, size=13, bold=True, color=C_ACCENT2)
    add_text(sl, l1,  0.65, y+0.33, 6.9, 0.3,  size=11, color=C_TEXT)
    add_text(sl, l2,  0.65, y+0.60, 6.9, 0.3,  size=11, color=C_SUBTEXT)

# 구현 계획
add_rect(sl, 8.15, 1.3, 4.85, 5.8, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "구현 계획", 8.35, 1.38, 4.5, 0.42, size=16, bold=True, color=C_ACCENT2)

plan = [
    (C_ACCENT,  "Step 1", "calc_Vd() 함수 작성"),
    (C_SUBTEXT, "",       "Ra, Rb, Vs 계산"),
    (C_SUBTEXT, "",       "u*, L은 MOST에서 가져옴"),
    (C_ACCENT2, "Step 2", "지면 침착 판정 추가"),
    (C_SUBTEXT, "",       "z < h_dep 조건"),
    (C_SUBTEXT, "",       "cuRAND 확률 판정"),
    (C_ACCENT3, "Step 3", "IBM 표면 침착 추가"),
    (C_SUBTEXT, "",       "p_type≥1000 입자 탐색"),
    (C_SUBTEXT, "",       "dist < h_IBM 조건"),
    (C_YELLOW,  "Step 4", "침착량 공간 맵 출력"),
    (C_SUBTEXT, "",       "VTK deposition_map"),
]
for i, (col, step, desc) in enumerate(plan):
    y = 1.92 + i*0.46
    if step:
        add_text(sl, step, 8.35, y, 1.1, 0.4, size=13, bold=True, color=col)
        add_text(sl, desc, 9.5,  y, 3.3, 0.4, size=13, color=C_TEXT)
    else:
        add_text(sl, "· " + desc, 9.3, y, 3.5, 0.4, size=12, color=C_SUBTEXT)

slide_number(sl, 6, 6)

# ══════════════════════════════════════════════════════════
out = "/home/kdh/Urban-SPH/dry_deposition_model.pptx"
prs.save(out)
print(f"Saved: {out}")
