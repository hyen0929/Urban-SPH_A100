from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

# ── 색상 팔레트 (기존 dry_deposition 스크립트와 동일) ────────────────
C_BG       = RGBColor(0x1E, 0x1E, 0x2E)
C_ACCENT   = RGBColor(0x89, 0xB4, 0xFA)
C_ACCENT2  = RGBColor(0xA6, 0xE3, 0xA1)
C_ACCENT3  = RGBColor(0xF3, 0x8B, 0xA8)
C_TEXT     = RGBColor(0xCD, 0xD6, 0xF4)
C_SUBTEXT  = RGBColor(0x93, 0x9A, 0xB0)
C_BOX      = RGBColor(0x31, 0x32, 0x44)
C_WHITE    = RGBColor(0xFF, 0xFF, 0xFF)
C_YELLOW   = RGBColor(0xF9, 0xE2, 0xAF)
C_ORANGE   = RGBColor(0xFA, 0xB3, 0x87)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]

TOTAL = 6

# ══════════════════════════════════════════════════════════
# 헬퍼 함수
# ══════════════════════════════════════════════════════════
def add_bg(slide, color=C_BG):
    w, h = prs.slide_width, prs.slide_height
    sh = slide.shapes.add_shape(1, 0, 0, w, h)
    sh.fill.solid(); sh.fill.fore_color.rgb = color; sh.line.fill.background()

def add_rect(slide, l, t, w, h, fill=C_BOX, line=None, line_color=C_ACCENT):
    sh = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    if line:
        sh.line.color.rgb = line_color; sh.line.width = Pt(line)
    else:
        sh.line.fill.background()
    return sh

def add_text(slide, text, l, t, w, h,
             size=18, bold=False, color=C_TEXT,
             align=PP_ALIGN.LEFT, italic=False):
    txb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    txb.word_wrap = True
    tf = txb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.alignment = align
    run = p.add_run(); run.text = text
    run.font.size = Pt(size); run.font.bold = bold
    run.font.italic = italic; run.font.color.rgb = color
    return txb

def hline(slide, l, t, w, color=C_SUBTEXT, thick=0.03):
    sh = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(thick))
    sh.fill.solid(); sh.fill.fore_color.rgb = color; sh.line.fill.background()

def accent_bar(slide, color=C_ACCENT):
    sh = slide.shapes.add_shape(1, 0, 0, prs.slide_width, Inches(0.12))
    sh.fill.solid(); sh.fill.fore_color.rgb = color; sh.line.fill.background()

def page_num(slide, n):
    add_text(slide, f"{n} / {TOTAL}", 12.3, 7.15, 0.8, 0.3,
             size=10, color=C_SUBTEXT, align=PP_ALIGN.RIGHT)


# ══════════════════════════════════════════════════════════
# Slide 1 — Title
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl)

# 상단 장식 바
for col, frac in [(C_ACCENT3, 0.45), (C_ACCENT2, 0.35), (C_ACCENT, 0.2)]:
    sh = sl.shapes.add_shape(1, 0, 0, int(prs.slide_width * frac), Inches(0.12))
    sh.fill.solid(); sh.fill.fore_color.rgb = col; sh.line.fill.background()

add_rect(sl, 1.2, 1.7, 10.9, 3.7, fill=C_BOX, line=1.5, line_color=C_ACCENT3)
add_text(sl, "방사성 입자 재부유 모델", 1.4, 1.9, 10.5, 1.1,
         size=36, bold=True, color=C_ACCENT3, align=PP_ALIGN.CENTER)
add_text(sl, "Stochastic Resuspension Model for Deposited Radioactive Particles", 1.4, 2.85, 10.5, 0.65,
         size=17, color=C_TEXT, align=PP_ALIGN.CENTER, italic=True)
add_text(sl, "Biasi et al. (2001)  ·  Loosmore (2003)  ·  SOPHIA_gpu LDM Module", 1.4, 3.58, 10.5, 0.55,
         size=14, color=C_SUBTEXT, align=PP_ALIGN.CENTER)
add_text(sl, "functions_DD.cuh — buffer_type 기반 확률적 침착/재부유 프레임워크", 1.4, 4.1, 10.5, 0.5,
         size=13, color=C_SUBTEXT, align=PP_ALIGN.CENTER)

add_rect(sl, 0, 6.8, 13.33, 0.7, fill=RGBColor(0x18, 0x18, 0x28))
add_text(sl, "ESLAB  ·  Seoul National University  ·  SOPHIA_gpu LDM Module", 0.3, 6.85, 12.7, 0.4,
         size=13, color=C_SUBTEXT, align=PP_ALIGN.CENTER)
page_num(sl, 1)


# ══════════════════════════════════════════════════════════
# Slide 2 — 재부유란 무엇인가 / 입자 상태 전이
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl); accent_bar(sl, C_ACCENT3)

add_text(sl, "재부유 (Resuspension) 개요", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT3)
add_text(sl, "Why Resuspension Matters for Radioactive Contamination", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 왼쪽: 개념 설명
add_rect(sl, 0.3, 1.35, 6.0, 5.7, fill=C_BOX, line=1, line_color=C_ACCENT3)
add_text(sl, "물리적 메커니즘", 0.5, 1.43, 5.6, 0.42, size=16, bold=True, color=C_ACCENT3)
hline(sl, 0.4, 1.9, 5.8, C_ACCENT3)

items = [
    ("풍속 의존성",    "마찰속도 u* 가 임계값 u*_t 초과 시 발생"),
    ("공기역학 리프팅","난류 버스트(turbulent burst)가 입자를 들어올림"),
    ("관련 표면",      "지면, 건물 지붕, 도로 등 건식 침착 표면"),
    ("방사성 세슘",    "Fukushima 이후 Cs-137 재이동 경로로 주목"),
    ("2차 오염",       "침착 후 재부유 → 인체 노출 경로 재활성화"),
]
for i, (title, desc) in enumerate(items):
    y = 2.0 + i * 0.88
    add_text(sl, title, 0.5, y, 2.0, 0.38, size=13, bold=True, color=C_ACCENT2)
    add_text(sl, desc,  2.55, y, 3.55, 0.75, size=12, color=C_TEXT)

# 오른쪽: buffer_type 상태 전이 다이어그램
add_rect(sl, 6.7, 1.35, 6.3, 5.7, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "입자 상태 전이 (buffer_type)", 6.9, 1.43, 5.9, 0.42,
         size=16, bold=True, color=C_ACCENT2)
hline(sl, 6.8, 1.9, 6.1, C_ACCENT2)

# 상태 박스 3개
states = [
    (C_SUBTEXT,  "i_type = 0",         "미방출 (Pre-release)",      "방출 대기 중"),
    (C_ACCENT2,  "buffer_type = 0",    "활성 (Active)",             "Lagrangian 추적 + 침착 체크"),
    (C_ACCENT3,  "buffer_type = 1",    "침착 (Deposited)",          "이동 없음, 재부유 체크"),
]
for i, (col, tag, label, desc) in enumerate(states):
    y = 2.1 + i * 1.55
    add_rect(sl, 7.0, y, 5.7, 1.2, fill=RGBColor(0x22, 0x22, 0x38), line=1.5, line_color=col)
    add_text(sl, tag,   7.15, y+0.07, 2.2, 0.38, size=11, color=C_SUBTEXT)
    add_text(sl, label, 7.15, y+0.4,  5.4, 0.38, size=15, bold=True, color=col)
    add_text(sl, desc,  7.15, y+0.76, 5.4, 0.35, size=11, color=C_TEXT)

# 화살표 텍스트
add_text(sl, "↓  침착 (P_dep)", 9.5, 3.32, 3.5, 0.38, size=12, color=C_ACCENT, bold=True)
add_text(sl, "↑  재부유 (P_res)", 9.5, 4.87, 3.5, 0.38, size=12, color=C_ORANGE, bold=True)

page_num(sl, 2)


# ══════════════════════════════════════════════════════════
# Slide 3 — 재부유율 모델 수식
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl); accent_bar(sl, C_ACCENT2)

add_text(sl, "재부유율 모델", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_ACCENT2)
add_text(sl, "Biasi et al. (2001)  ·  Loosmore (2003)  ·  First-order stochastic rate", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 핵심 수식 박스
add_rect(sl, 0.5, 1.32, 12.3, 0.95, fill=RGBColor(0x28, 0x30, 0x20), line=2, line_color=C_ACCENT2)
add_text(sl, "λ_res  =  C_res × max( 0,  u*²  −  u*_t² )    [s⁻¹]", 0.7, 1.42, 11.9, 0.75,
         size=26, bold=True, color=C_ACCENT2, align=PP_ALIGN.CENTER)

# 확률 식
add_rect(sl, 0.5, 2.42, 12.3, 0.85, fill=RGBColor(0x20, 0x28, 0x30), line=1.5, line_color=C_ACCENT)
add_text(sl, "P_res  =  1  −  exp( −λ_res × Δt )", 0.7, 2.52, 11.9, 0.65,
         size=22, bold=True, color=C_ACCENT, align=PP_ALIGN.CENTER)

# 파라미터 테이블
add_rect(sl, 0.4, 3.45, 8.2, 3.65, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "모델 파라미터", 0.6, 3.52, 7.8, 0.42, size=15, bold=True, color=C_ACCENT2)
hline(sl, 0.5, 3.97, 8.0, C_ACCENT2)

params = [
    ("C_res",    "1.0 × 10⁻³",  "s⁻¹/(m/s)²", "경험 계수 (Biasi 2001 기반)"),
    ("u*_t",     "0.2",          "m/s",         "재부유 임계 마찰속도"),
    ("u*",       "us_update",    "m/s",         "MOST 커널에서 실시간 갱신"),
    ("Δt",       "freq_LDM × dt","s",            "LDM 유효 타임스텝"),
    ("λ_res",    "0 ~ 10⁻³",    "s⁻¹",         "u* 에 따라 동적으로 계산"),
]
col_x = [0.6, 2.5, 4.3, 5.2]
for i, (sym, val, unit, note) in enumerate(params):
    y = 4.08 + i * 0.55
    add_text(sl, sym,  col_x[0], y, 1.8, 0.42, size=13, bold=True, color=C_YELLOW)
    add_text(sl, val,  col_x[1], y, 1.7, 0.42, size=13, color=C_TEXT)
    add_text(sl, unit, col_x[2], y, 0.9, 0.42, size=12, color=C_SUBTEXT)
    add_text(sl, note, col_x[3], y, 3.2, 0.42, size=12, color=C_TEXT)

# 오른쪽: u* vs λ 수치표
add_rect(sl, 8.9, 3.45, 4.1, 3.65, fill=C_BOX, line=1, line_color=C_ORANGE)
add_text(sl, "u* 별 λ_res 수치 예시", 9.1, 3.52, 3.7, 0.42, size=14, bold=True, color=C_ORANGE)
hline(sl, 9.0, 3.97, 3.9, C_ORANGE)

examples = [
    ("u* [m/s]", "λ_res [s⁻¹]", "반감기"),
    ("< 0.2",    "0",             "∞  (재부유 없음)"),
    ("0.3",      "5 × 10⁻⁵",    "~4 시간"),
    ("0.5",      "2.1 × 10⁻⁴",  "~55 분"),
    ("1.0",      "9.6 × 10⁻⁴",  "~12 분"),
]
for i, (a, b, c) in enumerate(examples):
    y = 4.08 + i * 0.55
    is_hdr = (i == 0)
    col_e = C_SUBTEXT if is_hdr else C_TEXT
    bold_e = is_hdr
    add_text(sl, a, 9.1,  y, 0.9,  0.42, size=12, bold=bold_e, color=col_e)
    add_text(sl, b, 10.1, y, 1.25, 0.42, size=12, bold=bold_e, color=col_e if is_hdr else C_ACCENT2)
    add_text(sl, c, 11.4, y, 1.4,  0.42, size=12, bold=bold_e, color=col_e if is_hdr else C_YELLOW)

page_num(sl, 3)


# ══════════════════════════════════════════════════════════
# Slide 4 — 모델 배경: Biasi & Loosmore
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl); accent_bar(sl, C_YELLOW)

add_text(sl, "모델 이론적 배경", 0.4, 0.2, 12.5, 0.6,
         size=28, bold=True, color=C_YELLOW)
add_text(sl, "Biasi et al. (2001)  &  Loosmore (2003) — 모델 비교 및 채택 근거", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 두 모델 박스 비교
for box_l, col, title, sub, lines in [
    (0.35, C_ACCENT2, "Biasi et al. (2001)",
     "Nuclear Safety / STORM experiment",
     [("대상",       "핵 시설 내 입자 재부유 실험 기반"),
      ("형태",       "λ = Λ₀ · (u*/u*₀)ⁿ"),
      ("특징",       "임계 마찰속도 개념 명확히 도입"),
      ("한계",       "실내 조건 기반 — 도시 야외 적용 시 보정 필요"),
      ("채택 이유",  "방사성 세슘 재부유에 가장 직접적으로 적용")]),
    (6.7, C_ORANGE, "Loosmore (2003)",
     "Outdoor resuspension / Long-term trends",
     [("대상",       "야외 오염 표면 장기 재부유 측정"),
      ("형태",       "λ = C · (u*² − u*_t²)   ← 본 구현에 채택"),
      ("특징",       "임계값 이하 완전 0: 물리적으로 명확"),
      ("장점",       "u*² 의존성 → 이차 풍속 민감도 반영"),
      ("채택 이유",  "SOPHIA_gpu 도시 흐름에서 u* 범위에 적합")]),
]:
    add_rect(sl, box_l, 1.35, 6.0, 5.6, fill=C_BOX, line=1.2, line_color=col)
    add_text(sl, title, box_l+0.2, 1.42, 5.6, 0.45, size=16, bold=True, color=col)
    add_text(sl, sub,   box_l+0.2, 1.85, 5.6, 0.35, size=11, color=C_SUBTEXT, italic=True)
    hline(sl, box_l+0.1, 2.23, 5.7, col)
    for i, (k, v) in enumerate(lines):
        y = 2.32 + i * 0.82
        add_text(sl, k+":", box_l+0.25, y,      1.3, 0.38, size=12, bold=True, color=col)
        add_text(sl, v,     box_l+1.65, y,      4.1, 0.7,  size=12, color=C_TEXT)

# 중앙 비교 표시
add_text(sl, "vs", 6.3, 3.7, 0.7, 0.5, size=22, bold=True, color=C_SUBTEXT, align=PP_ALIGN.CENTER)

page_num(sl, 4)


# ══════════════════════════════════════════════════════════
# Slide 5 — SOPHIA_gpu 구현 상세 (functions_DD.cuh)
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl); accent_bar(sl, C_ACCENT)

add_text(sl, "SOPHIA_gpu 구현: KERNEL_dry_deposition", 0.4, 0.2, 12.5, 0.6,
         size=26, bold=True, color=C_ACCENT)
add_text(sl, "functions_DD.cuh  —  buffer_type 기반 침착/재부유 판정", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# 코드 플로우 박스
add_rect(sl, 0.35, 1.32, 12.6, 5.25, fill=C_BOX, line=1, line_color=C_ACCENT)
add_text(sl, "KERNEL_dry_deposition (tdt_in, tt, LP1, z0_in)", 0.55, 1.38, 12.1, 0.42,
         size=14, bold=True, color=C_ACCENT)
hline(sl, 0.45, 1.83, 12.3, C_ACCENT)

code = [
    # (indent, text, color, bold)
    (0, "// ① RNG 초기화: Langevin step(offset=0)과 독립된 시퀀스(offset=+17)", C_SUBTEXT, False),
    (0, "curand_init(seed, i, tt*100 + 17, &ss);   r = curand_uniform(&ss);", C_TEXT, False),
    (0, "", C_TEXT, False),
    (0, "// ② 활성 입자 (buffer_type == 0) — 지면 침착 체크", C_SUBTEXT, False),
    (0, "if (buffer_type == 0  &&  dz < DD_h_dep)  {", C_ACCENT2, False),
    (1, "Vd    = DD_calc_Vd(z_particle, z0)     // Ra + Rb + Vs", C_TEXT, False),
    (1, "P_dep = DD_calc_P_dep(Vd, tdt_in, dz)  // 1 − exp(−Vd·dt/dz)", C_YELLOW, True),
    (1, "if (r < P_dep)  →  buffer_type = 1;  z = L1 + ε;  ux=uy=uz=0", C_ACCENT3, False),
    (0, "}", C_ACCENT2, False),
    (0, "", C_TEXT, False),
    (0, "// ③ 침착 입자 (buffer_type == 1) — 재부유 체크", C_SUBTEXT, False),
    (0, "if (buffer_type == 1)  {", C_ORANGE, False),
    (1, "λ_res = DD_calc_lambda_res(us_update)   // C_res·max(0, u*²−u*_t²)", C_TEXT, False),
    (1, "P_res = DD_calc_P_res(tdt_in)           // 1 − exp(−λ·dt)", C_YELLOW, True),
    (1, "if (r < P_res)  →  buffer_type = 0;  uz = us_update × 0.3", C_ACCENT2, False),
    (0, "}", C_ORANGE, False),
]
for i, (indent, text, col, bold) in enumerate(code):
    x = 0.55 + indent * 0.45
    add_text(sl, text, x, 1.95 + i * 0.295, 12.3 - indent*0.45, 0.32,
             size=11.5, color=col, bold=bold)

page_num(sl, 5)


# ══════════════════════════════════════════════════════════
# Slide 6 — 전체 침착-재부유 사이클 & 파라미터 요약
# ══════════════════════════════════════════════════════════
sl = prs.slides.add_slide(BLANK)
add_bg(sl); accent_bar(sl, C_ACCENT3)

add_text(sl, "전체 침착-재부유 사이클 & 파라미터 요약", 0.4, 0.2, 12.5, 0.6,
         size=26, bold=True, color=C_ACCENT3)
add_text(sl, "Dry Deposition–Resuspension Cycle  ·  Tuning Parameters", 0.4, 0.78, 12.5, 0.4,
         size=14, color=C_SUBTEXT, italic=True)

# ── 사이클 다이어그램 (왼쪽 절반) ──────────────────────────────────
add_rect(sl, 0.3, 1.32, 7.0, 5.7, fill=C_BOX, line=1, line_color=C_ACCENT3)
add_text(sl, "입자 생애 주기", 0.5, 1.4, 6.6, 0.42, size=15, bold=True, color=C_ACCENT3)
hline(sl, 0.4, 1.85, 6.8, C_ACCENT3)

# 상태 박스들
state_boxes = [
    (1.5, 2.1,  C_SUBTEXT, "방출 전\n(i_type = 0)",     "발생원 내 대기"),
    (1.5, 3.35, C_ACCENT2, "활성\n(buffer_type = 0)",   "Langevin 추적\n침착 체크 (P_dep)"),
    (1.5, 4.85, C_ACCENT3, "침착\n(buffer_type = 1)",   "지면 고정\n재부유 체크 (P_res)"),
]
for bx, by, col, label, desc in state_boxes:
    add_rect(sl, bx, by, 2.7, 0.95, fill=RGBColor(0x22, 0x22, 0x38), line=1.5, line_color=col)
    add_text(sl, label, bx+0.12, by+0.07, 1.3, 0.82, size=12, bold=True, color=col)
    add_text(sl, desc,  bx+1.5,  by+0.07, 1.1, 0.82, size=10, color=C_TEXT)

# 화살표 텍스트
add_text(sl, "↓  t_release 도달",                    4.35, 2.47, 2.8, 0.38, size=11, color=C_SUBTEXT)
add_text(sl, "↓  P_dep = 1−exp(−Vd·dt/dz)",         4.35, 3.72, 2.9, 0.38, size=11, color=C_ACCENT)
add_text(sl, "↑  P_res = 1−exp(−λ_res·dt)",         4.35, 4.55, 2.9, 0.38, size=11, color=C_ORANGE)
add_text(sl, "← 재부유 후 다시 활성 상태로",          1.2,  5.88, 4.2, 0.38, size=10, color=C_ORANGE, italic=True)

# ── 파라미터 요약 (오른쪽 절반) ─────────────────────────────────────
add_rect(sl, 7.6, 1.32, 5.45, 5.7, fill=C_BOX, line=1, line_color=C_ACCENT2)
add_text(sl, "튜닝 파라미터 (functions_DD.cuh)", 7.8, 1.4, 5.05, 0.42,
         size=14, bold=True, color=C_ACCENT2)
hline(sl, 7.7, 1.85, 5.25, C_ACCENT2)

tuning = [
    # (define, value, unit, description, color)
    ("DD_h_dep",        "2.0",        "m",               "침착 체크 높이 (지면 위)",   C_ACCENT),
    ("DD_ustar_thresh", "0.2",        "m/s",             "재부유 임계 마찰속도",        C_ORANGE),
    ("DD_C_res",        "1.0×10⁻³",  "s⁻¹/(m/s)²",     "재부유율 계수",               C_ACCENT2),
    ("DD_dp",           "1.0×10⁻⁶",  "m",               "입자 직경 (Cs-137)",          C_SUBTEXT),
    ("DD_tau_p",        "9.0×10⁻⁶",  "s",               "입자 이완시간",               C_SUBTEXT),
    ("DD_Vs_const",     "8.9×10⁻⁵",  "m/s",             "중력 침강속도",               C_SUBTEXT),
    ("DD_eps0",         "3.0",        "—",               "점성 경계층 경험 상수",        C_SUBTEXT),
]
for i, (name, val, unit, desc, col) in enumerate(tuning):
    y = 2.0 + i * 0.68
    add_rect(sl, 7.65, y, 5.3, 0.6, fill=RGBColor(0x22, 0x22, 0x38))
    add_text(sl, name, 7.78, y+0.07, 2.1, 0.44, size=11, bold=True, color=col)
    add_text(sl, val,  9.9,  y+0.07, 1.1, 0.44, size=11, color=C_YELLOW)
    add_text(sl, unit, 11.1, y+0.07, 1.0, 0.44, size=10, color=C_SUBTEXT)
    add_text(sl, desc, 7.78, y+0.38, 5.1, 0.25, size=9,  color=C_SUBTEXT)

page_num(sl, 6)


# ══════════════════════════════════════════════════════════
out = "/home/kdh/Urban-SPH/resuspension_model.pptx"
prs.save(out)
print(f"Saved: {out}")
