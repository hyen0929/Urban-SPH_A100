#!/usr/bin/env python3
"""NUTHOS-15 발표자료의 결과 슬라이드(11–13)와 시사점(14)을 새 결과로 교체한다.

원본은 건드리지 않고 `NUTHOS15_Dohyun Kim_v2.pptx` 로 저장한다.
크롬(대괄호·세로선·쪽번호·섹션 라벨)은 기존 14번 슬라이드를 복제해 그대로 물려받는다 —
그 슬라이드는 제목만 있고 그림·동영상 참조가 없어 rel 복사 없이 안전하게 복제된다.

  python3 make_slides.py
"""
import copy, os, sys
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_THEME_COLOR
from pptx.enum.text import PP_ALIGN

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'paper', 'figs')
MORPH = os.path.normpath(os.path.join(HERE, '..'))   # 형태 지표 CSV 위치
SRC = os.path.join(HERE, '..', '..', 'NUTHOS15_Dohyun Kim.pptx')
DST = os.path.join(HERE, '..', '..', 'NUTHOS15_Dohyun Kim_v2.pptx')

TITLE_FONT, BODY_FONT = 'Pretendard', 'Pretendard'
RED = RGBColor(0xC0, 0x00, 0x00)          # 기존 슬라이드의 강조색
GREY = RGBColor(0x7A, 0x8A, 0x92)         # 기존 패널 캡션 색
INK = RGBColor(0x2A, 0x2A, 0x2A)
RNS = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def clone_chrome(prs, src_idx):
    """크롬만 있는 슬라이드를 복제해 새 슬라이드로 돌려준다."""
    src = prs.slides[src_idx]
    dst = prs.slides.add_slide(src.slide_layout)
    for shp in list(dst.shapes):                      # 레이아웃 개체틀 제거
        shp._element.getparent().remove(shp._element)
    for shp in src.shapes:
        dst.shapes._spTree.append(copy.deepcopy(shp._element))
    return dst


def set_title(slide, parts):
    """제목 텍스트박스를 (문자열, 강조여부) 목록으로 다시 쓴다."""
    for sh in slide.shapes:
        if sh.name == 'TextBox 8':
            tf = sh.text_frame; tf.clear()
            pa = tf.paragraphs[0]
            for txt, hot in parts:
                r = pa.add_run(); r.text = txt
                r.font.name = TITLE_FONT; r.font.size = Pt(24); r.font.bold = True
                r.font.color.rgb = RED if hot else INK   # 흰 배경 슬라이드 → 본문색은 2A2A2A
            sh.width = Inches(11.0)
            return


def set_section(slide, num, label):
    for sh in slide.shapes:
        if sh.name == 'TextBox 1':
            tf = sh.text_frame; tf.clear(); pa = tf.paragraphs[0]
            a = pa.add_run(); a.text = num
            a.font.name = BODY_FONT; a.font.size = Pt(11); a.font.bold = True
            a.font.color.theme_color = MSO_THEME_COLOR.ACCENT_2
            b = pa.add_run(); b.text = '     '
            b.font.name = BODY_FONT; b.font.size = Pt(11); b.font.bold = True
            c = pa.add_run(); c.text = label
            c.font.name = BODY_FONT; c.font.size = Pt(10); c.font.bold = True
            c.font.color.theme_color = MSO_THEME_COLOR.BACKGROUND_2
            sh.width = Inches(2.0)
            return


def pic(slide, name, left, top, width):
    """폭을 주면 종횡비를 지켜 넣고, 실제 높이를 돌려준다."""
    p = os.path.join(FIGS, name)
    w, h = Image.open(p).size
    ht = width*h/w
    slide.shapes.add_picture(p, Inches(left), Inches(top), Inches(width), Inches(ht))
    return ht


def text(slide, left, top, width, lines, size=12.5, color=GREY, bold=True, space=6, font=None,
         align=None):
    """lines: 문자열 또는 (문자열, 색, 굵기, 크기) 튜플의 목록."""
    box = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(0.3))
    tf = box.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    first = True
    for ln in lines:
        pa = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        pa.space_after = Pt(space)
        if align is not None: pa.alignment = align
        segs = ln if isinstance(ln, list) else [(ln, color, bold, size)]
        for txt, col, bd, sz in segs:
            r = pa.add_run(); r.text = txt
            r.font.name = font or BODY_FONT; r.font.size = Pt(sz); r.font.bold = bd
            r.font.color.rgb = col
    return box


def rule(slide, left, top, width, color=RGBColor(0xD5, 0xDC, 0xE1)):
    from pptx.enum.shapes import MSO_SHAPE
    s = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top),
                               Inches(width), Emu(9525))
    s.fill.solid(); s.fill.fore_color.rgb = color; s.line.fill.background()
    s.shadow.inherit = False
    return s



# ---------------------------------------------------------------- 형태 지표 표
def _morph_rows():
    """분류별로 묶고, 그룹 안에서는 λp(150) 오름차순. (type, [(src, λ100, λ150, λ200, sector), ...])"""
    import csv as _csv
    lam = {}
    for r in _csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
        lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
    src = {}
    for r in _csv.DictReader(open(os.path.join(MORPH, 'morph_sources.csv'))):
        if int(float(r['R'])) == 150: src[int(r['id'])] = r
    def kind(i):
        if int(float(src[i]['court_closed'])): return 'enclosed courtyard'
        f = float(src[i]['free_sector_deg'])
        return 'open' if f >= 250 else ('semi-open' if f >= 100 else 'street canyon')
    groups = []
    for kd in ('open', 'semi-open', 'street canyon', 'enclosed courtyard'):
        mem = sorted([i for i in lam[150] if kind(i) == kd], key=lambda i: lam[150][i])
        groups.append((kd, [(i, lam[100][i], lam[150][i], lam[200][i],
                             float(src[i]['free_sector_deg'])) for i in mem]))
    return groups


def morph_table(s, left, top, width, height):
    groups = _morph_rows()
    n = sum(len(g[1]) for g in groups)
    hdr = ['Type', 'Source', 'λp\n(100 m)', 'λp\n(150 m)', 'λp\n(200 m)', 'Open\nsector']
    shp = s.shapes.add_table(n+1, len(hdr), Inches(left), Inches(top),
                             Inches(width), Inches(height))
    tbl = shp.table
    tbl.first_row = True; tbl.horz_banding = False
    for w, c in zip((1.58, 0.70, 0.70, 0.72, 0.70, 0.60), tbl.columns):
        c.width = Inches(w*width/5.0)
    for r in tbl.rows: r.height = Inches(height/(n+1))

    def put(cell, txt, bold=False, size=8.5, color=INK, align=PP_ALIGN.CENTER, fill=None):
        cell.text = ''
        tf = cell.text_frame; tf.word_wrap = False
        tf.margin_left = tf.margin_right = Inches(0.03)
        tf.margin_top = tf.margin_bottom = 0
        for k, line in enumerate(txt.split('\n')):
            pa = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
            pa.alignment = align
            run = pa.add_run(); run.text = line
            run.font.name = BODY_FONT; run.font.size = Pt(size); run.font.bold = bold
            run.font.color.rgb = color
        cell.fill.solid()
        cell.fill.fore_color.rgb = fill if fill else RGBColor(0xFF, 0xFF, 0xFF)
        cell.vertical_anchor = 3                                  # MSO_ANCHOR.MIDDLE

    HEAD = RGBColor(0xE8, 0xED, 0xF1); PICK = RGBColor(0xFB, 0xEC, 0xEC)
    TINT = [RGBColor(0xFF, 0xFF, 0xFF), RGBColor(0xF7, 0xF9, 0xFB)]   # 그룹 교대 음영
    for j, h in enumerate(hdr):
        put(tbl.cell(0, j), h, bold=True, size=8, fill=HEAD)
    r = 1
    for gi, (kd, mem) in enumerate(groups):
        r0 = r
        for (i, a, b, c, fs) in mem:
            bg = TINT[gi % 2]
            put(tbl.cell(r, 1), f'S{i}', bold=True, fill=bg)
            put(tbl.cell(r, 2), f'{a:.2f}', fill=bg)
            put(tbl.cell(r, 3), f'{b:.2f}', bold=True, fill=PICK)  # 채택한 반경
            put(tbl.cell(r, 4), f'{c:.2f}', fill=bg)
            put(tbl.cell(r, 5), f'{fs:.0f}°', fill=bg)
            r += 1
        cell = tbl.cell(r0, 0)
        if r-1 > r0: cell.merge(tbl.cell(r-1, 0))                  # 세로 병합(제자리 연산)
        put(cell, kd, bold=True, size=8.5, align=PP_ALIGN.LEFT, fill=TINT[gi % 2])
    return shp


# ---------------------------------------------------------------- 슬라이드 내용
def slide_method(s):
    set_section(s, '03', 'RESULTS')
    set_title(s, [('Quantifying urban form: ', False), ('plan area density', True)])
    text(s, 1.18, 1.64, 10.6,
         [[('λ', INK, True, 15), ('p', INK, True, 11),
           (' (x, y ; R)   =   building plan area within radius R   ÷   window area',
            INK, False, 15)]], align=PP_ALIGN.CENTER, space=0)
    text(s, 1.30, 2.18, 4.3, [r'λp field   (R = 150 m)'], size=12)
    pic(s, 'f31_lambda_field.png', 1.35, 2.44, 4.28)
    text(s, 6.30, 2.18, 5.4, ['Each release point, three window radii'], size=12)
    morph_table(s, 6.30, 2.44, 5.45, 3.18)
    rule(s, 1.18, 5.98, 10.35)
    b1 = [('Ranking is stable across the band', RED, True, 12),
          (' — rank correlation with λp(150) is 0.91 (100 m), 0.97 (200 m);  middle value used',
           INK, False, 12)]
    b2 = [('Open sector', RED, True, 12),
          (' — 360 rays are cast from the release point;  the widest angular gap that meets '
           'no building within 50 m', INK, False, 12)]
    b3 = [('Enclosed courtyards', RED, True, 12),
          (' — the open area around the release never reaches the street network;  '
           'kept out of the fits', INK, False, 12)]
    text(s, 1.18, 6.10, 10.6, [b1, b2, b3], space=1)


def slide11(s):
    set_section(s, '03', 'RESULTS')
    set_title(s, [('Same source term, ', False), ('four dispersion regimes', True)])
    text(s, 1.22, 1.72, 5.2, ['All 13 releases   ·   pedestrian level, z < 2 m'], size=12.5)
    pic(s, 'f25_domain_exposure_map.png', 1.22, 2.02, 5.25)
    text(s, 6.95, 1.72, 4.9, ['One release of each type, same colour scale'], size=12.5)
    pic(s, 'f34_closeups_by_type.png', 6.95, 2.02, 4.90)
    rule(s, 1.18, 6.14, 10.35)
    num = [('Footprint ', INK, False, 12), ('A', INK, True, 12), ('₉₀', INK, True, 9),
           ('   1,000 – 14,900 m²', RED, True, 12),
           ('   (15× among street-level sources; 600 m² in the enclosed courtyard)', INK, False, 11)]
    peak = [('Peak pedestrian exposure', INK, False, 12),
            ('   105 – 903 particle·s·m⁻³   (8.6×)', RED, True, 12),
            (',  trapped fraction 0.80 – 1.00', INK, False, 12)]
    text(s, 1.18, 6.26, 10.4, [num, peak], space=2)
    text(s, 1.18, 6.94, 10.4,
         ['13 releases in one city — identical rate, schedule and random seed (99,988 particles each).'],
         size=10.5, bold=False)


def slide12(s):
    set_section(s, '03', 'RESULTS')
    set_title(s, [('And it is predictable — ', False), ('form sets the footprint', True)])
    text(s, 1.18, 1.68, 10.6,
         [[('Deposition footprint against plan area density  λp (R = 150 m),  '
            '13 releases', INK, False, 13)]], align=PP_ALIGN.CENTER, space=0)
    text(s, 1.55, 2.20, 4.6, ['How far the material spreads'], size=12)
    pic(s, 'f01_footprint_vs_lambda.png', 1.55, 2.50, 4.75)
    text(s, 6.95, 2.20, 4.6, ['How much stays within 150 m'], size=12)
    pic(s, 'f02_trapped_vs_lambda.png', 6.95, 2.50, 4.75)
    rule(s, 1.18, 6.00, 10.35)
    b1 = [('Fitted over the 10 street-level sources', RED, True, 12),
          (' — the three enclosed courtyards are shown but excluded', INK, False, 12)]
    b2 = [('ρs = −0.84', RED, True, 12), (' footprint area,  ', INK, False, 12),
          ('+0.92', RED, True, 12), (' trapped fraction   (p < 0.05)', INK, False, 12)]
    b3 = [('Invariant', RED, True, 12),
          (' to release height (2.5 / 5 / 20 m) and to time step (0.005 / 0.01 s)',
           INK, False, 12)]
    text(s, 1.18, 6.12, 10.6, [b1, b2, b3], space=1)


def slide13(s):
    set_section(s, '03', 'RESULTS')
    set_title(s, [('Averaging keeps the plume, loses the ', False),
                  ('block-scale pattern', True)])
    text(s, 1.30, 1.72, 5.4, ['Transient LES  −  frozen mean flow   (log₁₀)'], size=12.5)
    pic(s, 'f08_depmap_difference.png', 1.30, 2.02, 4.60)
    text(s, 7.05, 1.72, 4.8, ['How far down in scale the agreement survives'], size=12.5)
    pic(s, 'f28_meanflow_scale.png', 7.05, 2.02, 4.55)
    rule(s, 1.18, 6.02, 10.35)
    b1 = [('Bulk reproduced', INK, False, 12.5),
          ('   deposited amount ×1.08 and ×0.99', RED, True, 12.5),
          (' — Pearson on counts stays 0.93–0.99 at every scale', INK, False, 12.5)]
    b2 = [('Pattern lost below the building scale', INK, False, 12.5),
          ('   rank correlation 0.93 → 0.45', RED, True, 12.5),
          (',  both geometries on one curve', INK, False, 12.5)]
    b3 = [('Footprint size off by 10–20 %', INK, False, 12.5),
          (',  sign set by the geometry   (×0.80 city,  ×1.09 array)', INK, False, 12.5)]
    text(s, 1.18, 6.16, 10.6, [b1, b2, b3], space=2)
    text(s, 1.18, 7.00, 10.6,
         ['Identical release and seeds; only the flow differs.  Sub-grid Langevin dispersion retained in both.'],
         size=10.5, bold=False)


def slide_diff(s):
    """전이 대 평균장 — 도메인 전체(차이 지도)와 분류별 국소(클로즈업)를 한 장에."""
    set_section(s, '03', 'RESULTS')
    set_title(s, [('Freeze the wind and the plume ', False), ('changes shape', True)])
    text(s, 1.18, 1.68, 10.6, [[('Pedestrian exposure below 2 m,  transient LES vs its own 150 s mean',
                                 INK, False, 13)]], align=PP_ALIGN.CENTER, space=0)
    text(s, 1.25, 2.34, 3.7, ['All 13 releases   ·   difference'], size=11.5)
    pic(s, 'f36_exposure_difference.png', 1.25, 2.64, 3.55)
    text(s, 5.15, 2.04, 6.9, ['Representative release of each type'], size=11.5)
    pic(s, 'f35_transient_vs_mean_closeups.png', 5.15, 2.34, 6.95)
    rule(s, 1.18, 5.58, 10.35)
    b1 = [('Bulk is reproduced', RED, True, 12),
          (' — total exposure ×1.03, Pearson 0.97 on cells;  but rank correlation is only 0.63',
           INK, False, 12)]
    b2 = [('The most-changed release of each type', RED, True, 12),
          (' falls to rank ρ = 0.23 – 0.52;  these four are also the four lowest of all 13 '
           '(range 0.23 – 0.92)', INK, False, 12)]
    b3 = [('Best agreement is the tightest courtyard', RED, True, 12),
          (' (ρ = 0.92) — averaging is right where the air is already trapped', INK, False, 12)]
    text(s, 1.18, 5.70, 10.6, [b1, b2, b3], space=1)
    text(s, 1.18, 6.40, 10.6,
         ['Identical release and seeds;  only the flow differs — transient LES vs its own 150 s mean.'],
         size=10.5, bold=False)


def slide14(s):
    set_section(s, '04', 'IMPLICATIONS')
    set_title(s, [('What this changes', False)])
    items = [
        ('01', 'Release position and building geometry dominate near-field consequence.',
         'Thirteen releases of identical strength, schedule and random seed differ only in position. '
         'The deposition footprint area varies by a factor of 15 and the peak pedestrian-level exposure '
         'by a factor of 8.6.'),
        ('02', 'Plan area density predicts the form of deposition and dispersion.',
         'λp evaluated over a 150 m window yields ρs = −0.84 for footprint area and +0.92 for the fraction '
         'retained within 150 m (p ≤ 0.004). The relation is preserved across release heights of 2.5 – 20 m '
         'and across time steps of 0.005 – 0.01 s.'),
        ('03', 'A time-averaged flow field alters the behaviour of the dispersing particles.',
         'With the release and the random seed held fixed, substituting the resolved flow by its own 150 s '
         'mean reduces the cell-wise rank correlation of pedestrian exposure to 0.63, and to 0.23 for '
         'individual releases. The footprint contracts to 80 % and the mean plume height falls from '
         '20.2 to 17.1 m.'),
        ('04', 'Urban form governs the spatial concentration of exposure, not its total.',
         'Total pedestrian exposure is independent of urban form (ρs = −0.12) while the peak concentration '
         'varies by a factor of 8.6; the deposited mass is likewise not predictable from form. '
         'For dose assessment the peak is the governing quantity.'),
    ]
    y = 2.02
    for num, head, body in items:
        text(s, 1.18, y, 0.6, [[(num, RED, True, 15)]], space=0)
        text(s, 1.88, y - 0.02, 10.0, [[(head, INK, True, 14)]], space=0)
        text(s, 1.88, y + 0.32, 10.0, [[(body, GREY, False, 11)]], space=0)
        y += 1.12


# ---------------------------------------------------------------- 조립
def main():
    prs = Presentation(SRC)
    n_before = len(prs.slides._sldIdLst)
    print(f'원본 슬라이드 {n_before}장')

    built = []
    for fn in (slide_method, slide11, slide12, slide13, slide14):
        s = clone_chrome(prs, 13)            # 14번(제목만 있는 시사점 슬라이드)을 템플릿으로
        fn(s); built.append(s)
    print(f'새 슬라이드 {len(built)}장 생성')

    lst = prs.slides._sldIdLst
    ids = list(lst)
    new_ids = ids[n_before:]                 # 방금 덧붙인 5장 (method, 11, 12, 13, 14)
    for i in (13, 12, 11, 10):               # 기존 11~14 제거(0-기반), 뒤에서부터
        el = ids[i]
        prs.part.drop_rel(el.get(RNS))
        lst.remove(el)
    # 형태 지표 슬라이드는 설정(9번) 바로 뒤 = 10번 자리.
    # 나머지 결과 슬라이드는 기존 10번(resolved turbulence) 다음인 11번 자리부터.
    meth, rest = new_ids[0], new_ids[1:]
    lst.remove(meth); lst.insert(9, meth)
    for k, el in enumerate(rest):
        lst.remove(el); lst.insert(11 + k, el)

    prs.save(DST)
    print(f'저장: {DST}  ({os.path.getsize(DST)/1e6:.0f} MB, 슬라이드 {len(prs.slides._sldIdLst)}장)')


def build_single(fn, out, template_idx=13):
    """슬라이드 한 장만 담은 pptx — 원본에서 크롬만 물려받고 나머지 슬라이드는 지운다."""
    prs = Presentation(SRC)
    n0 = len(prs.slides._sldIdLst)
    fn(clone_chrome(prs, template_idx))
    lst = prs.slides._sldIdLst
    for el in list(lst)[:n0]:                        # 기존 슬라이드 전부 제거
        prs.part.drop_rel(el.get(RNS)); lst.remove(el)
    prs.save(out)
    print(f'저장: {out}  ({os.path.getsize(out)/1e6:.1f} MB, 슬라이드 {len(lst)}장)')


if __name__ == '__main__':
    if '--only' in sys.argv:                      # 한 장만: --only slide12 [파일명]
        fn = globals()[sys.argv[sys.argv.index('--only')+1]]
        out = sys.argv[sys.argv.index('--only')+2] if len(sys.argv) > sys.argv.index('--only')+2 \
            else 'NUTHOS15_one_slide.pptx'
        build_single(fn, os.path.join(HERE, '..', '..', out))
    elif '--method-only' in sys.argv:
        build_single(slide_method, os.path.join(HERE, '..', '..', 'NUTHOS15_method_slide.pptx'))
    elif '--diff-only' in sys.argv:
        build_single(slide_diff, os.path.join(HERE, '..', '..', 'NUTHOS15_transient_vs_mean.pptx'))
    else:
        main()
