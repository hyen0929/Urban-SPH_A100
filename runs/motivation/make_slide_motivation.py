#!/usr/bin/env python3
"""도입부 슬라이드 두 판 — 그림만 다르고 본문 세 항목은 같다.

  A  NUTHOS15_motivation_ladder.pptx   거리 사다리(f50)를 위에, 항목 셋을 아래에
  B  NUTHOS15_motivation_maps.pptx     2패널 지도(f51)를 왼쪽, 항목 셋을 오른쪽에

  python3 make_slide_motivation.py
"""
import os, sys
from PIL import Image
from pptx.util import Inches

DECK = '/home/kdh/Urban-SPH/runs/lambdap/deposition'
sys.path.insert(0, DECK)
import make_slides as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, 'figs')
INK, RED, GREY = M.INK, M.RED, M.GREY

ITEMS = [
    ('01', 'The reactor moves next to the load',
     'District heating, hydrogen production and desalination place the i-SMR near consumption '
     'sites, not on a remote coast.'),
    ('02', 'The planning zone shrinks to the fence',
     'PAZ 3 – 5 km and UPZ 20 – 30 km were set for large reactors;  10 CFR 50.160 scales '
     'the plume zone down to the site boundary.'),
    ('03', 'The dispersion tool did not move with it',
     'The i-SMR EPZ methodology still takes a Gaussian plume as its primary tool, with building '
     'wakes entering only as a correction.'),
]


def pic(slide, name, left, top, width):
    p = os.path.join(FIGS, name)
    w, h = Image.open(p).size
    slide.shapes.add_picture(p, Inches(left), Inches(top), Inches(width), Inches(width*h/w))
    return width*h/w


def head(s):
    M.set_section(s, '01', 'MOTIVATION')
    M.set_title(s, [('The planning zone is shrinking ', False), ('onto the city', True)])


def items(s, x, y, w, step, num_pt, head_pt, body_pt, head_dy=0.30):
    for num, hd, body in ITEMS:
        M.text(s, x, y, 0.6, [[(num, RED, True, num_pt)]], space=0)
        M.text(s, x+0.62, y-0.02, w, [[(hd, INK, True, head_pt)]], space=0)
        M.text(s, x+0.62, y+head_dy, w, [[(body, GREY, False, body_pt)]], space=0)
        y += step


def slide_ladder(s):
    head(s)
    pic(s, 'f50_epz_ladder.png', 1.68, 1.70, 9.60)          # 높이 2.92
    M.rule(s, 1.18, 4.78, 10.35)
    items(s, 1.18, 4.92, 10.0, 0.62, 14, 13, 10.5, head_dy=0.29)


def slide_maps(s):
    head(s)
    pic(s, 'f51_epz_maps.png', 1.18, 1.85, 6.50)            # 높이 3.38
    items(s, 7.85, 1.98, 4.10, 1.12, 14, 12, 10, head_dy=0.27)
    M.text(s, 1.18, 5.66, 10.6,
           [[('At that range the transport is set by the street network — ', INK, False, 12),
             ('which a plume model does not carry.', RED, True, 12)]], space=0)


SHORT = [
    ('01', 'The reactor moves next to the load',
     'District heating, hydrogen and desalination put the i-SMR near consumption sites;  '
     'the nearest blocks stand a few hundred meters from the fence.'),
    ('02', 'Releases also start inside the city',
     'Radiological security events and transport accidents put the source at street level, '
     'inside the blocks.'),
    ('03', 'Buildings set the transport there',
     'Yet the i-SMR EPZ methodology still takes a Gaussian plume as its primary tool, with '
     'building wakes only as a correction.'),
]


def slide_city(s):
    M.set_section(s, '01', 'MOTIVATION')
    M.set_title(s, [('SMR siting puts the population ', False), ('in the near field', True)])
    pic(s, 'f52_smr_next_to_city.png', 1.10, 1.78, 6.85)     # 높이 4.36
    y = 1.95
    for num, hd, body in SHORT:
        M.text(s, 8.25, y, 0.6, [[(num, RED, True, 14)]], space=0)
        M.text(s, 8.87, y-0.02, 3.55, [[(hd, INK, True, 12)]], space=0)
        M.text(s, 8.87, y+0.27, 3.55, [[(body, GREY, False, 10)]], space=0)
        y += 1.25
    M.rule(s, 8.25, 5.44, 3.95)
    M.text(s, 8.25, 5.56, 3.95,
           [[('This work fixes the source term and varies only its position — ', INK, False, 11),
             ('13 releases in one city.', RED, True, 11)]], space=0)


PHOTO = '/home/kdh/Urban-SPH/runs/SMR near city.PNG'


def slide_photo(s):
    """실사 렌더를 크게 놓고 오른쪽에 설명 — 그림만 다르고 본문은 slide_city 와 같다."""
    M.set_section(s, '01', 'MOTIVATION')
    M.set_title(s, [('SMR siting puts the population ', False), ('in the near field', True)])

    from PIL import Image as _Im
    w, h = _Im.open(PHOTO).size
    W = 6.85
    s.shapes.add_picture(PHOTO, Inches(1.10), Inches(1.82), Inches(W), Inches(W*h/w))
    M.text(s, 1.10, 1.82 + W*h/w + 0.12, W,
           [[('An SMR at the edge of a city — the nearest residents a few hundred meters away.',
              INK, False, 11)]], space=0)
    M.text(s, 1.10, 1.82 + W*h/w + 0.40, W,
           ['Image: Westinghouse'], size=8, color=GREY, bold=False, space=0)

    y = 2.00
    for num, hd, body in SHORT:
        M.text(s, 8.25, y, 0.6, [[(num, RED, True, 14)]], space=0)
        M.text(s, 8.87, y-0.02, 3.55, [[(hd, INK, True, 12)]], space=0)
        M.text(s, 8.87, y+0.27, 3.55, [[(body, GREY, False, 10)]], space=0)
        y += 1.30
    M.rule(s, 8.25, 5.52, 3.95)
    M.text(s, 8.25, 5.64, 3.95,
           [[('This work fixes the source term and varies only its position — ', INK, False, 11),
             ('13 releases in one city.', RED, True, 11)]], space=0)


if __name__ == '__main__':
    M.build_single(slide_ladder, '/home/kdh/Urban-SPH/runs/NUTHOS15_motivation_ladder.pptx')
    M.build_single(slide_maps,   '/home/kdh/Urban-SPH/runs/NUTHOS15_motivation_maps.pptx')
    M.build_single(slide_city,   '/home/kdh/Urban-SPH/runs/NUTHOS15_motivation_city.pptx')
    M.build_single(slide_photo,  '/home/kdh/Urban-SPH/runs/NUTHOS15_motivation_photo.pptx')
