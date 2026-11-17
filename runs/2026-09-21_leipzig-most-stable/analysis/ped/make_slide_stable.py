#!/usr/bin/env python3
"""안정 성층 결과 한 장짜리 슬라이드 — 중립(f25) 대 안정(f41) 나란히.

NUTHOS-15 발표자료의 크롬·서체·색을 그대로 쓰되 `make_slides.py` 는 건드리지 않는다 —
import 해서 헬퍼만 빌려 쓴다. 두 지도는 별개 그림이고 색범위가 같다(map_stable.py 참조).

  python3 make_slide_stable.py
"""
import os, sys
from PIL import Image
from pptx.util import Inches

DECK = '/home/kdh/Urban-SPH/runs/lambdap/deposition'
sys.path.insert(0, DECK)
import make_slides as M                       # SRC/크롬/서체/색을 그대로 물려받는다

HERE = os.path.dirname(os.path.abspath(__file__))
NEU_FIG = os.path.join(DECK, 'paper', 'figs', 'f25_domain_exposure_map.png')
STA_FIG = os.path.join(HERE, 'figs', 'f41_domain_exposure_map_stable.png')
OUT = '/home/kdh/Urban-SPH/runs/NUTHOS15_stable_exposure.pptx'

INK, RED = M.INK, M.RED


def pic(slide, path, left, top, width):
    """M.pic 과 같되 절대경로를 받는다 — 두 그림이 서로 다른 디렉터리에 있다."""
    w, h = Image.open(path).size
    slide.shapes.add_picture(path, Inches(left), Inches(top), Inches(width),
                             Inches(width*h/w))


def slide_stable(s):
    M.set_section(s, '03', 'RESULTS')
    M.set_title(s, [('Stable air keeps it at street level — ', False), ('3× exposure', True)])

    W, GAP = 5.125, 0.35                       # 두 지도 같은 크기, 본문 폭 10.6 에 맞춘다
    L2 = 1.18 + W + GAP
    M.text(s, 1.18, 1.72, W, ['Neutral   ·   13 releases, pedestrian level z < 2 m'], size=12)
    M.text(s, L2,   1.72, W, ['Stable  (Ts = 285 K)   ·   same color scale'], size=12)
    pic(s, NEU_FIG, 1.18, 2.02, W)
    pic(s, STA_FIG, L2,   2.02, W)

    M.rule(s, 1.18, 6.00, 10.35)

    b1 = [('Exposure below 2 m  ×2.98', RED, True, 12),
          (' — median ×2.85 per source (2.09 – 5.26), all 13 rise;  '
           'peak 1,460 → 4,823 particle·s·m⁻³', INK, False, 12)]
    b2 = [('Same footprint, denser inside', RED, True, 12),
          (' — exposed cells 16,252 → 16,239, yet the median cell rises ×1.65 and 57 % are higher',
           INK, False, 12)]
    b3 = [('Neutral is the optimistic case', RED, True, 12),
          (' — a neutral-only assessment understates the night-time dose threefold',
           INK, False, 12)]
    M.text(s, 1.18, 6.12, 10.6, [b1, b2, b3], space=1)


if __name__ == '__main__':
    M.build_single(slide_stable, OUT)
