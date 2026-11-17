#!/usr/bin/env python3
"""설정 슬라이드(10번 Leipzig) — 비교 구성 블록 아래에 방출 물질 한 줄을 회색으로 덧붙인다.

'One city, three comparisons' 블록과 지면 경계 행 수정은 이미 발표자료에 반영돼 있다.
이 스크립트는 그 아래 회색 한 줄만 추가한다 — 발표자료 어디에도 무엇을 뿌렸는지가 없었다.
모델에 핵종은 들어가지 않는다(functions_DD.cuh 는 입경과 밀도만 쓴다). 그래서 에어로졸로 적는다.

  python3 make_slide_setup.py
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_slides as M
from pptx import Presentation

SRC = M.SRC
OUT = '/home/kdh/Urban-SPH/runs/NUTHOS15_setup_slide.pptx'
IDX = 9                                    # 0-기반 — 10번 "A real city geometry ... Leipzig"
LINE = 'Source term — a passive aerosol, with dry deposition and resuspension'


def main():
    prs = Presentation(SRC)
    slide = prs.slides[IDX]
    M.text(slide, 1.18, 6.52, 5.4, [LINE], size=8.5, color=M.GREY, bold=False, space=0)

    lst = prs.slides._sldIdLst
    for i, el in enumerate(list(lst)):
        if i != IDX:
            prs.part.drop_rel(el.get(M.RNS)); lst.remove(el)
    prs.save(OUT)
    print(f'저장: {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB, 슬라이드 {len(lst)}장)')


if __name__ == '__main__':
    main()
