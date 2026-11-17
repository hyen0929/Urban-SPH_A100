#!/usr/bin/env python3
"""결론 슬라이드 한 장 — 안정도 항목(05)을 더하고 04 에 단서를 붙인다.

현재 발표자료의 18번 슬라이드를 그대로 두고 나머지를 지운 뒤,
  · 04 머리글 끝에 '— in neutral air' 를 붙이고
  · 같은 서식의 05 항목을 복제해 추가하고
  · 다섯 항목이 들어가도록 간격을 1.12 → 0.95 in 으로 좁힌다.
'behaviour' 는 발표자료 나머지 표기(미국식)에 맞춰 'behavior' 로 고친다.

  python3 make_slide_conclusion.py
"""
import copy, os, sys
from pptx import Presentation
from pptx.util import Inches

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_slides as M

SRC = '/home/kdh/Urban-SPH/runs/NUTHOS15_Dohyun Kim.pptx'
OUT = '/home/kdh/Urban-SPH/runs/NUTHOS15_conclusion_5items.pptx'
IDX = 17                                   # 0-기반 — 18번 Summary & Conclusion
TOPS = [2.02, 2.97, 3.92, 4.87, 5.82]      # 항목 다섯 개의 머리 위치 [in]
ITEMS = [('TextBox 87', 'TextBox 89', 'TextBox 91'),
         ('TextBox 93', 'TextBox 95', 'TextBox 97'),
         ('TextBox 99', 'TextBox 101', 'TextBox 103'),
         ('TextBox 105', 'TextBox 107', 'TextBox 109')]

HEAD04 = 'Urban form governs the spatial concentration of exposure — in neutral air.'
HEAD05 = 'Atmospheric stability changes both the magnitude and the form dependence.'
BODY05 = ('A 285 K ground beneath 300 K inflow raises pedestrian exposure by a factor of 2.98 '
          'across all thirteen releases — 1.10 from reduced deposition and 2.70 from vertical '
          'trapping, with the airborne median height falling from 20.2 to 10.2 m. The '
          'independence of total exposure from urban form does not survive: ρs turns from '
          '−0.01 to +0.62.')


def by_name(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    raise KeyError(name)


def set_text(shape, text):
    """첫 run 의 서식을 살린 채 글자만 바꾼다."""
    runs = shape.text_frame.paragraphs[0].runs
    runs[0].text = text
    for r in runs[1:]:
        r.text = ''


def place(slide, names, top):
    by_name(slide, names[0]).top = Inches(top)
    by_name(slide, names[1]).top = Inches(top - 0.02)
    by_name(slide, names[2]).top = Inches(top + 0.32)


def main():
    prs = Presentation(SRC)
    slide = prs.slides[IDX]

    set_text(by_name(slide, 'TextBox 107'), HEAD04)
    b = by_name(slide, 'TextBox 101')                 # 영국식 표기 정리
    set_text(b, b.text_frame.paragraphs[0].runs[0].text.replace('behaviour', 'behavior'))

    tree = slide.shapes._spTree                        # 04 항목을 복제해 05 를 만든다
    new = []
    for nm in ITEMS[3]:
        el = copy.deepcopy(by_name(slide, nm)._element)
        tree.append(el)
        new.append(slide.shapes[-1])
    new[0].name, new[1].name, new[2].name = 'TextBox 205', 'TextBox 207', 'TextBox 209'
    set_text(new[0], '05'); set_text(new[1], HEAD05); set_text(new[2], BODY05)

    for names, top in zip(ITEMS, TOPS):                # 간격 재배치
        place(slide, names, top)
    place(slide, ('TextBox 205', 'TextBox 207', 'TextBox 209'), TOPS[4])

    lst = prs.slides._sldIdLst
    for i, el in enumerate(list(lst)):
        if i != IDX:
            prs.part.drop_rel(el.get(M.RNS)); lst.remove(el)
    prs.save(OUT)
    print(f'저장: {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB, 슬라이드 {len(lst)}장)')


if __name__ == '__main__':
    main()
