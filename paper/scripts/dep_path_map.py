#!/usr/bin/env python3
"""f46 — 침적물이 어디에 붙었는가: 대비되는 두 방출점의 평면도.

제염 대상이 지면이냐 건물 외피냐가 방출 위치로 갈린다는 것을 보이는 그림.
`buffer_type`(1 = 지면, 2 = 벽·지붕)을 직접 출력한 벽면 침적 ON 계열
(runs/2026-10-07_leipzig-src50w-rot0, 50소스 × 풍향 0°)에서 침적 입자 위치를 그린다.

고른 두 곳은 양쪽 다 주변에 건물이 있다. 건물이 아예 없는 방출점(λp = 0)을 쓰면
"건물이 없으니 당연하다"가 되어 비교가 성립하지 않는다.
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
import figures as F
from figstyle import C, finish
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = '/home/kdh/Urban-SPH'
st.OUT = os.environ.get('FIG_OUT', os.path.join(HERE, '..', 'figs', 'main'))
RUN = f'{ROOT}/runs/2026-10-07_leipzig-src50w-rot0'
PANELS = [(18, 'Semi-open'), (6, 'Enclosed courtyard')]
HALF = 130.0

ns = {'__file__': f'{ROOT}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
p, Fd = ns['read_vtk'](f'{RUN}/out/LDM_30000stp.vtk')
src = Fd['p_type'].astype(int)-100
bt = Fd['buffer_type'].astype(int)
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(HERE, '..', 'data/00_setup/sources_n50.csv')))}
# 건물 외곽선은 IBM 마커 입자에서 1 m 격자로 직접 만든다. 공용 heightfield 는 2.5 m 라
# 대각 벽이 계단 모양으로 그려진다(실제 건물 윤곽이 아니다).
BP, _ = ns['read_vtk'](f'{RUN}/out/building_0stp.vtk')
BCMAP = LinearSegmentedColormap.from_list('b', [F.BLDG, F.BLDG])
BCELL = 1.0

def _shifts(g):
    return [np.roll(np.roll(g, dx, 0), dy, 1) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]

def footprint(x0, x1, y0, y1, cell=BCELL):
    """창 안 건물 평면 윤곽. 입자 래스터 → 닫기 → 구멍 채우기."""
    m = (BP[:, 0] >= x0) & (BP[:, 0] < x1) & (BP[:, 1] >= y0) & (BP[:, 1] < y1)
    nx, ny = int((x1-x0)/cell), int((y1-y0)/cell)
    g = np.zeros((nx, ny), bool)
    ix = ((BP[m, 0]-x0)/cell).astype(int).clip(0, nx-1)
    iy = ((BP[m, 1]-y0)/cell).astype(int).clip(0, ny-1)
    g[ix, iy] = True
    g = np.logical_and.reduce(_shifts(np.logical_or.reduce(_shifts(g))))      # closing
    out = np.zeros_like(g); out[0, :] = out[-1, :] = out[:, 0] = out[:, -1] = True
    out &= ~g
    for _ in range(4000):                                                     # 구멍 채우기
        n = np.logical_or.reduce(_shifts(out)) & ~g
        if n.sum() == out.sum(): break
        out = n
    g |= ~g & ~out
    return g, x0+np.arange(nx+1)*cell, y0+np.arange(ny+1)*cell

PW, GAP, L, R, B, T = 2.05, 0.05, 0.05, 0.05, 0.30, 0.34
FW, FH = L+2*PW+GAP+R, B+PW+T
fig = plt.figure(figsize=(FW, FH))
for k, (i, lab) in enumerate(PANELS):
    ax = fig.add_axes([(L+k*(PW+GAP))/FW, B/FH, PW/FW, PW/FH])
    sx, sy = S[i]
    fp, fx, fy = footprint(sx-HALF, sx+HALF, sy-HALF, sy+HALF)
    ax.pcolormesh(fx, fy, np.where(fp.T, 1.0, np.nan), cmap=BCMAP,
                  shading='flat', rasterized=True, zorder=1)
    m = src == i
    g, w = m & (bt == 1), m & (bt == 2)
    ax.plot(p[g, 0], p[g, 1], '.', color=C['blue'], ms=0.7, alpha=0.30,
            mec='none', rasterized=True, zorder=2)
    ax.plot(p[w, 0], p[w, 1], '.', color=C['orange'], ms=0.9, alpha=0.75,
            mec='none', rasterized=True, zorder=3)
    ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.0, mew=.8, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(sx-HALF, sx+HALF); ax.set_ylim(sy-HALF, sy+HALF)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_linewidth(.7)
    shr = w.sum()/max((g | w).sum(), 1)*100
    ax.set_title(f'{chr(97+k)}) {lab}', fontsize=8, pad=4)
    if k == 0:
        ax.annotate('', xy=(sx-45, sy+108), xytext=(sx-110, sy+108),
                    arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.0), zorder=7)
        ax.text(sx-110, sy+114, 'wind', fontsize=6.5, zorder=7)
    else:
        ax.plot([sx+20, sx+120], [sy-112, sy-112], '-', c=C['ink'], lw=1.6, zorder=7)
        ax.text(sx+70, sy-106, '100 m', fontsize=6.5, ha='center', zorder=7)
    print(f"  src {i}: 지면 {g.sum():,}  벽·지붕 {w.sum():,}  →  벽 비중 {shr:.1f} %")
LEGFS = 7.0
h = [plt.Line2D([], [], ls='none', marker='o', ms=3.4, mec='none', color=c, label=t)
     for c, t in ((C['blue'], 'deposited on the ground'), (C['orange'], 'deposited on walls and roofs'))]
# 글자 크기의 절반만큼 위로 (72 pt/in 기준, figure 높이 비율로 환산)
fig.legend(handles=h, loc='lower center',
           bbox_to_anchor=(0.5, 0.005 + (LEGFS/2)/(FH*72)), ncol=2, fontsize=LEGFS,
           frameon=False, handletextpad=0.35, columnspacing=1.4)
from matplotlib.transforms import Bbox
os.makedirs(st.OUT, exist_ok=True)
kw = dict(bbox_inches=Bbox([[0, 0], fig.get_size_inches()]), pad_inches=0, dpi=1200)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(st.OUT, f'f46_deposition_path_map.{ext}'), **kw)
plt.close(fig)
