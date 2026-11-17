#!/usr/bin/env python3
"""안정 성층 도메인 노출 지도 — 중립 f25 와 같은 양식, 같은 색범위.

f25 를 왼쪽, 이 그림을 오른쪽에 나란히 놓으려는 것이므로 색범위를 f25 의 것
(= 중립 최대값)에 그대로 맞춘다. 그래야 두 장의 색이 같은 뜻을 갖는다.
안정 런에서 중립 최대를 넘는 셀은 그려지는 16,239 개 중 6 개(0.04 %)라 포화는 무시할 수준이다.

  python3 map_stable.py        →  figs/f41_domain_exposure_map_stable.{pdf,png}
"""
import os, sys
import numpy as np

sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st

HERE = os.path.dirname(os.path.abspath(__file__))
st.OUT = os.path.join(HERE, 'figs')          # figures 보다 먼저 — finish 가 이 값을 본다
import figures as F
from figstyle import finish

NEU = '/home/kdh/Urban-SPH/runs/lambdap/deposition/ped_grid_fine.npz'
STA = os.path.join(HERE, 'ped_grid_fine.npz')


def field(p):
    G = np.load(p); cell = float(G['cell'])
    return sum(G[f'S{i}'] for i in F.S)/(cell*cell*2.0), G['lim'], cell


EN, _, _      = field(NEU)
ES, lim, cell = field(STA)
gx = lim[0]+(np.arange(ES.shape[0])+.5)*cell
gy = lim[2]+(np.arange(ES.shape[1])+.5)*cell
one = 1.0/(cell*cell*2.0)                    # 1 입자·s 에 해당하는 밀도

fig, ax, cax = F._mapfig()
F._buildings(ax)
Em = np.ma.masked_where(ES < 0.5*one, ES)
m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=F.MAPCMAP,
                  vmin=np.log10(one), vmax=np.log10(EN.max()),      # ← f25 와 동일
                  shading='gouraud', rasterized=True, zorder=2)
cb = fig.colorbar(m, cax=cax)
cb.set_label(r'exposure below 2 m   ($\log_{10}$, particle s m$^{-3}$)', fontsize=7.5)
cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
for i, (sx, sy) in F.S.items():
    ax.plot(sx, sy, 'o', mfc='white', mec=st.C['ink'], ms=4.4, mew=.8, zorder=6)
F._mapdeco(ax)
F.label_points(ax, [F.S[i][0] for i in F.IDS], [F.S[i][1] for i in F.IDS],
               [str(i) for i in F.IDS], fontsize=7)
finish(fig, ax, 'f41_domain_exposure_map_stable', minor=False, tight=False)

n_over = int((ES > EN.max()).sum())
print(f'색범위 {np.log10(one):.2f} – {np.log10(EN.max()):.2f} (f25 와 동일), 포화 셀 {n_over} 개')
