#!/usr/bin/env python3
"""중립 대 안정 성층 보행자 노출 — 평면도 비교.

  f39  나란히 (중립 | 안정), 색 범위 공유
  f40  한 장 차이 지도  log10(안정 / 중립)

격자·건물·바람화살표·축 규격은 f25/f34/f36 과 동일하다 (figures.py 헬퍼 재사용).
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
if os.environ.get('FIG_OUT'): st.OUT = os.environ['FIG_OUT']
import figures as F
from figstyle import C, finish

HERE = os.path.dirname(os.path.abspath(__file__))
NEU = '/home/kdh/Urban-SPH/runs/lambdap/deposition/ped_grid_fine.npz'
STA = os.path.join(HERE, 'ped_grid_fine.npz')

def field(p):
    G = np.load(p); cell = float(G['cell'])
    return sum(G[f'S{i}'] for i in F.S)/(cell*cell*2.0), G['lim'], cell

EN, lim, cell = field(NEU)
ES, _, _      = field(STA)
one = 1.0/(cell*cell*2.0)                       # 1 입자·s 에 해당하는 밀도
gx = lim[0]+(np.arange(EN.shape[0])+.5)*cell
gy = lim[2]+(np.arange(EN.shape[1])+.5)*cell


def _wind(ax):
    ax.annotate('', xy=(360, 420), xytext=(240, 420),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.1), zorder=7)
    ax.text(240, 432, 'wind', fontsize=7.5, zorder=7)

def _scale(ax):
    ax.plot([950, 1150], [-440, -440], '-', c=C['ink'], lw=1.8, zorder=7)
    ax.text(1050, -432, '200 m', fontsize=7.5, ha='center', zorder=7)


# ---------------------------------------------------------------- f39 나란히
PW, PH = 3.30, 3.20                              # 패널 [in] (도메인 종횡비 960/930)
L, GAP, CBG, CBW, R, B, T = 0.58, 0.26, 0.13, 0.16, 0.60, 0.50, 0.26
FW = L + PW + GAP + PW + CBG + CBW + R
FH = B + PH + T
fig = plt.figure(figsize=(FW, FH))
axes = [fig.add_axes([(L+k*(PW+GAP))/FW, B/FH, PW/FW, PH/FH]) for k in (0, 1)]
cax = fig.add_axes([(L+2*PW+GAP+CBG)/FW, B/FH, CBW/FW, PH/FH])
vmax = np.log10(max(EN.max(), ES.max()))
for ax, E, lab in zip(axes, (EN, ES), ('Neutral', r'Stable  ($T_\mathrm{s}$ = 285 K)')):
    F._buildings(ax)
    Em = np.ma.masked_where(E < 0.5*one, E)
    m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=F.MAPCMAP, vmin=np.log10(one),
                      vmax=vmax, shading='gouraud', rasterized=True, zorder=2)
    for sx, sy in F.S.values():
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=3.6, mew=.7, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(*F.MAP_XLIM); ax.set_ylim(*F.MAP_YLIM)
    ax.set_xlabel('x   (m)')
    ax.set_title(lab, fontsize=8.5, pad=4)
_wind(axes[0]); _scale(axes[1])
axes[0].set_ylabel('y   (m)')
axes[1].tick_params(labelleft=False)
cb = fig.colorbar(m, cax=cax)
cb.set_label(r'exposure below 2 m   ($\log_{10}$, particle s m$^{-3}$)', fontsize=7.5)
cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
finish(fig, axes[0], 'f39_exposure_map_neutral_vs_stable', minor=False, tight=False)


# ---------------------------------------------------------------- f40 차이 한 장
D = np.log10((ES + 0.5*one)/(EN + 0.5*one))      # 한쪽만 0 일 때 발산 방지
D = np.ma.masked_where((ES < 0.5*one) & (EN < 0.5*one), D)
fig, ax, cax = F._mapfig()
F._buildings(ax)
m = ax.pcolormesh(gx, gy, D.T, cmap='RdBu_r', vmin=-1.0, vmax=1.0,
                  shading='gouraud', rasterized=True, zorder=2)
cb = fig.colorbar(m, cax=cax)
cb.set_label(r'stable $-$ neutral   ($\log_{10}$ exposure ratio)', fontsize=7.5)
cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
for sx, sy in F.S.values():
    ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.4, mew=.8, zorder=6)
F._mapdeco(ax)
finish(fig, ax, 'f40_exposure_difference_stable_vs_neutral', minor=False, tight=False)

sel = (ES >= 0.5*one) | (EN >= 0.5*one)
print(f'셀 {sel.sum():,} 개 중 안정이 높은 곳 {100*(ES[sel] > EN[sel]).mean():.1f} %')
print(f'log10 비 중앙값 {np.median(D.compressed()):+.3f}  (= ×{10**np.median(D.compressed()):.2f})')
print(f'노출 셀 수  중립 {int((EN>=one).sum()):,}  안정 {int((ES>=one).sum()):,}')
