#!/usr/bin/env python3
"""도입부(01 MOTIVATION) 그림 두 종 — EPZ 가 도시 규모로 내려온다는 한 장.

  f50  거리 사다리 — 로그 거리축에 규제 구역과 이 연구의 대역을 함께 얹는다
  f51  2패널 지도  — 좌: 대형로 EPZ 가 덮는 범위, 우: 부지경계 EPZ 가 블록 위에 놓인 모습

서체·색·저장 규격은 발표자료의 다른 그림과 같다(figstyle 재사용).
부지경계 반경은 부지마다 다르므로 예시값(300 m)이며 그림에 그렇게 적는다.

  python3 make_motivation_figs.py
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.colors import LinearSegmentedColormap

DECK = '/home/kdh/Urban-SPH/runs/lambdap/deposition'
sys.path.insert(0, DECK)
import figstyle as st
HERE = os.path.dirname(os.path.abspath(__file__))
st.OUT = os.path.join(HERE, 'figs')
from figstyle import C, finish

BLDG = '#DCE1E6'
HF = np.load('/home/kdh/Urban-SPH/runs/lambdap/heightfield.npz')
H, X0, Y0, HC = HF['H'], float(HF['x0']), float(HF['y0']), float(HF['cell'])
CX = X0 + H.shape[1]*HC/2.0                      # 도시 블록 중심 [m]
CY = Y0 + H.shape[0]*HC/2.0
SITE_R = 0.30                                    # 부지경계 EPZ 예시 반경 [km]


def buildings(ax, ox=0.0, oy=0.0, scale=1.0):
    """건물 발자국. scale=1e-3 이면 m -> km."""
    bx = (X0 + np.arange(H.shape[1]+1)*HC - ox)*scale
    by = (Y0 + np.arange(H.shape[0]+1)*HC - oy)*scale
    ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                  cmap=LinearSegmentedColormap.from_list('b', [BLDG, BLDG]),
                  shading='flat', rasterized=True, zorder=2)


# ------------------------------------------------------------------ f50 사다리
def ladder():
    fig = plt.figure(figsize=(8.40, 2.55))
    ax = fig.add_axes([0.045, 0.25, 0.94, 0.70])
    ax.set_xscale('log'); ax.set_xlim(0.075, 48)
    ax.set_ylim(0.05, 4.45)
    ax.set_yticks([]); ax.set_xlabel('distance from the release   (km)')
    ax.tick_params(top=False, which='both')
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)

    def entry(y, x0, x1, color, label, sub, lw=5.5):
        if x1 > x0:
            ax.plot([x0, x1], [y, y], '-', c=color, lw=lw, solid_capstyle='butt', zorder=3)
        else:                                          # 단일 거리 — 짧은 세로 눈금
            ax.plot([x0, x0], [y-0.13, y+0.13], '-', c=color, lw=1.6, zorder=3)
        xc = np.sqrt(x0*x1)
        ax.text(xc, y+0.46, label, ha='center', va='bottom',
                fontsize=8.5, fontweight='bold', color=color, zorder=4)
        ax.text(xc, y+0.15, sub, ha='center', va='bottom',
                fontsize=7.2, color=C['grey'], zorder=4)

    entry(3.45, 20, 30, C['grey'], 'UPZ  20 – 30 km', 'Korea — set for large LWRs')
    entry(2.45, 16.1, 16.1, C['grey'], '10-mile EPZ', 'US — large LWR')
    entry(1.45, 3, 5, C['grey'], 'PAZ  3 – 5 km', 'Korea — set for large LWRs')
    entry(2.45, 0.10, 0.50, C['blue'], 'Site boundary EPZ',
          'SMR — 10 CFR 50.160 / NuScale')
    entry(0.45, 0.20, 1.00, C['orange'], 'This work',
          'building-resolved, street by street')
    finish(fig, ax, 'f50_epz_ladder', minor=False, tight=False)


# ------------------------------------------------------------------ f51 2패널
def maps():
    PW, PH, GAP, L, B, T = 3.30, 3.30, 0.62, 0.46, 0.42, 0.34
    FW, FH = L + PW + GAP + PW + 0.14, B + PH + T
    fig = plt.figure(figsize=(FW, FH))
    axes = [fig.add_axes([(L+k*(PW+GAP))/FW, B/FH, PW/FW, PH/FH]) for k in (0, 1)]

    # --- 좌: 대형로 EPZ
    ax = axes[0]
    for r, lab, col in ((30, 'UPZ  30 km', C['grey']), (16.1, '10 mile', C['grey']),
                        (5, 'PAZ  5 km', C['grey'])):
        ax.add_patch(Circle((0, 0), r, fill=False, ec=col, lw=1.0,
                            ls='--' if r == 16.1 else '-', zorder=3))
        ax.text(0, r, '  '+lab, ha='left', va='bottom', fontsize=7.2, color=col, zorder=4)
    buildings(ax, CX, CY, 1e-3)                                  # 실제 축척의 도시(점만하다)
    ax.add_patch(Rectangle((-0.46, -0.46), 0.92, 0.92, fill=False, ec=C['orange'],
                           lw=1.0, zorder=5))
    ax.annotate('the city, to scale\n0.9 km across', xy=(0.5, 0.5), xytext=(8.5, 9.5),
                fontsize=7.2, color=C['orange'], zorder=6,
                arrowprops=dict(arrowstyle='-', color=C['orange'], lw=.8))
    ax.set_xlim(-33, 33); ax.set_ylim(-33, 33); ax.set_aspect('equal')
    ax.set_title('Large light-water reactor', fontsize=8.5, pad=4)
    ax.set_xlabel('km'); ax.tick_params(labelsize=7)

    # --- 우: SMR 부지경계 EPZ
    ax = axes[1]
    buildings(ax, CX, CY, 1e-3)
    ax.add_patch(Circle((0, 0), SITE_R, fill=False, ec=C['blue'], lw=1.3, zorder=5))
    ax.text(0, SITE_R, '  site boundary EPZ\n  (illustrative, 300 m)', ha='left',
            va='bottom', fontsize=7.2, color=C['blue'], zorder=6)
    ax.plot([0.18, 0.38], [-0.42, -0.42], '-', c=C['ink'], lw=1.8, zorder=7)
    ax.text(0.28, -0.40, '200 m', ha='center', fontsize=7.2, zorder=7)
    ax.set_xlim(-0.50, 0.50); ax.set_ylim(-0.50, 0.50); ax.set_aspect('equal')
    ax.set_title('Small modular reactor', fontsize=8.5, pad=4)
    ax.set_xlabel('km'); ax.tick_params(labelsize=7)

    finish(fig, axes[0], 'f51_epz_maps', minor=False, tight=False)




# ------------------------------------------------------------------ f52 한 장
def city():
    """원자로가 도시 바로 옆에 앉은 평면도 — 부지경계 EPZ 가 블록에 닿는다.

    배치는 예시다(실제 부지 아님). 건물 발자국만 Leipzig 실측이다.
    """
    XL, XR, YB, YT = -520, 1180, -472, 472
    PW = 6.10                                    # 축 박스 [in]
    PH = PW*(YT-YB)/(XR-XL)
    L, B, T, R = 0.10, 0.34, 0.30, 0.10
    FW, FH = L+PW+R, B+PH+T
    fig = plt.figure(figsize=(FW, FH))
    ax = fig.add_axes([L/FW, B/FH, PW/FW, PH/FH])

    RX, RY, EPZ = -150.0, 0.0, 300.0             # 원자로 위치 [m], 부지경계 반경 [m]
    buildings(ax)                                # m 단위 그대로

    ax.add_patch(Circle((RX, RY), EPZ, fill=True, fc='#EAF0F5', ec=C['blue'],
                        lw=1.2, ls='--', zorder=1))
    ax.text(RX, RY-EPZ-18, 'site boundary EPZ   300 m', ha='center', va='top',
            fontsize=8, color=C['blue'], zorder=6)

    ax.add_patch(Rectangle((RX-85, RY-85), 170, 170, fc='white', ec=C['ink'],
                           lw=.9, zorder=4))
    ax.add_patch(Circle((RX, RY), 42, fc=C['blue'], ec=C['blue'], lw=0, zorder=5))
    ax.text(RX, RY+100, 'SMR', ha='center', va='bottom', fontsize=8.5,
            fontweight='bold', color=C['ink'], zorder=6)

    NB = (201.0, -121.0)                         # 실측 최근접 건물 — 원자로에서 372 m
    ax.annotate('', xy=NB, xytext=(RX, RY),
                arrowprops=dict(arrowstyle='<->', color=C['orange'], lw=1.0), zorder=6)
    ax.text((RX+NB[0])/2, (RY+NB[1])/2 - 75, '370 m to the nearest block', ha='center',
            va='top', fontsize=8, color=C['orange'], zorder=6)
    import csv as _csv
    src = list(_csv.DictReader(open('/home/kdh/Urban-SPH/runs/lambdap/sources13.csv')))
    ax.plot([float(r['x']) for r in src], [float(r['y']) for r in src], 'o',
            mfc=C['red'], mec='white', ms=5.0, mew=.8, ls='none', zorder=6)
    ax.text(760, 405, '13 release positions inside the blocks', ha='center', va='bottom',
            fontsize=8, color=C['red'], zorder=6)

    ax.annotate('', xy=(-180, 405), xytext=(-380, 405),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.1), zorder=6)
    ax.text(-380, 417, 'wind', fontsize=8, zorder=6)

    ax.plot([820, 1020], [-430, -430], '-', c=C['ink'], lw=1.8, zorder=7)
    ax.text(920, -422, '200 m', ha='center', fontsize=8, zorder=7)

    ax.set_xlim(XL, XR); ax.set_ylim(YB, YT); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(C['light'])
    finish(fig, ax, 'f52_smr_next_to_city', minor=False, tight=False)


if __name__ == '__main__':
    ladder(); maps(); city()
    print('저장:', st.OUT)
