#!/usr/bin/env python3
"""그림 9 — 도메인 전체 보행자 노출 지도와, 거리 통제 셀 단위 λp 관계.

(a) 13 소스 합산 노출장 (z<=2 m, t=150~300 s 적분) + 건물 평면 밑그림
(b) 소스에서 50~150 m 떨어진 셀만 골라 국소 λp(R=150 m) 로 계층화한 노출 분포
    — 소스 13점이 아니라 수천 개 셀로 같은 관계를 보는 것. 거리 효과를 밴드로 통제한다.
"""
import csv, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
exec(open(os.path.join(HERE, 'dep_paper_figs.py')).read().split('# ---------------- 데이터')[0], globals())

G = np.load(os.path.join(HERE, 'ped_grid.npz'))
lim, cell = G['lim'], float(G['cell'])
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
TOT = sum(G[f'S{i}'] for i in S)                      # 입자·s
E = TOT/(cell*cell*2.0)                               # 입자·s/m³
nx, ny = E.shape
gx = lim[0] + (np.arange(nx)+.5)*cell
gy = lim[2] + (np.arange(ny)+.5)*cell

# --- λp(R=150) 장: 높이장에 원형 커널 FFT 합성곱 ---
hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
B = (H > 1.0).astype(np.float64)
R = 150.0
rad = int(R/hc)
yy, xx = np.mgrid[-rad:rad+1, -rad:rad+1]
K = ((xx*xx + yy*yy) * hc*hc <= R*R).astype(np.float64)
K /= K.sum()
pad = np.zeros_like(B)
Fb = np.fft.rfft2(B)
Kp = np.zeros_like(B)
Kp[:K.shape[0], :K.shape[1]] = K
Kp = np.roll(np.roll(Kp, -(K.shape[0]//2), 0), -(K.shape[1]//2), 1)
LAM = np.fft.irfft2(Fb*np.fft.rfft2(Kp), B.shape)     # 셀별 λp

# 노출 격자(10 m)에 λp 와 건물 비율을 옮긴다
GX, GY = np.meshgrid(gx, gy, indexing='ij')
ix = np.floor((GX-x0)/hc).astype(int); iy = np.floor((GY-y0)/hc).astype(int)
ok = (ix >= 0) & (ix < H.shape[1]) & (iy >= 0) & (iy < H.shape[0])
lam_g = np.full(E.shape, np.nan); bld_g = np.zeros(E.shape, bool)
lam_g[ok] = LAM[iy[ok], ix[ok]]
bld_g[ok] = H[iy[ok], ix[ok]] > 1.0
dmin = np.full(E.shape, 1e9)
for i, (sx, sy) in S.items():
    dmin = np.minimum(dmin, np.hypot(GX-sx, GY-sy))

# ================= Fig 9 =================
fig = plt.figure(figsize=(W2, 82*MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.85, 1.0], wspace=.40)

ax = fig.add_subplot(gs[0, 0])
grey = LinearSegmentedColormap.from_list('b', [(1, 1, 1, 0), (0.45, 0.47, 0.5, 1)])
ext = [lim[0], lim[1], lim[2], lim[3]]
ax.imshow(bld_g.T, origin='lower', extent=ext, cmap=grey, vmin=0, vmax=1, aspect='equal', interpolation='nearest')
Em = np.ma.masked_where(E <= 0, np.maximum(E, 1e-3))
im = ax.imshow(np.log10(Em.T), origin='lower', extent=ext, cmap='inferno', aspect='equal',
               vmin=-1.0, vmax=np.log10(E.max()), interpolation='nearest')
cb = plt.colorbar(im, ax=ax, fraction=.026, pad=.012)
cb.set_label(r'$\log_{10}$ exposure (particle s m$^{-3}$)', fontsize=6.5)
cb.ax.tick_params(labelsize=6, width=.5, length=2)
for i, (sx, sy) in S.items():
    ax.plot(sx, sy, 'o', mfc='none', mec='#7FD4FF', ms=4.5, mew=.9)
    ax.annotate(f'{i}', (sx, sy), color='#7FD4FF', fontsize=5.5, xytext=(4, 3), textcoords='offset points')
ax.set_xlim(180, 1180); ax.set_ylim(-480, 480)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)'); ax.tick_params(labelsize=6.5)
panel(ax, 'a', dx=-0.09, dy=1.05)

ax2 = fig.add_subplot(gs[0, 1])
sel = (dmin > 50) & (dmin <= 150) & (~bld_g) & np.isfinite(lam_g) & (E > 0)
lv = lam_g[sel]; ev = E[sel]
edges = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.60])
mid, med, q1, q3, n = [], [], [], [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (lv >= a) & (lv < b)
    if m.sum() < 20: continue
    mid.append((a+b)/2); med.append(np.median(ev[m]))
    q1.append(np.percentile(ev[m], 25)); q3.append(np.percentile(ev[m], 75)); n.append(int(m.sum()))
mid, med, q1, q3 = map(np.array, (mid, med, q1, q3))
ax2.fill_between(mid, q1, q3, color=OI['blue'], alpha=.18, lw=0)
ax2.plot(mid, med, 'o-', c=OI['blue'], ms=3.5, lw=1.1)
for x_, y_, nn in zip(mid, q3, n):
    ax2.annotate(f'n={nn:,}', (x_, y_), fontsize=5.5, color='#555', xytext=(0, 3),
                 textcoords='offset points', ha='center')
ax2.set_yscale('log')
ax2.set_xlabel(r'local $\lambda_\mathrm{p}$ (R = 150 m)')
ax2.set_ylabel(r'exposure (particle s m$^{-3}$)')
rs = np.corrcoef(np.argsort(np.argsort(lv)), np.argsort(np.argsort(np.log10(ev))))[0, 1]
ax2.text(.04, .06, f'cells 50–150 m from a source\n' + rf'$\rho_\mathrm{{s}} = {rs:+.2f}$ (n = {sel.sum():,})',
         transform=ax2.transAxes, fontsize=6.5, va='bottom')
panel(ax2, 'b', dx=-0.26, dy=1.05)
save(fig, 'fig9_domain_exposure')
print(f'paper/fig9 생성 · 셀 {int(sel.sum()):,}개, ρs={rs:+.2f}')
print('λp 구간별 중앙 노출:', {f'{a:.2f}-{b:.2f}': round(float(m), 1) for (a, b), m in zip(zip(edges[:-1], edges[1:]), med)})
