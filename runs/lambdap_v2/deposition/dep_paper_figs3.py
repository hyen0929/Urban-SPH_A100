#!/usr/bin/env python3
"""논문 그림 7·8 — 보행자 노출(지도·상관)과 벽면 침적(민감도)."""
import csv, json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
exec(open(os.path.join(HERE, 'dep_paper_figs.py')).read().split('# ---------------- 데이터')[0], globals())
COURT = {3, 6, 10}

lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
ped = {int(r['src']): {k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(HERE, 'ped_exposure.csv')))}
bands = list(csv.DictReader(open(os.path.join(HERE, 'ped_exposure_bands.csv'))))
G = np.load(os.path.join(HERE, 'ped_grid.npz'))
lim, cell = G['lim'], float(G['cell'])
ids = sorted(S)


def sp(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


# ================= Fig 7: 보행자 노출 =================
fig = plt.figure(figsize=(W2, 95*MM))
gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1.0], hspace=.45, wspace=.38)
show = [(8, 'open'), (13, 'semi-open'), (4, 'dense'), (10, 'courtyard')]
vmax = max(np.log10(G[f'S{i}'].max()/(cell*cell*2.0) + 1) for i, _ in show)
for c, (i, tag) in enumerate(show):
    ax = fig.add_subplot(gs[0, c])
    E = G[f'S{i}']/(cell*cell*2.0)
    sx, sy = S[i]
    half = 220
    x0i = int((sx-half-lim[0])/cell); x1i = int((sx+half-lim[0])/cell)
    y0i = int((sy-half-lim[2])/cell); y1i = int((sy+half-lim[2])/cell)
    sub = E[max(x0i, 0):x1i, max(y0i, 0):y1i]
    ax.imshow(np.log10(sub.T + 1), origin='lower', cmap='magma', vmin=0, vmax=vmax,
              extent=[sx-half, sx+half, sy-half, sy+half], aspect='equal')
    ax.plot(sx, sy, 'o', mfc='none', mec='#7FD4FF', ms=4, mew=.8)
    ax.set_title(f'S{i} · {tag}\n' + rf'$\lambda_\mathrm{{p}}$ = {lam[150][i]:.2f}', fontsize=6.5, pad=3)
    ax.set_xticks([sx-200, sx, sx+200]); ax.set_xticklabels(['−200', '0', '200'])
    ax.set_yticks([sy-200, sy, sy+200]); ax.set_yticklabels(['−200', '0', '200'] if c == 0 else [])
    ax.tick_params(labelsize=6)
    ax.set_xlabel('x − x$_\\mathrm{src}$ (m)', fontsize=6.5)
    if c == 0: ax.set_ylabel('y − y$_\\mathrm{src}$ (m)', fontsize=6.5)
panel(fig.axes[0], 'a', dx=-0.30, dy=1.30)

x = np.array([lam[150][i] for i in ids]); keep = np.array([i not in COURT for i in ids])
for c, (key, ylab, lg) in enumerate([('area_1pct_m2', r'Exposure area (m$^2$)', True),
                                     ('peak_s_per_m3', r'Peak exposure (particle s m$^{-3}$)', False),
                                     ('exposure_s_per_particle', 'Exposure per released particle (s)', False)]):
    ax = fig.add_subplot(gs[1, c])
    y = np.array([ped[i][key] for i in ids])
    ax.scatter(x[keep], y[keep], s=14, c=OI['blue'], lw=0, zorder=3)
    ax.scatter(x[~keep], y[~keep], s=16, facecolors='none', edgecolors=OI['red'], lw=.9, zorder=3)
    yy = np.log10(y) if lg else y
    b_, a_ = np.polyfit(x[keep], yy[keep], 1)
    xs = np.linspace(x.min()*.95, x.max()*1.03, 20)
    ax.plot(xs, 10**(a_+b_*xs) if lg else a_+b_*xs, '-', c='k', lw=.7)
    if lg: ax.set_yscale('log')
    for i, xi, yi in zip(ids, x, y):
        ax.annotate(f'{i}', (xi, yi), fontsize=5.5, xytext=(3, 2.5), textcoords='offset points', color='#444')
    pos = {'area_1pct_m2': (.96, .96, 'right', 'top'), 'peak_s_per_m3': (.04, .96, 'left', 'top'),
           'exposure_s_per_particle': (.04, .06, 'left', 'bottom')}[key]
    ax.text(pos[0], pos[1], rf'$\rho_\mathrm{{s}} = {sp(x[keep], yy[keep]):+.2f}$', transform=ax.transAxes,
            ha=pos[2], va=pos[3], fontsize=6.5)
    ax.set_xlabel(r'$\lambda_\mathrm{p}$ (R = 150 m)'); ax.set_ylabel(ylab)
    panel(ax, 'bcd'[c])
# 고도대 민감도
ax = fig.add_subplot(gs[1, 3])
bl = ['0-1.5', '0-2', '0-4', '2-5', '5-10', '10-20']
for key, lab, c_ in [('area_1pct_m2', 'area', OI['blue']), ('frac_within_150m', 'trapped', OI['green']),
                     ('exposure_s_per_particle', 'total', OI['orange'])]:
    v = []
    for bb in bl:
        rows = [r for r in bands if r['band'] == bb]
        y = np.array([float(r[key]) for r in sorted(rows, key=lambda r: int(r['src']))])
        xs_ = np.array([lam[150][int(r['src'])] for r in sorted(rows, key=lambda r: int(r['src']))])
        k2 = np.array([int(r['src']) not in COURT for r in sorted(rows, key=lambda r: int(r['src']))])
        v.append(sp(xs_[k2], y[k2]))
    ax.plot(range(len(bl)), v, 'o-', ms=2.6, lw=.9, c=c_, label=lab)
ax.axhline(0, c='k', lw=.5); ax.set_xticks(range(len(bl))); ax.set_xticklabels(bl, rotation=45, fontsize=6)
ax.set_xlabel('Height band (m)'); ax.set_ylabel(r'$\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
ax.set_ylim(-1.15, 1.25); ax.legend(fontsize=6, loc='lower left', ncol=3, columnspacing=.8, handlelength=1.2)
panel(ax, 'e', dx=-0.34)
save(fig, 'fig7_pedestrian_exposure')

# ================= Fig 8: 벽면 침적 =================
runs = [('2 m', '2026-09-17_leipzig-dt01-wall-meanflow'),
        ('3 m', '2026-09-17_leipzig-dt01-wall3m-meanflow'),
        ('5 m', '2026-09-17_leipzig-dt01-wall5m-meanflow')]
ns2 = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns2)
read = ns2['read_vtk']
W = {}
for lab, d in runs:
    p, F = read(os.path.join(RUNS, d, 'out', 'LDM_30000stp.vtk'))
    b = F['buffer_type'].astype(int); s = F['p_type'].astype(int) - 100
    W[lab] = (np.array([int(((s == i) & (b == 1)).sum()) for i in ids]),
              np.array([int(((s == i) & (b == 2)).sum()) for i in ids]),
              p[b == 2, 2])
fig, axs = plt.subplots(1, 3, figsize=(W2, 55*MM))
mk = dict((l, m) for l, m in zip([r[0] for r in runs], ['o', 's', '^']))
for lab, _ in runs:
    g, w, _z = W[lab]
    ratio = w/np.maximum(g, 1)
    axs[0].scatter(x[keep], ratio[keep], s=14, marker=mk[lab], c=OI['blue'] if lab == '2 m' else
                   (OI['green'] if lab == '3 m' else OI['purple']), lw=0, label=f'$d_\\mathrm{{w}}$ = {lab}')
    axs[0].scatter(x[~keep], ratio[~keep], s=16, marker=mk[lab], facecolors='none',
                   edgecolors=OI['red'], lw=.8)
axs[0].set_xlabel(r'$\lambda_\mathrm{p}$ (R = 150 m)'); axs[0].set_ylabel('wall+roof / ground deposition')
axs[0].legend(loc='upper left'); axs[0].set_yscale('log')
g2, w2, z2 = W['2 m']
axs[1].hist(z2, bins=np.arange(0, 60, 2.5), color=OI['blue'], alpha=.85)
axs[1].axvline(18.8, ls=':', c='k', lw=.7)
axs[1].text(19.8, axs[1].get_ylim()[1]*.85, r'$H_\mathrm{built}$', fontsize=6, color='#444')
axs[1].set_xlabel('Deposit height z (m)'); axs[1].set_ylabel('Wall + roof deposits')
lab_corr = ['ground\nonly', '+wall\n2 m', '+wall\n3 m', '+wall\n5 m']
N = 99988
vals = [sp(x[keep], (W['2 m'][0]/N)[keep])] + [sp(x[keep], ((W[l][0]+W[l][1])/N)[keep]) for l, _ in runs]
axs[2].bar(range(4), vals, color=[OI['orange']] + [OI['blue']]*3, width=.6)
axs[2].axhline(0, c='k', lw=.5)
axs[2].set_xticks(range(4)); axs[2].set_xticklabels(lab_corr, fontsize=6)
axs[2].set_ylabel(r'$\rho_\mathrm{s}$ (deposited fraction vs $\lambda_\mathrm{p}$)')
axs[2].set_ylim(-0.6, 0.05)
for i_, v in enumerate(vals):
    axs[2].text(i_, v-0.03, f'{v:+.2f}', ha='center', fontsize=6)
for a, l in zip(axs, 'abc'): panel(a, l)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig8_wall_deposition')
print('paper/fig7, fig8 생성')
