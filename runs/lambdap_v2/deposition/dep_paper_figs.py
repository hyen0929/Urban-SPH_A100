#!/usr/bin/env python3
"""논문 투고용 그림 (PDF 벡터 + 600 dpi PNG).

규격: 단일단 89 mm / 양단 183 mm, 본문 7 pt, 패널문자 a·b·c (굵게, 좌상단),
색은 Okabe–Ito 색맹안전 팔레트, 제목 없음(캡션은 captions.md), 눈금 바깥쪽.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
OUT = os.path.join(HERE, 'paper')
os.makedirs(OUT, exist_ok=True)
COURT = {3, 6, 10}
MM = 1/25.4
W1, W2 = 89*MM, 183*MM
OI = dict(blue='#0072B2', orange='#E69F00', green='#009E73', red='#D55E00',
          purple='#CC79A7', sky='#56B4E9', yellow='#F0E442', black='#000000')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial'],
    'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5,
    'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.01,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})


def panel(ax, letter, dx=-0.22, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='left')


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.pdf'))
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=600)
    plt.close(fig)


# ---------------- 데이터 ----------------
met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
corr = [r for r in csv.DictReader(open(os.path.join(HERE, 'dep_lambda_corr.csv')))]
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
ts = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_timeseries.csv'))):
    ts.setdefault((r['run'], int(r['src'])), []).append((float(r['time']), int(r['n_released']), int(r['n_dep'])))
for k in ts: ts[k].sort()
ids = sorted(S)
d = np.load(os.path.join(HERE, 'grid_dep.npz'))
GA, GB, lim = d['GA'], d['GB'], d['lim']


def sc(ax, x, y, log=False):
    """중정은 빈 기호, 나머지는 채운 기호. 중정 제외 회귀선."""
    keep = np.array([i not in COURT for i in ids])
    ax.scatter(x[keep], y[keep], s=14, c=OI['blue'], zorder=3, lw=0)
    ax.scatter(x[~keep], y[~keep], s=16, facecolors='none', edgecolors=OI['red'], lw=.9, zorder=3)
    yy = np.log10(y) if log else y
    b, a = np.polyfit(x[keep], yy[keep], 1)
    xs = np.linspace(x.min()*.95, x.max()*1.03, 20)
    ax.plot(xs, 10**(a + b*xs) if log else a + b*xs, '-', c=OI['black'], lw=.7, zorder=2)
    rs = np.corrcoef(np.argsort(np.argsort(x[keep])), np.argsort(np.argsort(yy[keep])))[0, 1]
    for i, xi, yi in zip(ids, x, y):
        ax.annotate(f'{i}', (xi, yi), textcoords='offset points', xytext=(3, 2.5), fontsize=5.5, color='#444')
    if log: ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda_\mathrm{p}$ (R = 75 m)')
    return rs


# ---------------- Fig 1: λp 대 침적 형상 ----------------
fig, axs = plt.subplots(1, 3, figsize=(W2, 55*MM))
x = np.array([lam[75][i] for i in ids])
for ax, key, ylab, lg in [(axs[0], 'area90_m2', r'Footprint area $A_{90}$ (m$^2$)', True),
                          (axs[1], 'f_r150', r'Trapped fraction ($r \leq 150$ m)', False),
                          (axs[2], 'f_dep', 'Deposited fraction', False)]:
    y = np.array([met['A'][i][key] for i in ids])
    rs = sc(ax, x, y, lg)
    ax.set_ylabel(ylab)
    ax.text(.97, .95 if key != 'f_r150' else .06, rf'$\rho_\mathrm{{s}} = {rs:+.2f}$',
            transform=ax.transAxes, ha='right', va='top' if key != 'f_r150' else 'bottom', fontsize=6.5)
for a, l in zip(axs, 'abc'): panel(a, l)
axs[1].legend(handles=[Line2D([], [], ls='', marker='o', ms=3.5, mfc=OI['blue'], mec='none', label='street-level source'),
                       Line2D([], [], ls='', marker='o', ms=3.5, mfc='none', mec=OI['red'], label='enclosed courtyard')],
              loc='upper left', handletextpad=.4)
fig.tight_layout(w_pad=1.8)
save(fig, 'fig1_lambda_scatter')

# ---------------- Fig 2: 반경 주사 + 시계열 ----------------
fig, axs = plt.subplots(1, 2, figsize=(W2, 55*MM))
style = [('area90_m2', 'Footprint area', OI['blue'], 'o-'), ('f_r150', 'Trapped fraction', OI['green'], 's-'),
         ('shift_m', 'Centroid shift', OI['purple'], '^-'), ('f_dep', 'Deposited fraction', OI['orange'], 'd-'),
         ('elong', 'Elongation', OI['sky'], 'v-')]
for mt, lab, c, m in style:
    rr = sorted((int(float(r['R'])), float(r['spearman_nocourt'])) for r in corr if r['metric'] == mt and r['run'] == 'A')
    axs[0].plot([a for a, _ in rr], [b for _, b in rr], m, c=c, ms=2.6, lw=.9, label=lab)
axs[0].axhline(0, c='k', lw=.5)
axs[0].axvspan(75, 100, color='#eeeeee', zorder=0, lw=0)
axs[0].set_xlabel('Window radius R (m)')
axs[0].set_ylabel(r'Spearman $\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
axs[0].set_ylim(-1.05, 1.05); axs[0].legend(loc='center right', ncol=1)
for i in ids:
    a = np.array([(t, n/max(rel, 1)) for t, rel, n in ts[('A', i)]])
    axs[1].plot(a[:, 0], a[:, 1], lw=.8, c=OI['red'] if i in COURT else OI['blue'],
                alpha=.95 if i in COURT else .45)
    axs[1].annotate(f'{i}', (a[-1, 0], a[-1, 1]), fontsize=5.5, xytext=(2, -1.5), textcoords='offset points',
                    color=OI['red'] if i in COURT else '#555')
axs[1].set_xlabel('Time (s)'); axs[1].set_ylabel('Deposited fraction')
axs[1].set_xlim(150, 312)
for a, l in zip(axs, 'ab'): panel(a, l)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig2_scale_and_time')

# ---------------- Fig 3: 침적 지도 A / B / 차이 ----------------
fig, axs = plt.subplots(3, 1, figsize=(W1, 150*MM))
ext = [lim[0], lim[1], lim[2], lim[3]]
vmax = np.log10(max(GA.max(), GB.max()) + 1)
for ax, G, lab in [(axs[0], GA, 'transient'), (axs[1], GB, 'mean flow')]:
    im = ax.imshow(np.log10(G.T + 1), origin='lower', extent=ext, cmap='magma', aspect='equal', vmin=0, vmax=vmax)
    cb = plt.colorbar(im, ax=ax, fraction=.022, pad=.015)
    cb.set_label(r'$\log_{10}(N_\mathrm{dep}$ per 100 m$^2$ + 1)', fontsize=6)
    cb.ax.tick_params(labelsize=6, width=.5, length=2)
D = np.log10(GA.T + 1) - np.log10(GB.T + 1)
im = axs[2].imshow(D, origin='lower', extent=ext, cmap='RdBu_r', aspect='equal', vmin=-1.2, vmax=1.2)
cb = plt.colorbar(im, ax=axs[2], fraction=.022, pad=.015)
cb.set_label('transient − mean flow', fontsize=6); cb.ax.tick_params(labelsize=6, width=.5, length=2)
for ax in axs:
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='none', mec='w' if ax is not axs[2] else 'k',
                ms=3.2, mew=.7)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    ax.tick_params(labelsize=6)
for a, l in zip(axs, 'abc'): panel(a, l, dx=-0.13, dy=1.06)
fig.tight_layout(h_pad=1.2)
save(fig, 'fig3_deposition_maps')

# ---------------- Fig 4: A vs B ----------------
fig, axs = plt.subplots(1, 2, figsize=(W2*0.66, 55*MM))
for ax, key, lab, lg in [(axs[0], 'f_dep', 'Deposited fraction', False),
                         (axs[1], 'area90_m2', r'Footprint area $A_{90}$ (m$^2$)', True)]:
    a = np.array([met['A'][i][key] for i in ids]); b = np.array([met['B'][i][key] for i in ids])
    keep = np.array([i not in COURT for i in ids])
    ax.scatter(a[keep], b[keep], s=14, c=OI['blue'], lw=0, zorder=3)
    ax.scatter(a[~keep], b[~keep], s=16, facecolors='none', edgecolors=OI['red'], lw=.9, zorder=3)
    lo, hi = min(a.min(), b.min())*.8, max(a.max(), b.max())*1.2
    ax.plot([lo, hi], [lo, hi], '-', c='k', lw=.6)
    if lg: ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    for i, xi, yi in zip(ids, a, b):
        ax.annotate(f'{i}', (xi, yi), fontsize=5.5, xytext=(3, 2.5), textcoords='offset points', color='#444')
    ax.set_xlabel(f'{lab}, transient'); ax.set_ylabel(f'{lab}, mean flow')
    med = np.median(b/a)
    ax.text(.04, .95, f'median ratio {med:.2f}', transform=ax.transAxes, fontsize=6.5, va='top')
for a, l in zip(axs, 'ab'): panel(a, l)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig4_transient_vs_mean')
print('paper/fig1~fig4 (.pdf/.png) 생성')
