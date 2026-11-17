#!/usr/bin/env python3
"""dep_metrics.py 산출물로 그림 생성 (레이블은 영문 — 서버에 한글 폰트가 CJK JP 뿐)."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}
plt.rcParams.update({'font.size': 9, 'figure.dpi': 130, 'savefig.bbox': 'tight'})

met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    met.setdefault(r['run'], {})[int(r['src'])] = {k: (float(v) if v else np.nan) for k, v in r.items() if k not in ('run',)}
lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
corr = [r for r in csv.DictReader(open(os.path.join(HERE, 'dep_lambda_corr.csv')))]
ids = sorted(lam[150])


def scatter(ax, R, metric, run, ylab, logy=False):
    x = np.array([lam[R][i] for i in ids]); y = np.array([met[run][i][metric] for i in ids])
    c = np.array(['#c0392b' if i in COURT else '#2c6fbb' for i in ids])
    ax.scatter(x, y, c=c, s=42, zorder=3)
    for i, xi, yi in zip(ids, x, y):
        ax.annotate(f'S{i}', (xi, yi), textcoords='offset points', xytext=(4, 3), fontsize=7.5)
    keep = np.array([i not in COURT for i in ids])
    if keep.sum() > 2:                      # 중정 제외 회귀선
        b, a = np.polyfit(x[keep], y[keep], 1)
        xs = np.linspace(x.min(), x.max(), 20)
        ax.plot(xs, a + b*xs, '--', c='#555', lw=1, zorder=2)
        rs = np.corrcoef(np.argsort(np.argsort(x[keep])), np.argsort(np.argsort(y[keep])))[0, 1]
        ax.set_title(f'{ylab}  vs  $\\lambda_p$(R={R} m), run {run}\nSpearman (10 sources, courtyards excluded) = {rs:+.2f}', fontsize=9)
    if logy: ax.set_yscale('log')
    ax.set_xlabel('$\\lambda_p$'); ax.set_ylabel(ylab); ax.grid(alpha=.3)


# --- 그림 1: 핵심 산점도 4종 (Run A, R=150) ---
fig, axs = plt.subplots(2, 2, figsize=(9, 7))
scatter(axs[0, 0], 150, 'f_dep', 'A', 'deposited fraction')
scatter(axs[0, 1], 150, 'area90_m2', 'A', 'footprint area (90% mass) [m$^2$]', logy=True)
scatter(axs[1, 0], 150, 'f_r150', 'A', 'trapped fraction (r $\\leq$ 150 m)')
scatter(axs[1, 1], 150, 'elong', 'A', 'elongation (major/minor)')
fig.suptitle('Transient run (A), t = 300 s — red = enclosed courtyard sources (S3, S6, S10)', fontsize=10)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(HERE, "fig1_scatter_R150.png")); plt.close(fig)

# --- 그림 2: 최적 반경 탐색 (Spearman vs R) ---
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, run in zip(axs, ('A', 'B')):
    for mt, lab in [('f_dep', 'deposited fraction'), ('area90_m2', 'footprint area'),
                    ('f_r150', 'trapped fraction'), ('elong', 'elongation'), ('shift_m', 'centroid shift')]:
        rr = [(int(float(r['R'])), float(r['spearman_nocourt'])) for r in corr if r['metric'] == mt and r['run'] == run]
        rr.sort(); ax.plot([a for a, _ in rr], [b for _, b in rr], 'o-', ms=3.5, label=lab)
    ax.axhline(0, c='k', lw=.8); ax.axvspan(75, 200, color='#f0f0f0', zorder=0)
    ax.set_xlabel('window radius R [m]'); ax.set_title(f'run {run}'); ax.grid(alpha=.3)
axs[0].set_ylabel('Spearman vs $\\lambda_p$ (courtyards excluded)')
axs[1].legend(fontsize=7.5, loc='lower right')
fig.suptitle('Which spatial scale of urban form does deposition respond to?', fontsize=10)
fig.savefig(os.path.join(HERE, 'fig2_radius_scan.png')); plt.close(fig)

# --- 그림 3: 침적 지도 A / B / 차이 ---
d = np.load(os.path.join(HERE, 'grid_dep.npz'))
GA, GB, lim, cell = d['GA'], d['GB'], d['lim'], float(d['cell'])
ext = [lim[0], lim[1], lim[2], lim[3]]
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
fig, axs = plt.subplots(3, 1, figsize=(11, 11))
for ax, G, t in [(axs[0], GA, 'A  transient LES'), (axs[1], GB, 'B  frozen mean flow')]:
    im = ax.imshow(np.log10(G.T + 1), origin='lower', extent=ext, cmap='inferno', aspect='equal', vmin=0, vmax=np.log10(max(GA.max(), GB.max())+1))
    plt.colorbar(im, ax=ax, label='log$_{10}$(particles per 10 m cell + 1)', fraction=.025)
    ax.set_title(f'run {t} — deposited {int(G.sum()):,}')
D = np.log10(GA.T + 1) - np.log10(GB.T + 1)
im = axs[2].imshow(D, origin='lower', extent=ext, cmap='RdBu_r', aspect='equal', vmin=-1.5, vmax=1.5)
plt.colorbar(im, ax=axs[2], label='log$_{10}$ ratio  A / B', fraction=.025)
axs[2].set_title('difference (red = more in transient, blue = more in mean flow)')
for ax in axs:
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='none', mec='#00d0ff' if i not in COURT else '#00ff70', ms=7, mew=1.2)
        ax.annotate(f'S{i}', (sx, sy), color='w', fontsize=7, xytext=(5, 4), textcoords='offset points')
    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig3_dep_maps.png")); plt.close(fig)

# --- 그림 4: A vs B 소스별 비교 ---
ab = {r['scope']: r for r in csv.DictReader(open(os.path.join(HERE, 'ab_compare.csv')))}
fig, axs = plt.subplots(1, 3, figsize=(12, 3.5))
fa = [met['A'][i]['f_dep'] for i in ids]; fb = [met['B'][i]['f_dep'] for i in ids]
axs[0].scatter(fa, fb, c=['#c0392b' if i in COURT else '#2c6fbb' for i in ids], s=40)
lims = [0, max(fa+fb)*1.1]; axs[0].plot(lims, lims, 'k--', lw=.8)
for i, a, b in zip(ids, fa, fb): axs[0].annotate(f'S{i}', (a, b), fontsize=7, xytext=(4, 3), textcoords='offset points')
axs[0].set_xlabel('deposited fraction, A'); axs[0].set_ylabel('deposited fraction, B'); axs[0].set_title('per-source deposition'); axs[0].grid(alpha=.3)
r_ab = [float(ab[f'S{i}']['pearson']) for i in ids]; s_ab = [float(ab[f'S{i}']['spearman']) for i in ids]
axs[1].bar(np.arange(13)-0.2, r_ab, .4, label='Pearson'); axs[1].bar(np.arange(13)+0.2, s_ab, .4, label='Spearman')
axs[1].set_xticks(range(13)); axs[1].set_xticklabels([f'S{i}' for i in ids], fontsize=7)
axs[1].set_ylabel('map correlation A vs B'); axs[1].set_ylim(0, 1.05); axs[1].legend(fontsize=7); axs[1].grid(alpha=.3, axis='y')
aa = [met['A'][i]['area90_m2'] for i in ids]; bb = [met['B'][i]['area90_m2'] for i in ids]
axs[2].scatter(aa, bb, c=['#c0392b' if i in COURT else '#2c6fbb' for i in ids], s=40)
m = max(aa+bb)*1.2; axs[2].plot([100, m], [100, m], 'k--', lw=.8); axs[2].set_xscale('log'); axs[2].set_yscale('log')
for i, a, b in zip(ids, aa, bb): axs[2].annotate(f'S{i}', (a, b), fontsize=7, xytext=(4, 3), textcoords='offset points')
axs[2].set_xlabel('footprint area, A [m$^2$]'); axs[2].set_ylabel('footprint area, B [m$^2$]'); axs[2].set_title('footprint area'); axs[2].grid(alpha=.3)
fig.suptitle('Transient (A) vs frozen mean flow (B)', fontsize=10)
fig.savefig(os.path.join(HERE, 'fig4_AB_compare.png')); plt.close(fig)
print('fig1~fig4 생성')
