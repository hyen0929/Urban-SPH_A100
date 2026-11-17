#!/usr/bin/env python3
"""논문 그림 1~9 일괄 생성 (paperstyle.py 스타일). 이전 dep_paper_figs{,2,3,4}.py 를 대체한다.

  python3 make_figures.py [그림번호 ...]      번호 생략 시 전부
"""
import csv, json, os, re, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from paperstyle import (MM, W1, W2, C, panel, save, statbox, place_labels, minor_ticks, scatter_sources)

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A_RUN = os.path.join(RUNS, '2026-10-20_leipzig-s13-most-wall')
COURT = {3, 6, 10}
R_MAIN = 150                      # 주 반경(형태 분석 사전 정의값)

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']


def csvd(path, key='src'):
    return {int(r[key]): {k: (float(v) if v not in ('', None) else np.nan) for k, v in r.items()}
            for r in csv.DictReader(open(path))}


lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
IDS = sorted(S)
KEEP = np.array([i not in COURT for i in IDS])
X = np.array([lam[R_MAIN][i] for i in IDS])
met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
corr = list(csv.DictReader(open(os.path.join(HERE, 'dep_lambda_corr.csv'))))
ped = csvd(os.path.join(HERE, 'ped_exposure.csv'))
bands = list(csv.DictReader(open(os.path.join(HERE, 'ped_exposure_bands.csv'))))
flow = csvd(os.path.join(HERE, 'flow_predictors.csv'))
XLAB = rf'$\lambda_\mathrm{{p}}$  (R = {R_MAIN} m)'


def sp(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


def fit_line(ax, x, y, k, log=False):
    yy = np.log10(y) if log else y
    b, a = np.polyfit(x[k], yy[k], 1)
    xs = np.linspace(x.min()-.02, x.max()+.02, 30)
    ax.plot(xs, 10**(a+b*xs) if log else a+b*xs, '-', c=C['ink'], lw=1.0, zorder=3)
    return sp(x[k], yy[k])


def legend_sources(ax, loc='best'):
    ax.legend(handles=[Line2D([], [], ls='', marker='o', ms=4.5, mfc=C['blue'], mec='white', mew=.6,
                              label='street-level source'),
                       Line2D([], [], ls='', marker='o', ms=4.5, mfc='white', mec=C['red'], mew=1.3,
                              label='enclosed courtyard')], loc=loc)


# ============================== Fig 1 ==============================
def fig1():
    fig, axs = plt.subplots(1, 3, figsize=(W2, 62*MM))
    spec = [('area90_m2', r'Footprint area $A_{90}$  (m$^2$)', True, 'upper right', (4e2, 3.2e4)),
            ('f_r150', r'Trapped fraction  ($r \leq 150$ m)', False, 'lower right', (0.78, 1.07)),
            ('f_dep', 'Deposited fraction', False, 'upper right', (0.09, 0.46))]
    for ax, (key, ylab, log, loc, ylim) in zip(axs, spec):
        y = np.array([met['A'][i][key] for i in IDS])
        # 축 한계를 먼저 정해야 라벨 자동배치가 축 밖으로 나가지 않는다
        if log: ax.set_yscale('log')
        ax.set_xlim(0.03, 0.55); ax.set_ylim(*ylim)
        ax.set_xlabel(XLAB); ax.set_ylabel(ylab); minor_ticks(ax)
        rs = fit_line(ax, X, y, KEEP, log)
        scatter_sources(ax, X, y, IDS)
        statbox(ax, rf'$\rho_\mathrm{{s}} = {rs:+.2f}$', loc)
    legend_sources(axs[1], loc='upper left')
    for a, l in zip(axs, 'abc'): panel(a, l)
    fig.tight_layout(w_pad=2.4)
    save(fig, 'fig1_lambda_scatter')


# ============================== Fig 2 ==============================
def fig2():
    ts = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'dep_timeseries.csv'))):
        if r['run'] == 'A':
            ts.setdefault(int(r['src']), []).append((float(r['time']), int(r['n_released']), int(r['n_dep'])))
    for k in ts: ts[k].sort()
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62*MM))
    style = [('area90_m2', 'Footprint area', C['blue'], 'o'), ('f_r150', 'Trapped fraction', C['green'], 's'),
             ('shift_m', 'Centroid shift', C['purple'], '^'), ('f_dep', 'Deposited fraction', C['orange'], 'D'),
             ('elong', 'Elongation', C['sky'], 'v')]
    for mt, lab, c, m in style:
        rr = sorted((int(float(r['R'])), float(r['spearman_nocourt']))
                    for r in corr if r['metric'] == mt and r['run'] == 'A')
        axs[0].plot([a for a, _ in rr], [b for _, b in rr], marker=m, c=c, ms=4, lw=1.3, label=lab, mec='white', mew=.5)
    axs[0].axhline(0, c=C['ink'], lw=.8)
    axs[0].axvspan(75, 100, color='#EDF1F4', zorder=0, lw=0)
    axs[0].set_xlabel('Window radius R  (m)')
    axs[0].set_ylabel(r'Spearman $\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
    axs[0].set_ylim(-1.05, 1.15); axs[0].set_xlim(0, 320); minor_ticks(axs[0])
    axs[0].legend(loc='center right', ncol=1)
    for i in IDS:
        a = np.array([(t, n/max(rel, 1)) for t, rel, n in ts[i]])
        axs[1].plot(a[:, 0], a[:, 1], lw=1.3 if i in COURT else 1.0,
                    c=C['red'] if i in COURT else C['sky'], alpha=1.0 if i in COURT else .75, zorder=4 if i in COURT else 3)
    place_labels(axs[1], [ts[i][-1][0] for i in IDS],
                 [ts[i][-1][2]/max(ts[i][-1][1], 1) for i in IDS], [str(i) for i in IDS], pad=6)
    axs[1].set_xlabel('Time  (s)'); axs[1].set_ylabel('Deposited fraction')
    axs[1].set_xlim(150, 316); axs[1].set_ylim(0, .45); minor_ticks(axs[1])
    axs[1].legend(handles=[Line2D([], [], c=C['sky'], lw=1.2, label='street-level'),
                           Line2D([], [], c=C['red'], lw=1.4, label='courtyard')], loc='upper left')
    for a, l in zip(axs, 'ab'): panel(a, l)
    fig.tight_layout(w_pad=2.6)
    save(fig, 'fig2_scale_and_time')


# ============================== Fig 3 ==============================
def fig3():
    d = np.load(os.path.join(HERE, 'grid_dep.npz'))
    GA, GB, lim = d['GA'], d['GB'], d['lim']
    ext = [lim[0], lim[1], lim[2], lim[3]]
    vmax = np.log10(max(GA.max(), GB.max())+1)
    fig, axs = plt.subplots(3, 1, figsize=(W1, 152*MM))
    for ax, G, lab in [(axs[0], GA, 'transient LES'), (axs[1], GB, 'frozen mean flow')]:
        im = ax.imshow(np.log10(G.T+1), origin='lower', extent=ext, cmap='magma', aspect='equal', vmin=0, vmax=vmax)
        cb = plt.colorbar(im, ax=ax, fraction=.024, pad=.012)
        cb.set_label(r'$\log_{10}(N + 1)$ per 100 m$^2$', fontsize=7)
        cb.ax.tick_params(labelsize=7, width=.8, length=2.6)
        statbox(ax, lab, 'upper left', 7.5)
    D = np.log10(GA.T+1) - np.log10(GB.T+1)
    im = axs[2].imshow(D, origin='lower', extent=ext, cmap='RdBu_r', aspect='equal', vmin=-1.2, vmax=1.2)
    cb = plt.colorbar(im, ax=axs[2], fraction=.024, pad=.012)
    cb.set_label('transient − mean flow', fontsize=7); cb.ax.tick_params(labelsize=7, width=.8, length=2.6)
    statbox(axs[2], 'difference', 'upper left', 7.5)
    for ax in axs:
        for i, (sx, sy) in S.items():
            ax.plot(sx, sy, 'o', mfc='none', mec='#7FD4FF' if ax is not axs[2] else C['ink'], ms=4, mew=.9)
        ax.set_xlabel('x  (m)'); ax.set_ylabel('y  (m)'); minor_ticks(ax)
    for a, l in zip(axs, 'abc'): panel(a, l, dx=-0.12, dy=1.07)
    fig.tight_layout(h_pad=1.4)
    save(fig, 'fig3_deposition_maps')


# ============================== Fig 4 ==============================
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(W2*.68, 62*MM))
    for ax, key, lab, log in [(axs[0], 'f_dep', 'Deposited fraction', False),
                              (axs[1], 'area90_m2', r'Footprint area $A_{90}$  (m$^2$)', True)]:
        a = np.array([met['A'][i][key] for i in IDS]); b = np.array([met['B'][i][key] for i in IDS])
        lo, hi = min(a.min(), b.min())*.75, max(a.max(), b.max())*1.35
        if log: ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.plot([lo, hi], [lo, hi], '-', c=C['grey'], lw=1.0, zorder=2)
        scatter_sources(ax, a, b, IDS)
        ax.set_xlabel(f'{lab}, transient'); ax.set_ylabel(f'{lab}, mean flow')
        statbox(ax, f'median ratio {np.median(b/a):.2f}', 'upper left')
        minor_ticks(ax)
    for a, l in zip(axs, 'ab'): panel(a, l)
    fig.tight_layout(w_pad=2.6)
    save(fig, 'fig4_transient_vs_mean')


# ============================== Fig 5 ==============================
def fig5():
    prof = list(csv.DictReader(open(os.path.join(HERE, 'canopy_profile.csv'))))
    z = np.array([float(r['z_m']) for r in prof]); U = np.array([float(r['U_mean']) for r in prof])
    K = np.array([float(r['k_turb']) for r in prof])
    J = json.load(open(os.path.join(HERE, 'canopy_aero.json')))
    fit, mac, mo = J['fit'], J['macdonald'], J['morphology']
    Hb = mo['H_built']
    fig, axs = plt.subplots(1, 2, figsize=(W2*.68, 66*MM))
    ax = axs[0]
    ax.axhspan(fit['fit_lo'], fit['fit_hi'], color='#EDF1F4', zorder=0, lw=0)
    ax.plot(U, z, 'o', ms=3.4, c=C['blue'], mec='white', mew=.4, zorder=4, label='LES 150-s mean')
    zz = np.linspace(fit['fit_lo'], fit['fit_hi'], 60)
    ax.plot(fit['ustar']/0.4*np.log((zz-fit['d'])/fit['z0']), zz, '-', c=C['ink'], lw=1.4, zorder=5,
            label='log law fit')
    ax.axhline(Hb, ls=(0, (4, 2)), c=C['grey'], lw=1.0, zorder=3)
    ax.set_xlabel(r'$\overline{U}$  (m s$^{-1}$)'); ax.set_ylabel('z  (m)')
    ax.set_ylim(0, 165); ax.set_xlim(0, 9.8); minor_ticks(ax)
    ax.legend(loc='lower right')
    statbox(ax, f"$d$ = {fit['d']:.1f} m\n$z_0$ = {fit['z0']:.2f} m\n$u_*$ = {fit['ustar']:.2f} m s$^{{-1}}$\n"
                f"Macdonald: {mac['d']:.1f}, {mac['z0']:.2f}", 'upper left')
    ax2 = axs[1]
    ax2.plot(K, z, 's', ms=3.2, c=C['green'], mec='white', mew=.4, zorder=4)
    ax2.axhline(Hb, ls=(0, (4, 2)), c=C['grey'], lw=1.0)
    ax2.annotate(r'$H_\mathrm{built}$ = %.1f m' % Hb, (K.max()*.96, Hb+4), fontsize=7, color=C['grey'], ha='right')
    ax2.set_xlabel(r'$\overline{k}$  (m$^2$ s$^{-2}$)'); ax2.set_ylabel('z  (m)')
    ax2.set_ylim(0, 165); minor_ticks(ax2)
    for a, l in zip(axs, 'ab'): panel(a, l)
    fig.tight_layout(w_pad=2.6)
    save(fig, 'fig5_canopy_profile')


# ============================== Fig 6 ==============================
def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(W2*.68, 62*MM))
    y = np.array([met['A'][i]['f_dep'] for i in IDS])
    for ax, xv, xlab in [(axs[0], X, XLAB),
                         (axs[1], np.array([flow[i]['k50'] for i in IDS]),
                          r'$\overline{k}$ within 50 m, z < 10 m  (m$^2$ s$^{-2}$)')]:
        ax.set_xlim(xv.min()-0.06*(xv.max()-xv.min()), xv.max()+0.10*(xv.max()-xv.min()))
        ax.set_ylim(0.09, 0.46)
        b, a = np.polyfit(xv[KEEP], y[KEEP], 1)
        xs = np.linspace(xv.min()*.95, xv.max()*1.04, 30)
        ax.plot(xs, a+b*xs, '-', c=C['ink'], lw=1.0, zorder=3)
        scatter_sources(ax, xv, y, IDS)
        ax.set_xlabel(xlab); ax.set_ylabel('Deposited fraction')
        statbox(ax, rf'$\rho_\mathrm{{s}} = {sp(xv[KEEP], y[KEEP]):+.2f}$', 'upper left')
        minor_ticks(ax)
    for a, l in zip(axs, 'ab'): panel(a, l)
    fig.tight_layout(w_pad=2.6)
    save(fig, 'fig6_flow_predictor')


# ============================== Fig 7 ==============================
def fig7():
    G = np.load(os.path.join(HERE, 'ped_grid.npz'))
    lim, cell = G['lim'], float(G['cell'])
    fig = plt.figure(figsize=(W2, 112*MM))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 1.0], hspace=.62, wspace=.68)
    show = [(8, 'open'), (13, 'semi-open'), (4, 'dense'), (10, 'courtyard')]
    vmax = max(np.log10(G[f'S{i}'].max()/(cell*cell*2.0)+1) for i, _ in show)
    for c, (i, tag) in enumerate(show):
        ax = fig.add_subplot(gs[0, c])
        E = G[f'S{i}']/(cell*cell*2.0); sx, sy = S[i]; half = 220
        xi0 = max(int((sx-half-lim[0])/cell), 0); xi1 = int((sx+half-lim[0])/cell)
        yi0 = max(int((sy-half-lim[2])/cell), 0); yi1 = int((sy+half-lim[2])/cell)
        im = ax.imshow(np.log10(E[xi0:xi1, yi0:yi1].T+1), origin='lower', cmap='magma', vmin=0, vmax=vmax,
                       extent=[sx-half, sx+half, sy-half, sy+half], aspect='equal')
        ax.plot(sx, sy, 'o', mfc='none', mec='#7FD4FF', ms=5, mew=1.0)
        ax.set_title(f'S{i} · {tag}   ' + rf'$\lambda_\mathrm{{p}}$ {lam[R_MAIN][i]:.2f}', fontsize=7.5, pad=4)
        ax.set_xticks([sx-200, sx, sx+200]); ax.set_xticklabels(['−200', '0', '200'])
        ax.set_yticks([sy-200, sy, sy+200]); ax.set_yticklabels(['−200', '0', '200'] if c == 0 else [])
        ax.set_xlabel(r'$x - x_\mathrm{src}$  (m)', fontsize=7.5)
        if c == 0: ax.set_ylabel(r'$y - y_\mathrm{src}$  (m)', fontsize=7.5)
        if c == 3:
            cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.03)
            cb.set_label(r'$\log_{10}$ exposure', fontsize=7); cb.ax.tick_params(labelsize=7, width=.8, length=2.6)
    panel(fig.axes[0], 'a', dx=-0.26, dy=1.28)
    spec = [('area_1pct_m2', r'Area  (m$^2$)', True, 'lower left'),
            ('peak_s_per_m3', r'Peak  (particle s m$^{-3}$)', False, 'upper left'),
            ('exposure_s_per_particle', r'Total  (s per particle)', False, 'lower left')]
    for c, (key, ylab, log, loc) in enumerate(spec):
        ax = fig.add_subplot(gs[1, c])
        y = np.array([ped[i][key] for i in IDS])
        if log: ax.set_yscale('log'); ax.set_ylim(y.min()*.55, y.max()*2.2)
        else: ax.set_ylim(y.min()-0.18*(y.max()-y.min()), y.max()+0.22*(y.max()-y.min()))
        ax.set_xlim(0.03, 0.55)
        ax.set_xlabel(XLAB); ax.set_ylabel(ylab); minor_ticks(ax)
        rs = fit_line(ax, X, y, KEEP, log)
        scatter_sources(ax, X, y, IDS)
        statbox(ax, rf'$\rho_\mathrm{{s}} = {rs:+.2f}$', loc)
        panel(ax, 'bcd'[c], dx=-0.34)
    ax = fig.add_subplot(gs[1, 3])
    bl = ['0-1.5', '0-2', '0-4', '2-5', '5-10', '10-20']
    for key, lab, c_, m in [('area_1pct_m2', 'area', C['blue'], 'o'), ('frac_within_150m', 'trapped', C['green'], 's'),
                            ('exposure_s_per_particle', 'total', C['orange'], 'D')]:
        v = []
        for bb in bl:
            rows = sorted([r for r in bands if r['band'] == bb], key=lambda r: int(r['src']))
            y = np.array([float(r[key]) for r in rows]); xs_ = np.array([lam[R_MAIN][int(r['src'])] for r in rows])
            k2 = np.array([int(r['src']) not in COURT for r in rows])
            v.append(sp(xs_[k2], y[k2]))
        ax.plot(range(len(bl)), v, marker=m, c=c_, ms=4, lw=1.3, label=lab, mec='white', mew=.5)
    ax.axhline(0, c=C['ink'], lw=.8)
    ax.set_xticks(range(len(bl))); ax.set_xticklabels(bl, rotation=45, fontsize=7)
    ax.set_xlabel('Height band  (m)'); ax.set_ylabel(r'$\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
    ax.set_ylim(-1.55, 1.15)
    ax.legend(loc='lower center', ncol=3, columnspacing=.6, handlelength=1.1, fontsize=6.8, borderpad=.3)
    panel(ax, 'e', dx=-0.34)
    save(fig, 'fig7_pedestrian_exposure')


# ============================== Fig 8 ==============================
def fig8():
    runs = [('2 m', '2026-09-17_leipzig-dt01-wall-meanflow', C['blue'], 'o'),
            ('3 m', '2026-09-17_leipzig-dt01-wall3m-meanflow', C['green'], 's'),
            ('5 m', '2026-09-17_leipzig-dt01-wall5m-meanflow', C['purple'], '^')]
    W = {}
    for lab, d, _, _ in runs:
        p, F = read_vtk(os.path.join(RUNS, d, 'out', 'LDM_30000stp.vtk'))
        b = F['buffer_type'].astype(int); s = F['p_type'].astype(int)-100
        W[lab] = (np.array([int(((s == i) & (b == 1)).sum()) for i in IDS]),
                  np.array([int(((s == i) & (b == 2)).sum()) for i in IDS]), p[b == 2, 2])
    fig, axs = plt.subplots(1, 3, figsize=(W2, 62*MM))
    for lab, _, col, mk in runs:
        g, w, _z = W[lab]
        ratio = np.maximum(w/np.maximum(g, 1), 1e-4)
        axs[0].plot(X[KEEP], ratio[KEEP], ls='', marker=mk, ms=4.6, mfc=col, mec='white', mew=.5,
                    label=rf'$d_\mathrm{{w}}$ = {lab}')
        axs[0].plot(X[~KEEP], ratio[~KEEP], ls='', marker=mk, ms=5.2, mfc='white', mec=C['red'], mew=1.2)
    axs[0].set_yscale('log'); axs[0].set_xlabel(XLAB); axs[0].set_ylabel('wall + roof / ground deposition')
    axs[0].set_xlim(0.03, 0.55); axs[0].legend(loc='lower right')
    statbox(axs[0], r'$\rho_\mathrm{s}$ = +0.84 to +0.89', 'lower left')
    g2, w2, z2 = W['2 m']
    axs[1].hist(z2, bins=np.arange(0, 60, 2.5), color=C['blue'], edgecolor='white', lw=.4)
    axs[1].axvline(18.8, ls=(0, (4, 2)), c=C['grey'], lw=1.0)
    axs[1].annotate(r'$H_\mathrm{built}$', (17.6, axs[1].get_ylim()[1]*.55), fontsize=7.5,
                    color=C['grey'], ha='right')
    axs[1].set_xlabel('Deposit height z  (m)'); axs[1].set_ylabel('Wall + roof deposits')
    statbox(axs[1], '63 % above 10 m\n(roofs)', 'upper right')
    minor_ticks(axs[1])
    N = 99988
    vals = [sp(X[KEEP], (W['2 m'][0]/N)[KEEP])] + [sp(X[KEEP], ((W[l][0]+W[l][1])/N)[KEEP]) for l, _, _, _ in runs]
    labs = ['ground\nonly', '+ wall\n2 m', '+ wall\n3 m', '+ wall\n5 m']
    axs[2].bar(range(4), vals, color=[C['orange']]+[C['blue']]*3, width=.62, edgecolor='white', lw=.6)
    axs[2].axhline(0, c=C['ink'], lw=.8)
    axs[2].set_xticks(range(4)); axs[2].set_xticklabels(labs, fontsize=7)
    axs[2].set_ylabel(r'$\rho_\mathrm{s}$ (deposited fraction vs $\lambda_\mathrm{p}$)')
    axs[2].set_ylim(-0.62, 0.08)
    for i_, v in enumerate(vals):
        axs[2].annotate(f'{v:+.2f}', (i_, v), xytext=(0, 5), textcoords='offset points',
                        ha='center', va='top', fontsize=7.5, color='white', fontweight='bold')
    for a, l in zip(axs, 'abc'): panel(a, l)
    fig.tight_layout(w_pad=2.6)
    save(fig, 'fig8_wall_deposition')


# ============================== Fig 9 ==============================
def fig9():
    G = np.load(os.path.join(HERE, 'ped_grid.npz'))
    lim, cell = G['lim'], float(G['cell'])
    E = sum(G[f'S{i}'] for i in S)/(cell*cell*2.0)
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    B = (H > 1.0).astype(np.float64)
    gx = lim[0]+(np.arange(E.shape[0])+.5)*cell; gy = lim[2]+(np.arange(E.shape[1])+.5)*cell
    GX, GY = np.meshgrid(gx, gy, indexing='ij')
    ix = np.floor((GX-x0)/hc).astype(int); iy = np.floor((GY-y0)/hc).astype(int)
    ok = (ix >= 0) & (ix < H.shape[1]) & (iy >= 0) & (iy < H.shape[0])
    bld = np.zeros(E.shape, bool); bld[ok] = H[iy[ok], ix[ok]] > 1.0
    dmin = np.full(E.shape, 1e9)
    for i, (sx, sy) in S.items():
        dmin = np.minimum(dmin, np.hypot(GX-sx, GY-sy))

    def lamfield(R):
        rad = int(R/hc); yy, xx = np.mgrid[-rad:rad+1, -rad:rad+1]
        K = ((xx*xx+yy*yy)*hc*hc <= R*R).astype(np.float64); K /= K.sum()
        Kp = np.zeros_like(B); Kp[:K.shape[0], :K.shape[1]] = K
        Kp = np.roll(np.roll(Kp, -(K.shape[0]//2), 0), -(K.shape[1]//2), 1)
        L = np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp), B.shape)
        out = np.full(E.shape, np.nan); out[ok] = L[iy[ok], ix[ok]]; return out

    LAM = lamfield(R_MAIN)
    fig = plt.figure(figsize=(W2, 88*MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0], wspace=.52)
    ax = fig.add_subplot(gs[0, 0])
    grey = LinearSegmentedColormap.from_list('b', [(1, 1, 1, 0), (0.42, 0.45, 0.49, 1)])
    ext = [lim[0], lim[1], lim[2], lim[3]]
    ax.imshow(bld.T, origin='lower', extent=ext, cmap=grey, vmin=0, vmax=1, aspect='equal', interpolation='nearest')
    Em = np.ma.masked_where(E <= 0, np.maximum(E, 1e-3))
    im = ax.imshow(np.log10(Em.T), origin='lower', extent=ext, cmap='inferno', aspect='equal',
                   vmin=-1, vmax=np.log10(E.max()), interpolation='nearest')
    cb = plt.colorbar(im, ax=ax, fraction=.027, pad=.012)
    cb.set_label(r'$\log_{10}$ exposure  (particle s m$^{-3}$)', fontsize=7.5)
    cb.ax.tick_params(labelsize=7, width=.8, length=2.6)
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='none', mec='#7FD4FF', ms=5, mew=1.0, zorder=6)
    place_labels(ax, [S[i][0] for i in IDS], [S[i][1] for i in IDS], [str(i) for i in IDS],
                 fontsize=7, color='#BFE6FF')
    ax.annotate('', xy=(330, 445), xytext=(215, 445),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.3))
    ax.text(215, 398, 'wind', fontsize=7.5, color=C['ink'])
    ax.plot([980, 1180], [-455, -455], '-', c=C['ink'], lw=2.2)
    ax.text(1080, -447, '200 m', fontsize=7.5, ha='center', color=C['ink'])
    ax.set_xlim(180, 1180); ax.set_ylim(-480, 480)
    ax.set_xlabel('x  (m)'); ax.set_ylabel('y  (m)'); minor_ticks(ax)
    panel(ax, 'a', dx=-0.085, dy=1.05)

    ax2 = fig.add_subplot(gs[0, 1])
    edges = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.60])
    mids = (edges[:-1]+edges[1:])/2
    for (lo, hi), col, mk in [((25, 75), C['orange'], 'o'), ((50, 150), C['blue'], 's'),
                              ((100, 200), C['green'], '^'), ((150, 300), C['purple'], 'D')]:
        sel = (dmin > lo) & (dmin <= hi) & (~bld) & (E > 0) & np.isfinite(LAM)
        lv, ev = LAM[sel], E[sel]
        med = [np.median(ev[(lv >= a) & (lv < b)]) if ((lv >= a) & (lv < b)).sum() >= 20 else np.nan
               for a, b in zip(edges[:-1], edges[1:])]
        rs = sp(lv, np.log10(ev))
        ax2.plot(mids, med, marker=mk, c=col, ms=4.4, lw=1.4, mec='white', mew=.5,
                 label=rf'{lo}–{hi} m   $\rho_\mathrm{{s}}$ {rs:+.2f}')
    ax2.set_yscale('log')
    ax2.set_xlabel(rf'local $\lambda_\mathrm{{p}}$  (R = {R_MAIN} m)')
    ax2.set_ylabel(r'exposure  (particle s m$^{-3}$)')
    ax2.set_xlim(0.05, 0.56); ax2.set_ylim(0.018, 4.0)
    ax2.legend(title='distance from source', title_fontsize=7.5, loc='lower left',
               ncol=2, columnspacing=.7, handlelength=1.4, fontsize=6.8, borderpad=.35)
    minor_ticks(ax2, y=False)
    panel(ax2, 'b', dx=-0.22, dy=1.05)
    save(fig, 'fig9_domain_exposure')


if __name__ == '__main__':
    todo = sys.argv[1:] or list('123456789')
    for n in todo:
        globals()[f'fig{n}']()
        print(f'fig{n} 완료', flush=True)
