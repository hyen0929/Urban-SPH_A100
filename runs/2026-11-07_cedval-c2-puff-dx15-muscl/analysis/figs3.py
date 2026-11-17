#!/usr/bin/env python3
"""CEDVAL 중간 비교 그림: 1차(dx 2.5, IBM_length 1.0) vs 2차(dx 1.5, IBM_length 0.5) vs 풍동.
2차는 완결 퍼프만 (t_release + T98 <= t_last). 모델은 각 런의 맞춤 상수 k 를 곱한다.
사용: python3 -I figs_cmp.py <exp_dts_dir> <out_dir>
"""
import glob, os, re, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
from figstyle import C
import matplotlib.pyplot as plt

ROOT = '/home/kdh/Urban-SPH'
EXP, OUT = sys.argv[1], sys.argv[2]
R1 = f'{ROOT}/runs/2026-10-28_cedval-c2-puff'
R2 = f'{ROOT}/runs/2026-11-03_cedval-c2-puff-dx15'
SCR = os.path.dirname(os.path.abspath(__file__))
MM = 1 / 25.4
W1, W2 = 90 * MM, 190 * MM
plt.rcParams.update({'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 8,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                     'lines.markersize': 3.5})
SRCX = -278.55
LINES = [-234, -162, -54, 198]
EXPC, C1, C2 = C['ink'], '#9AA3AD', C['blue']


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print(name)


def plabel(ax, s):
    ax.text(0.03, 0.97, s, transform=ax.transAxes, ha='left', va='top', fontsize=8, fontweight='bold')


def exp_puffs():
    res = {}
    for f in sorted(glob.glob(EXP + '/*.dts')):
        t = open(f, 'rb').read().decode('latin1')
        g = lambda k: float(re.search(r'([-\d.]+) \[' + re.escape(k), t).group(1))
        key = (g('measurement location X - full'), g('measurement location Y - full'),
               round(g('measurement location Z - full'), 2))
        for z in t.split('ZONE')[1:]:
            a = np.array([[float(v) for v in l.split()[:2]] for l in z.split('\n')[1:] if l.strip()])
            if len(a) >= 3:
                res.setdefault(key, []).append(a[:, 1].sum() * np.median(np.diff(a[:, 0])))
    return {k: np.array(v) for k, v in res.items()}


def load(run, cmpname, partial):
    ts = np.load(f'{run}/analysis/receptor_ts.npz')
    cm = np.load(f'{run}/analysis/{cmpname}')
    t, Cc = ts['t'], ts['C']
    if partial:
        use = ts['t_release'] + float(cm['T98']) <= t[-1]
        Cc = Cc[:, use]
    return ts['rec'], Cc.sum(axis=2) * np.median(np.diff(t)) * float(cm['k'])


R3 = f'{ROOT}/runs/2026-11-07_cedval-c2-puff-dx15-muscl'
C3 = C['orange']
NP = 7
EP = exp_puffs()


def load7(run):
    ts = np.load(f'{run}/analysis/receptor_ts.npz')
    t = ts['t']
    return ts['rec'], ts['C'][:, :NP].sum(axis=2) * np.median(np.diff(t))


REC, D1 = load7(R1); _, D2 = load7(R2); _, D3 = load7(R3)
key = lambda r: (REC[r, 0], REC[r, 1], round(REC[r, 2], 2))
O = np.array([EP[key(r)].mean() if key(r) in EP else np.nan for r in range(len(REC))])
HAS = ~np.isnan(O)


def scale(D):          # 측정점 전체 기하평균으로 상수 하나 맞춤
    m = D.mean(1); z = HAS & (m > 0)
    return D * np.exp(np.mean(np.log(O[z] / m[z])))


D1, D2, D3 = scale(D1), scale(D2), scale(D3)
RUNS = [(D1, C1, '--', 'Δx 2.5 m, central'), (D2, C2, '-', 'Δx 1.5 m, central'), (D3, C3, '-', 'Δx 1.5 m, MUSCL')]
q = lambda a: (a.mean(), np.percentile(a, 25), np.percentile(a, 75))


def legend_below(fig, ax, ncol=4):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=ncol, frameon=False, fontsize=6.5)


def fig_lateral():
    fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.27))
    for i, (ax, x) in enumerate(zip(axs, LINES)):
        rs = sorted([r for r in range(len(REC)) if REC[r, 0] == x and REC[r, 2] < 3], key=lambda r: REC[r, 1])
        y = REC[rs, 1]
        for D, c, ls, lab in RUNS:
            ax.plot(y, D[rs].mean(1), ls, c=c, lw=1.0, label=lab)
        e = np.array([q(EP[key(r)]) for r in rs])
        ax.vlines(y, e[:, 1], e[:, 2], color=EXPC, lw=.6)
        ax.plot(y, e[:, 0], 'o', c=EXPC, ms=3, label='wind tunnel (puff IQR)')
        ax.set_title(f'{x - SRCX:.0f} m downwind', pad=3)
        ax.set_xlabel('y (m)'); ax.set_xlim(-145, 145); ax.set_ylim(bottom=0)
        plabel(ax, '(' + 'abcd'[i] + ')')
    axs[0].set_ylabel(r'mean dosage $D^*$, z = 2 m')
    fig.tight_layout(w_pad=.8); legend_below(fig, axs[0]); save(fig, 'cmp3_lateral')


def fig_vertical():
    fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.32), sharey=True)
    for i, (ax, x) in enumerate(zip(axs, LINES)):
        rs = sorted([r for r in range(len(REC)) if REC[r, 0] == x and REC[r, 1] == 0], key=lambda r: REC[r, 2])
        z = REC[rs, 2]
        for D, c, ls, lab in RUNS:
            ax.plot(D[rs].mean(1), z, ls, c=c, lw=1.0, label=lab)
        e = np.array([q(EP[key(r)]) for r in rs])
        ax.hlines(z, e[:, 1], e[:, 2], color=EXPC, lw=.6)
        ax.plot(e[:, 0], z, 'o', c=EXPC, ms=3, label='wind tunnel (puff IQR)')
        for h in (18, 27):
            ax.axhline(h, c='#C9CED4', lw=.5, ls=':', zorder=0)
        ax.set_title(f'{x - SRCX:.0f} m downwind', pad=3)
        ax.set_xlabel(r'mean dosage $D^*$'); ax.set_xlim(left=0); ax.set_ylim(0, 100)
        plabel(ax, '(' + 'abcd'[i] + ')')
    axs[0].set_ylabel('z (m)')
    fig.tight_layout(w_pad=.6); legend_below(fig, axs[0]); save(fig, 'cmp3_vertical')


def metrics(m, o):
    f2 = np.mean((m / o >= .5) & (m / o <= 2))
    nmse = np.mean((o - m) ** 2) / (o.mean() * m.mean())
    return f2, nmse


def fig_scatter():
    fig, axs = plt.subplots(1, 3, figsize=(W2, W2 * 0.36))
    lo, hi = 1, 5000; xx = np.array([lo, hi]); g = REC[:, 2] < 3
    for i, (ax, (D, c, _, lab)) in enumerate(zip(axs, RUNS)):
        m = D.mean(1); ok = HAS & (m > 0)
        ax.fill_between(xx, xx / 2, xx * 2, color='#E3E7EC', lw=0, zorder=0)
        ax.plot(xx, xx, '-', c='#9AA3AD', lw=.6, zorder=1)
        ax.plot(O[g & ok], m[g & ok], 'o', ms=3.2, mfc=c, mec='white', mew=.3, label='z = 2 m', zorder=3)
        ax.plot(O[~g & ok], m[~g & ok], '^', ms=3.4, mfc='white', mec=c, mew=.7, label='vertical profiles', zorder=3)
        z0 = HAS & (m <= 0)
        if z0.any():
            ax.plot(O[z0], np.full(z0.sum(), lo * 1.15), 'x', c=c, ms=3, mew=.7, label='model = 0', zorder=3)
        f2, nm = metrics(m[HAS], O[HAS])
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
        ax.set_xlabel(r'wind tunnel $\overline{D}^*$')
        ax.set_title(f'{lab}\nFAC2 = {f2:.2f}, NMSE = {nm:.2f}', pad=3)
        plabel(ax, '(' + 'abc'[i] + ')')
    axs[0].set_ylabel(r'SPH $\overline{D}^*$')
    fig.tight_layout(w_pad=1.0)
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False, fontsize=6.5)
    save(fig, 'cmp3_scatter')


def fig_turb():
    P = np.load(f'{SCR}/turb_prof.npz'); zb = P['zb']; zc = (zb[:-1] + zb[1:]) / 2
    fig, axs = plt.subplots(1, 3, figsize=(W2, W2 * 0.34), sharey=True)
    for k, c, lab in (('c', C2, 'Δx 1.5 m, central'), ('M', C3, 'Δx 1.5 m, MUSCL')):
        A, U = P[f'{k}_arr'], P[f'{k}_up']
        axs[0].plot(A[:, 0], zc, '-', c=c, lw=1.0, label=lab + ', array')
        axs[0].plot(U[:, 0], zc, ':', c=c, lw=1.0, label=lab + ', upstream')
        axs[1].plot(A[:, 1], zc, '-', c=c, lw=1.0); axs[1].plot(A[:, 3], zc, '-.', c=c, lw=.8)
        axs[1].plot(U[:, 1], zc, ':', c=c, lw=1.0)
        axs[2].plot(A[:, 4], zc, '-', c=c, lw=1.0); axs[2].plot(U[:, 4], zc, ':', c=c, lw=1.0)
    axs[1].plot([], [], '-', c='#6B7480', lw=1.0, label=r'$\sigma_{u_x}$')
    axs[1].plot([], [], '-.', c='#6B7480', lw=.8, label=r'$\sigma_{u_z}$')
    axs[1].legend(loc='lower right', frameon=False, fontsize=6)
    for ax in axs:
        for h in (18, 27):
            ax.axhline(h, c='#C9CED4', lw=.5, ls=':', zorder=0)
        ax.set_xlim(left=0)
    axs[0].set_xlabel(r'time-mean $\overline{u}_x$ (m s$^{-1}$)')
    axs[1].set_xlabel(r'resolved std (m s$^{-1}$)')
    axs[2].set_xlabel(r'resolved TKE (m$^2$ s$^{-2}$)')
    axs[0].set_ylabel('z (m)'); axs[0].set_ylim(0, 50)
    for i, ax in enumerate(axs):
        plabel(ax, '(' + 'abc'[i] + ')')
    fig.tight_layout(w_pad=.8); legend_below(fig, axs[0], ncol=2); save(fig, 'cmp3_turbulence')


def fig_summary():
    g = REC[:, 2] < 3
    def line(D, x):
        if D is None:
            return np.mean([O[r] for r in range(len(REC)) if REC[r, 0] == x and g[r]])
        return np.mean([D[r].mean() for r in range(len(REC)) if REC[r, 0] == x and g[r]])
    labs = ['wind tunnel'] + [r[3] for r in RUNS]; cols = [EXPC] + [r[1] for r in RUNS]
    Ds = [None, D1, D2, D3]
    dec1 = [line(D, -162) / line(D, 198) for D in Ds]
    dec2 = [line(D, -54) / line(D, 198) for D in Ds]
    cv = [np.median([EP[key(r)].std() / EP[key(r)].mean() for r in range(len(REC)) if HAS[r]])] + \
         [np.median([D[r].std() / D[r].mean() for r in range(len(REC)) if HAS[r] and D[r].mean() > 0]) for D in (D1, D2, D3)]
    fig, axs = plt.subplots(1, 3, figsize=(W2, W2 * 0.28))
    xs = np.arange(4)
    for ax, v, yl in zip(axs, (dec1, dec2, cv), ('ground dosage ratio\n117 m / 477 m', 'ground dosage ratio\n225 m / 477 m',
                                                 'puff-to-puff CV of dosage\n(median over receptors)')):
        ax.bar(xs, v, color=cols, width=.62)
        ax.axhline(v[0], c=EXPC, lw=.6, ls=':')
        ax.set_xticks(xs); ax.set_xticklabels(['WT', 'Δx 2.5\ncentral', 'Δx 1.5\ncentral', 'Δx 1.5\nMUSCL'], fontsize=6.5)
        ax.set_ylabel(yl); ax.set_ylim(0, max(v) * 1.18)
        for xi, vi in zip(xs, v):
            ax.text(xi, vi, f'{vi:.2f}' if vi < 1 else f'{vi:.1f}', ha='center', va='bottom', fontsize=6)
    for i, ax in enumerate(axs):
        ax.text(0.97, 0.97, '(' + 'abc'[i] + ')', transform=ax.transAxes, ha='right', va='top', fontsize=8, fontweight='bold')
    fig.tight_layout(w_pad=1.0); save(fig, 'cmp3_summary')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    fig_lateral(); fig_vertical(); fig_scatter(); fig_turb(); fig_summary()
