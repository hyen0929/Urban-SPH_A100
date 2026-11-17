#!/usr/bin/env python3
"""1편 3.3 CEDVAL 검증 그림.

  FigV1  배치: 건물(지붕 종류별), 방출점, 측정점
  FigV2  지면(z=2 m) 평균 dosage 횡분포, 측정선 4 개 — 실험 vs 모델, 퍼프 간 25–75 %
  FigV3  y=0 연직 분포, 측정선 4 개
  FigV4  퍼프 간 변동: (a–c) 대표 3 지점 dosage/평균 분포  (d) 모든 측정점 CV 실험 vs 모델
  FigV5  산점도: 평균 dosage 모델 vs 실험, 1:1 과 factor-2 띠

모델 dosage 는 상수 k (compare_dose.py, 기하평균 맞춤) 를 곱해 실험 단위(10⁶·C*·T*)로 그린다.
자료: runs/<run>/analysis/receptor_ts.npz, compare_dose.npz, 실험 .dts
사용: python3 -I validation/cedval/figs.py runs/<run> <exp_dts_dir>
"""
import csv, glob, os, re, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
from figstyle import C
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = '/home/kdh/Urban-SPH'
RUN, EXP = os.path.join(ROOT, sys.argv[1]), sys.argv[2]
OUT = f'{ROOT}/paper/figs/paper1'
MM = 1 / 25.4
W1, W2 = 90 * MM, 190 * MM
plt.rcParams.update({'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 8,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                     'lines.markersize': 3.5})
SRCX = -278.55
LINES = [-234, -162, -54, 198]
EXPC, MODC = C['ink'], C['blue']


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print(name)


def plabel(ax, s):
    ax.text(0.03, 0.97, s, transform=ax.transAxes, ha='left', va='top', fontsize=8, fontweight='bold')


# ── 자료 ───────────────────────────────────────────────────────────────
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


ts = np.load(f'{RUN}/analysis/receptor_ts.npz')
cmp_ = np.load(f'{RUN}/analysis/compare_dose.npz')
K = float(cmp_['k'])
REC = ts['rec']
DM = ts['C'].sum(axis=2) * np.median(np.diff(ts['t'])) * K     # [nr, npuff], 실험 단위
EP = exp_puffs()


def key(r):
    return (REC[r, 0], REC[r, 1], round(REC[r, 2], 2))


def stats(a):
    return a.mean(), np.percentile(a, 25), np.percentile(a, 75)


# ── V1 배치 ────────────────────────────────────────────────────────────
def fig_layout():
    bl = list(csv.DictReader(open(f'{ROOT}/input/CEDVAL_C2/buildings.csv')))
    sty = {'F': ('#D5DAE0', None, 'flat 18 m'), 'X': ('#8C96A3', None, 'flat 27 m'),
           'H': ('#D5DAE0', '////', 'pitched, ridge ⊥ wind'), 'V': ('#D5DAE0', '----', 'pitched, ridge ∥ wind')}
    fig, ax = plt.subplots(figsize=(W1, W1 * 0.95))
    for b in bl:
        x, y, t = float(b['x_cedval']), float(b['y']), b['type']
        fc, h, _ = sty[t]
        ax.add_patch(Rectangle((x - 9, y - 9), 18, 18, fc=fc, ec='#6B7480', lw=.25, hatch=h))
    for t, (fc, h, lab) in sty.items():
        ax.add_patch(Rectangle((0, 0), 0, 0, fc=fc, ec='#6B7480', lw=.25, hatch=h, label=lab))
    g = REC[:, 2] < 3
    ax.plot(REC[g, 0], REC[g, 1], 'o', ms=2.6, mfc=C['red'], mec='white', mew=.3, label='receptor, z = 2 m', zorder=5)
    v = np.unique(REC[~g, :2], axis=0)
    ax.plot(v[:, 0], v[:, 1], 's', ms=3.6, mfc='none', mec=C['red'], mew=.7, label='vertical profile', zorder=5)
    ax.plot(SRCX, 0, '*', ms=7, mfc=C['orange'], mec=C['ink'], mew=.3, label='source', zorder=6)
    for x in LINES:
        ax.text(x, 314, f'{x - SRCX:.0f} m', ha='center', va='bottom', fontsize=6)
    ax.annotate('', xy=(-215, -300), xytext=(-300, -300), arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=.8))
    ax.text(-300, -292, 'wind', fontsize=6.5, va='bottom')
    ax.set_xlim(-310, 285); ax.set_ylim(-310, 310); ax.set_aspect('equal')
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, handlelength=1.4)
    save(fig, 'FigV1_cedval_layout')


# ── V2 횡분포 ──────────────────────────────────────────────────────────
def fig_lateral():
    fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.27), sharey=False)
    for i, (ax, x) in enumerate(zip(axs, LINES)):
        rs = [r for r in range(len(REC)) if REC[r, 0] == x and REC[r, 2] < 3]
        rs.sort(key=lambda r: REC[r, 1])
        y = REC[rs, 1]
        e = np.array([stats(EP[key(r)]) for r in rs]); m = np.array([stats(DM[r]) for r in rs])
        ax.fill_between(y, m[:, 1], m[:, 2], color=MODC, alpha=.18, lw=0)
        ax.plot(y, m[:, 0], '-', c=MODC, lw=1.2, label='SPH (puff IQR shaded)')
        ax.vlines(y, e[:, 1], e[:, 2], color=EXPC, lw=.6)
        ax.plot(y, e[:, 0], 'o', c=EXPC, ms=3, label='wind tunnel (puff IQR)')
        ax.set_title(f'{x - SRCX:.0f} m downwind', pad=3)
        ax.set_xlabel('y (m)'); ax.set_xlim(-145, 145); ax.set_ylim(bottom=0)
        plabel(ax, '(' + 'abcd'[i] + ')')
    axs[0].set_ylabel(r'mean dosage $D^*$')
    fig.tight_layout(w_pad=.8)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False, fontsize=6.5)
    save(fig, 'FigV2_cedval_lateral')


# ── V3 연직 ────────────────────────────────────────────────────────────
def fig_vertical():
    fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.32), sharey=True)
    for i, (ax, x) in enumerate(zip(axs, LINES)):
        rs = [r for r in range(len(REC)) if REC[r, 0] == x and REC[r, 1] == 0]
        rs.sort(key=lambda r: REC[r, 2])
        z = REC[rs, 2]
        e = np.array([stats(EP[key(r)]) for r in rs]); m = np.array([stats(DM[r]) for r in rs])
        ax.fill_betweenx(z, m[:, 1], m[:, 2], color=MODC, alpha=.18, lw=0)
        ax.plot(m[:, 0], z, '-', c=MODC, lw=1.2, label='SPH')
        ax.hlines(z, e[:, 1], e[:, 2], color=EXPC, lw=.6)
        ax.plot(e[:, 0], z, 'o', c=EXPC, ms=3, label='wind tunnel')
        for h in (18, 27):
            ax.axhline(h, c='#9AA3AD', lw=.5, ls=':', zorder=0)
        ax.set_title(f'{x - SRCX:.0f} m downwind', pad=3)
        ax.set_xlabel(r'mean dosage $D^*$'); ax.set_xlim(left=0); ax.set_ylim(0, 100)
        plabel(ax, '(' + 'abcd'[i] + ')')
    axs[0].set_ylabel('z (m)')
    axs[0].legend(loc='upper right', frameon=False, fontsize=6)
    fig.tight_layout(w_pad=.6)
    save(fig, 'FigV3_cedval_vertical')


# ── V4 퍼프 간 변동 ────────────────────────────────────────────────────
def fig_variability():
    pts = [(-234, -36, 2.02), (-54, 0, 2.02), (198, 0, 2.02)]
    fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.26), gridspec_kw={'width_ratios': [1, 1, 1, 1.05]})
    bins = np.linspace(0, 4, 21)
    for i, (ax, p) in enumerate(zip(axs[:3], pts)):
        r = [j for j in range(len(REC)) if tuple(REC[j]) == p][0]
        e, m = EP[key(r)], DM[r]
        ax.hist(e / e.mean(), bins=bins, density=True, histtype='step', color=EXPC, lw=1.0,
                label='wind tunnel (≈340 puffs)')
        ax.hist(m / m.mean(), bins=bins, density=True, color=MODC, alpha=.35, lw=0,
                label=f'SPH ({len(m)} puffs)')
        ax.set_title(f'{p[0] - SRCX:.0f} m, y = {p[1]:.0f} m', pad=3)
        ax.set_xlabel(r'$D / \overline{D}$'); ax.set_xlim(0, 4)
        plabel(ax, '(' + 'abc'[i] + ')')
    axs[0].set_ylabel('probability density')

    ax = axs[3]
    core = cmp_['exp'] > 50                    # 평균 D* > 50: 플룸 안쪽
    ax.plot(cmp_['cv_exp'][core], cmp_['cv_model'][core], 'o', ms=3, mfc=MODC, mec='white', mew=.3,
            label=r'$\overline{D}^* > 50$')
    ax.plot(cmp_['cv_exp'][~core], cmp_['cv_model'][~core], 'o', ms=3, mfc='white', mec=MODC, mew=.7,
            label='plume edge')
    ax.plot([0, 4], [0, 4], '-', c='#9AA3AD', lw=.6, zorder=0)
    ax.set_xlim(0, 4); ax.set_ylim(0, 4); ax.set_aspect('equal')
    ax.set_xlabel('CV, wind tunnel'); ax.set_ylabel('CV, SPH')
    ax.legend(loc='upper left', frameon=False, fontsize=6, handletextpad=.2, bbox_to_anchor=(0.0, 0.9))
    plabel(ax, '(d)')
    fig.tight_layout(w_pad=.8)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.38, 0.0), ncol=2, frameon=False, fontsize=6.5)
    save(fig, 'FigV4_cedval_variability')


# ── V5 산점도 ──────────────────────────────────────────────────────────
def fig_scatter():
    o, m, rec = cmp_['exp'], cmp_['model'] * K, cmp_['rec']
    g = rec[:, 2] < 3
    fig, ax = plt.subplots(figsize=(W1 * 0.75, W1 * 0.75))
    lo, hi = 1, 5000
    xx = np.array([lo, hi])
    ax.fill_between(xx, xx / 2, xx * 2, color='#E3E7EC', lw=0, zorder=0)
    ax.plot(xx, xx, '-', c='#9AA3AD', lw=.6, zorder=1)
    ax.plot(o[g], m[g], 'o', ms=3.2, mfc=MODC, mec='white', mew=.3, label='z = 2 m', zorder=3)
    ax.plot(o[~g], m[~g], '^', ms=3.4, mfc='white', mec=MODC, mew=.7, label='vertical profiles', zorder=3)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax.set_xlabel(r'wind tunnel $\overline{D}^*$'); ax.set_ylabel(r'SPH $\overline{D}^*$')
    ax.legend(loc='upper left', frameon=False, fontsize=6)
    save(fig, 'FigV5_cedval_scatter')


def metrics():
    o, m, rec = cmp_['exp'], cmp_['model'] * K, cmp_['rec']
    for lab, s in [('all', np.ones(len(o), bool)), ('z=2 m', rec[:, 2] < 3), ('vertical', rec[:, 2] >= 3)]:
        oo, pp = o[s], m[s]
        f2 = np.mean((pp / oo >= .5) & (pp / oo <= 2))
        fb = 2 * (oo.mean() - pp.mean()) / (oo.mean() + pp.mean())
        nm = np.mean((oo - pp) ** 2) / (oo.mean() * pp.mean())
        print(f'{lab:9s} n={s.sum():2d}  FAC2={f2:.2f}  FB={fb:+.2f}  NMSE={nm:.2f}')


if __name__ == '__main__':
    fig_layout(); fig_lateral(); fig_vertical(); fig_variability(); fig_scatter(); metrics()
