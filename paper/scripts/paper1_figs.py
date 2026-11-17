#!/usr/bin/env python3
"""1편(B&E) 4장 그림 — Leipzig 13 소스, 중립, MOST 지표, 벽·지붕 침적.

  Fig07  계산 영역과 13 개 방출점 (regime 별 기호)
  Fig08  반송류: (a) 순시 |ω| z=11.25 m  (b) 시간평균 |u_h| z=1.25 m  (c) 수평평균 연직 프로파일
  Fig09  보행자층(0–2 m) 시간적분 농도 — 13 소스 합
  Fig10  regime 대표 3 소스 클로즈업 (이류 / 채널 / 갇힘)
  Fig11  소스별 지표 대 λp(75) — 노출 3 + 침적 3
  Fig12  침적 경로: 지면 대 건물 외피 (개방 대 중정)
  Fig13  전이 대 시간평균: (a,b) 침적 지도  (c) 소스별 비

자료: paper/data/10_paper1/cache_*.npz (paper1_cache.py),
      runs/lambdap_v2/deposition/ped_grid_fine*.npz, sources13.csv, lambda_p_by_R.csv
출력: paper/figs/paper1/FigNN_*.pdf/.png (600 dpi)
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
import figures as F
from figstyle import C
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.transforms import Bbox

ROOT = '/home/kdh/Urban-SPH'
D10 = f'{ROOT}/paper/data/10_paper1'
V2 = f'{ROOT}/runs/lambdap_v2/deposition'
OUT = f'{ROOT}/paper/figs/paper1'
os.makedirs(OUT, exist_ok=True)

MM = 1/25.4
W1, W15, W2 = 90*MM, 140*MM, 190*MM          # Elsevier 단일 / 1.5 / 2 단
plt.rcParams.update({'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 8,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                     'lines.markersize': 3.5})

BLDG = '#D5DAE0'
BCMAP = LinearSegmentedColormap.from_list('b', [BLDG, BLDG])
MAPCMAP = F.MAPCMAP
U_REF, Z_REF, ALPHA = 10.0, 371.25, 0.19

# ── 공통 자료 ──────────────────────────────────────────────────────────
SRC = {int(r['id']): (float(r['x']), float(r['y'])) for r in
       csv.DictReader(open(f'{ROOT}/runs/lambdap/sources13.csv'))}
IDS = sorted(SRC)
LAM = {}
for r in csv.DictReader(open(f'{ROOT}/runs/lambdap/lambda_p_by_R.csv')):
    if int(float(r['R_m'])) == 75:
        LAM = {int(k[1:]): float(v) for k, v in r.items() if k.startswith('S')}

# regime — 5 m 노출 지도(그림 9)의 플룸 형태로 분류
REGIME = {**{i: 'advected' for i in (8, 9, 11, 13)},
          **{i: 'channelled' for i in (1, 3, 4, 5, 7, 12)},
          **{i: 'trapped' for i in (2, 6, 10)}}
RSTYLE = {'advected': (C['blue'], 'o'), 'channelled': (C['orange'], 's'),
          'trapped': (C['red'], '^')}

PG = np.load(f'{V2}/ped_grid_fine.npz')
PGB = np.load(f'{V2}/ped_grid_fine_B.npz')
PCELL, PLIM, NREL = float(PG['cell']), PG['lim'], PG['nrel']
PVOL = PCELL*PCELL*2.0
PGX = PLIM[0] + (np.arange(PG['S1'].shape[0])+.5)*PCELL
PGY = PLIM[2] + (np.arange(PG['S1'].shape[1])+.5)*PCELL


def expo(i, G=PG):
    """소스 i 의 보행자층 시간적분 농도, 방출 입자 1개당 [s m^-3]."""
    return G[f'S{i}']/NREL[IDS.index(i)]/PVOL


# ── 건물 윤곽: IBM 마커 1 m 래스터 + 닫기 (중정은 보존) ───────────────
_BP = None
def _markers():
    global _BP
    if _BP is None:
        ns = {'__file__': f'{ROOT}/runs/lambdap/deposition/dep_metrics.py'}
        exec(open(ns['__file__']).read().split('def load(')[0], ns)
        p, _ = ns['read_vtk'](f'{ROOT}/runs/2026-10-20_leipzig-s13-most-wall/out/building_0stp.vtk')
        _BP = p[p[:, 2] > 1.5]          # 지면 높이 마커 제외 — 넣으면 중정이 메워진다
    return _BP


def _sh(g):
    return [np.roll(np.roll(g, a, 0), b, 1) for a in (-1, 0, 1) for b in (-1, 0, 1)]


def footprint(x0, x1, y0, y1, cell=1.0):
    P = _markers()
    m = (P[:, 0] >= x0) & (P[:, 0] < x1) & (P[:, 1] >= y0) & (P[:, 1] < y1)
    nx, ny = int(round((x1-x0)/cell)), int(round((y1-y0)/cell))
    g = np.zeros((nx, ny), bool)
    g[((P[m, 0]-x0)/cell).astype(int).clip(0, nx-1), ((P[m, 1]-y0)/cell).astype(int).clip(0, ny-1)] = True
    for _ in range(2): g = np.logical_or.reduce(_sh(g))      # 닫기: 2 회 팽창 …
    for _ in range(2): g = np.logical_and.reduce(_sh(g))     # … 2 회 침식. 구멍 채우기는 안 한다(중정)
    return g, x0+np.arange(nx+1)*cell, y0+np.arange(ny+1)*cell


def buildings(ax, x0, x1, y0, y1, z=1):
    g, bx, by = footprint(x0, x1, y0, y1)
    ax.pcolormesh(bx, by, np.where(g.T, 1.0, np.nan), cmap=BCMAP, shading='flat',
                  rasterized=True, zorder=z)


def mapframe(ax, x0, x1, y0, y1, ticks=True):
    ax.set_aspect('equal'); ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    if not ticks:
        ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_linewidth(.6)


# 지도 위 글씨는 흰 바탕을 깐다 — 색장·건물 위에서도 읽히게
WBOX = dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85)


def wind(ax, x, y, L=90, fs=6.5):
    ax.annotate('', xy=(x+L, y), xytext=(x, y),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=.9), zorder=8)
    ax.text(x, y+0.02*(ax.get_ylim()[1]-ax.get_ylim()[0]), 'wind', fontsize=fs, zorder=8,
            bbox=WBOX)


def scalebar(ax, x, y, L, lab, fs=6.5):
    ax.plot([x, x+L], [y, y], '-', c=C['ink'], lw=1.4, zorder=8, solid_capstyle='butt')
    ax.text(x+L/2, y+0.015*(ax.get_ylim()[1]-ax.get_ylim()[0]), lab, fontsize=fs,
            ha='center', va='bottom', zorder=8, bbox=WBOX)


def plabel(ax, s, x=0.02, y=0.98, color=None, box=True):
    ax.text(x, y, s, transform=ax.transAxes, ha='left', va='top', fontsize=8,
            fontweight='bold', color=color or C['ink'], zorder=9, bbox=WBOX if box else None)


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print(name)


DOM = (200.0, 1160.0, -460.0, 460.0)


# ══════════════════════════════════════════════════════════════════════
def fig07():
    fig, ax = plt.subplots(figsize=(W15*0.82, W15*0.82*0.98))
    x0, x1, y0, y1 = DOM
    buildings(ax, x0, x1, y0, y1)
    for i in IDS:
        col, mk = RSTYLE[REGIME[i]]
        ax.plot(*SRC[i], mk, mfc=col, mec='white', mew=.6, ms=6.0, zorder=6)
        ax.text(SRC[i][0]+12, SRC[i][1]+10, str(i), fontsize=6.5, zorder=7, color=C['ink'])
    mapframe(ax, x0, x1, y0, y1)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    ax.set_title('(c)', loc='left', fontsize=9, fontweight='bold', pad=3)
    wind(ax, 230, 405); scalebar(ax, 940, -435, 200, '200 m')
    h = [plt.Line2D([], [], ls='none', marker=RSTYLE[k][1], mfc=RSTYLE[k][0], mec='white',
                    mew=.6, ms=5.5, label=k) for k in ('advected', 'channelled', 'trapped')]
    ax.legend(handles=h, loc='upper right', frameon=True, fontsize=6.5, handletextpad=.3,
              borderpad=.4, labelspacing=.3)
    save(fig, 'Fig07_domain_sources')


# ══════════════════════════════════════════════════════════════════════
def fig08():
    I = np.load(f'{D10}/cache_inst.npz'); M = np.load(f'{D10}/cache_meanfield.npz')
    x0, x1, y0, y1 = DOM
    fig = plt.figure(figsize=(W2, W2*0.40))
    axa = fig.add_axes([0.00, 0.08, 0.36, 0.86]); axb = fig.add_axes([0.37, 0.08, 0.36, 0.86])
    axc = fig.add_axes([0.82, 0.16, 0.17, 0.74])

    def pmap(ax, x, y, f, cmap, norm, lab):
        cell = 2.5
        nx, ny = int((x1-x0)/cell), int((y1-y0)/cell)
        G = np.full((nx, ny), np.nan)
        i = ((x-x0)/cell).astype(int); j = ((y-y0)/cell).astype(int)
        ok = (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
        G[i[ok], j[ok]] = f[ok]
        m = ax.pcolormesh(x0+np.arange(nx+1)*cell, y0+np.arange(ny+1)*cell, G.T,
                          cmap=cmap, norm=norm, shading='flat', rasterized=True, zorder=1)
        buildings(ax, x0, x1, y0, y1, z=2)
        mapframe(ax, x0, x1, y0, y1, ticks=False)
        cb = fig.colorbar(m, ax=ax, orientation='horizontal', fraction=.045, pad=.03, aspect=35)
        cb.set_label(lab, fontsize=7); cb.ax.tick_params(labelsize=6.5, width=.6, length=2)
        cb.outline.set_linewidth(.6)

    pmap(axa, I['z11.25_x'], I['z11.25_y'], I['z11.25_w'], 'magma_r',
         LogNorm(0.05, 5), r'instantaneous $|\omega|$ (s$^{-1}$),  z = 11.25 m')
    uh = np.hypot(M['z1.25_f'][:, 0], M['z1.25_f'][:, 1])
    pmap(axb, M['z1.25_x'], M['z1.25_y'], uh/U_REF, MAPCMAP.reversed() if hasattr(MAPCMAP, 'reversed') else MAPCMAP,
         plt.Normalize(0, 0.35), r'time-averaged $|u_h|/U_\mathrm{ref}$,  z = 1.25 m')
    plabel(axa, '(a)'); plabel(axb, '(b)')
    wind(axa, 212, 25, L=70); scalebar(axb, 940, -430, 200, '200 m')

    P = M['prof']; z, u = P[:, 0], P[:, 1]
    axc.plot(u/U_REF, z, '-', c=C['blue'], lw=1.1, zorder=3, label='LES')
    zz = np.linspace(1, 371, 300)
    axc.plot((zz/Z_REF)**ALPHA, zz, ':', c=C['grey'], lw=.9, zorder=2, label='inflow')
    axc.axhline(18.9, c=C['ink'], lw=.5, ls=(0, (2, 2)))
    axc.text(0.02, 22, r'$\bar H$', fontsize=6.5)
    axc.set_yscale('log'); axc.set_ylim(1, 380); axc.set_xlim(0, 1.05)
    axc.set_xlabel(r'$\langle u_h\rangle/U_\mathrm{ref}$'); axc.set_ylabel('z (m)')
    axc.legend(loc='upper left', fontsize=6, frameon=False, handlelength=1.6)
    plabel(axc, '(c)', x=0.70, y=0.20)
    save(fig, 'Fig08_carrier_flow')


def fit_loglaw(zlo=40.0, zhi=150.0):
    P = np.load(f'{D10}/cache_meanfield.npz')['prof']; z, u = P[:, 0], P[:, 1]
    m = (z >= zlo) & (z <= zhi); best = None
    for d in np.arange(0.0, 18.01, 0.1):
        zz = z[m]-d
        A = np.polyfit(np.log(zz), u[m], 1)
        r = np.sqrt(np.mean((u[m]-np.polyval(A, np.log(zz)))**2))
        if best is None or r < best[0]: best = (r, d, A)
    r, d, A = best
    return dict(d=d, us=0.4*A[0], z0=float(np.exp(-A[1]/A[0])), rmse=r)


# ══════════════════════════════════════════════════════════════════════
def fig09():
    E = sum(expo(i) for i in IDS)
    one = 1.0/(NREL.max()*PVOL)
    x0, x1, y0, y1 = DOM
    fig, ax = plt.subplots(figsize=(W15, W15*0.80))
    Em = np.ma.masked_where(E < one, E)
    m = ax.pcolormesh(PGX, PGY, Em.T, cmap=MAPCMAP, norm=LogNorm(one, E.max()),
                      shading='gouraud', rasterized=True, zorder=1)
    buildings(ax, x0, x1, y0, y1, z=2)
    for i in IDS:
        ax.plot(*SRC[i], 'o', mfc='white', mec=C['ink'], ms=3.6, mew=.7, zorder=6)
        ax.text(SRC[i][0]+10, SRC[i][1]+10, str(i), fontsize=6.3, zorder=7)
    mapframe(ax, x0, x1, y0, y1)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    wind(ax, 230, 405); scalebar(ax, 940, -435, 200, '200 m')
    cb = fig.colorbar(m, ax=ax, fraction=.040, pad=.02)
    cb.set_label(r'time-integrated concentration, $z \leq 2$ m (s m$^{-3}$ per released tracer)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5, width=.6, length=2); cb.outline.set_linewidth(.6)
    save(fig, 'Fig09_pedestrian_exposure_all')


# ══════════════════════════════════════════════════════════════════════
REP = [(8, 'advected'), (4, 'channelled'), (10, 'trapped')]
def _win(i):
    sx, sy = SRC[i]; return (sx-90, sx+330, sy-210, sy+210)


def fig10():
    one = 1.0/(NREL.max()*PVOL)
    vmax = max(expo(i).max() for i, _ in REP)
    fig = plt.figure(figsize=(W2, W2*0.36))
    w = 0.295
    axes = [fig.add_axes([0.005+k*(w+0.012), 0.06, w, 0.88]) for k in range(3)]
    cax = fig.add_axes([0.005+3*(w+0.012)+0.004, 0.10, 0.012, 0.80])
    for k, (ax, (i, lab)) in enumerate(zip(axes, REP)):
        x0, x1, y0, y1 = _win(i)
        E = expo(i); Em = np.ma.masked_where(E < one, E)
        m = ax.pcolormesh(PGX, PGY, Em.T, cmap=MAPCMAP, norm=LogNorm(one, vmax),
                          shading='gouraud', rasterized=True, zorder=1)
        buildings(ax, x0, x1, y0, y1, z=2)
        ax.plot(*SRC[i], 'o', mfc='white', mec=C['ink'], ms=4.2, mew=.8, zorder=6)
        mapframe(ax, x0, x1, y0, y1, ticks=False)
        plabel(ax, f'({"abc"[k]}) source {i}, {lab}')
        if k == 0: wind(ax, x0+20, y1-60, L=70)
        if k == 2: scalebar(ax, x1-120, y0+20, 100, '100 m')
    cb = fig.colorbar(m, cax=cax)
    cb.set_label('time-integrated concentration\n' r'$z \leq 2$ m (s m$^{-3}$)', fontsize=7, linespacing=1.3)
    cb.ax.tick_params(labelsize=6.5, width=.6, length=2); cb.outline.set_linewidth(.6)
    save(fig, 'Fig10_regimes')


# ══════════════════════════════════════════════════════════════════════
def metrics():
    L = np.load(f'{D10}/cache_ldm.npz')
    xyz, src, bt = L['A_xyz'], L['A_src'], L['A_bt']
    out = {}
    GX, GY = np.meshgrid(PGX, PGY, indexing='ij')
    for i in IDS:
        E = expo(i); A_ = PG[f'S{i}']
        r = np.hypot(GX-SRC[i][0], GY-SRC[i][1])
        m = src == i; g = (m & (bt == 1)).sum(); w = (m & (bt == 2)).sum()
        dep = m & (bt > 0)
        dx = xyz[dep, 0]; dy = xyz[dep, 1]
        out[i] = dict(D=float(A_.sum()/NREL[IDS.index(i)]), peak=float(E.max()),
                      area=float((E > 0.01*E.max()).sum()*PCELL*PCELL),
                      rmean=float((A_*r).sum()/A_.sum()),
                      fdep=float((g+w)/m.sum()), env=float(w/max(g+w, 1)),
                      ndep=int(g+w))
    return out


def sp(a, b):
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def perm(a, b, n=20000, seed=1):
    rng = np.random.default_rng(seed); o = abs(sp(a, b)); b = np.asarray(b)
    return float(np.mean([abs(sp(a, rng.permutation(b))) >= o for _ in range(n)]))


def fig11():
    Mx = metrics()
    lam = np.array([LAM[i] for i in IDS])
    panels = [('D', 'residence per tracer, $z\\leq2$ m (s)', False),
              ('area', 'exposed area, $>1\\%$ of peak (m$^2$)', True),
              ('peak', 'peak concentration (s m$^{-3}$)', True),
              ('fdep', 'deposited fraction', True),
              ('env', 'envelope share of deposition', False),
              ('rmean', 'mean distance of exposure (m)', False)]
    fig, axs = plt.subplots(2, 3, figsize=(W2, W2*0.56))
    stats = {}
    for k, (ax, (key, lab, logy)) in enumerate(zip(axs.flat, panels)):
        y = np.array([Mx[i][key] for i in IDS])
        if key == 'env': y = y*100; lab = lab + ' (%)'
        if key == 'fdep': y = y*100; lab = 'deposited fraction (%)'
        for i, yy in zip(IDS, y):
            col, mk = RSTYLE[REGIME[i]]
            ax.plot(LAM[i], yy, mk, mfc=col, mec='white', mew=.5, ms=5.0, zorder=3)
        if logy: ax.set_yscale('log')
        ax.set_xlim(-0.02, 0.50); ax.set_xlabel(r'$\lambda_p$ (R = 75 m)'); ax.set_ylabel(lab)
        r = sp(lam, y); p = perm(lam, y); stats[key] = (r, p)
        ax.margins(y=0.10)
        # 패널 라벨·상관계수는 축 밖 제목줄에 — 축 안에 두면 점과 겹친다
        ax.set_title(f'({"abcdef"[k]})', loc='left', fontsize=8, fontweight='bold', pad=3)
        ax.set_title(rf'$\rho_s$ = {r:+.2f}', loc='right', fontsize=7, pad=3)
    h = [plt.Line2D([], [], ls='none', marker=RSTYLE[k][1], mfc=RSTYLE[k][0], mec='white',
                    mew=.5, ms=5, label=k) for k in ('advected', 'channelled', 'trapped')]
    fig.legend(handles=h, loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02),
               fontsize=7, handletextpad=.3, columnspacing=1.5)
    fig.tight_layout(rect=(0, 0, 1, 0.97), w_pad=1.2, h_pad=1.0)
    save(fig, 'Fig11_metrics_vs_lambda')
    return Mx, stats


# ══════════════════════════════════════════════════════════════════════
def fig12(pair=((8, 'open'), (6, 'enclosed courtyard'))):
    L = np.load(f'{D10}/cache_ldm.npz'); xyz, src, bt = L['A_xyz'], L['A_src'], L['A_bt']
    fig = plt.figure(figsize=(W15, W15*0.54))
    w = 0.47
    axes = [fig.add_axes([0.01+k*(w+0.02), 0.10, w, 0.84]) for k in range(2)]
    for k, (ax, (i, lab)) in enumerate(zip(axes, pair)):
        sx, sy = SRC[i]
        if k == 0: x0, x1, y0, y1 = sx-50, sx+350, sy-200, sy+200   # 개방: 풍하 쪽으로 창을 민다
        else:      x0, x1, y0, y1 = sx-140, sx+140, sy-140, sy+140
        buildings(ax, x0, x1, y0, y1)
        m = src == i; g = m & (bt == 1); wl = m & (bt == 2)
        ax.plot(xyz[g, 0], xyz[g, 1], '.', c=C['blue'], ms=0.9, alpha=0.30, mec='none', rasterized=True, zorder=2)
        ax.plot(xyz[wl, 0], xyz[wl, 1], '.', c=C['orange'], ms=1.6, alpha=0.9, mec='none', rasterized=True, zorder=3)
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.2, mew=.8, zorder=6)
        mapframe(ax, x0, x1, y0, y1, ticks=False)
        plabel(ax, f'({"ab"[k]}) source {i}, {lab}')
        if k == 0: wind(ax, x0+15, y1-55, L=70)
        scalebar(ax, x1-120, y0+15, 100, '100 m')
    h = [plt.Line2D([], [], ls='none', marker='o', ms=3.5, mec='none', color=c, label=t)
         for c, t in ((C['blue'], 'deposited on the ground'), (C['orange'], 'deposited on walls and roofs'))]
    fig.legend(handles=h, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01),
               fontsize=7, handletextpad=.3, columnspacing=1.6)
    save(fig, 'Fig12_deposition_pathway')


# ══════════════════════════════════════════════════════════════════════
def fig13(i=11):
    L = np.load(f'{D10}/cache_ldm.npz')
    sx, sy = SRC[i]; x0, x1, y0, y1 = sx-100, sx+400, sy-150, sy+250
    cell = 5.0
    nx, ny = int((x1-x0)/cell), int((y1-y0)/cell)
    grids = []
    for tag in ('A', 'B'):
        xyz, src, bt = L[f'{tag}_xyz'], L[f'{tag}_src'], L[f'{tag}_bt']
        m = (src == i) & (bt > 0)
        Hh, _, _ = np.histogram2d(xyz[m, 0], xyz[m, 1], bins=[nx, ny], range=[[x0, x1], [y0, y1]])
        grids.append(Hh/(cell*cell))
    vmax = max(g.max() for g in grids); vmin = 1/(cell*cell)
    fig = plt.figure(figsize=(W2, W2*0.33))
    w = 0.285
    axa = fig.add_axes([0.00, 0.08, w, 0.86]); axb = fig.add_axes([w+0.01, 0.08, w, 0.86])
    cax = fig.add_axes([2*w+0.018, 0.12, 0.010, 0.78])
    axc = fig.add_axes([0.745, 0.20, 0.25, 0.72])
    for ax, Gd, lab in ((axa, grids[0], '(a) transient'), (axb, grids[1], '(b) time-averaged')):
        Gm = np.ma.masked_where(Gd <= 0, Gd)
        mm = ax.pcolormesh(x0+np.arange(nx+1)*cell, y0+np.arange(ny+1)*cell, Gm.T, cmap=MAPCMAP,
                           norm=LogNorm(vmin, vmax), shading='flat', rasterized=True, zorder=1)
        buildings(ax, x0, x1, y0, y1, z=2)
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.0, mew=.8, zorder=6)
        mapframe(ax, x0, x1, y0, y1, ticks=False); plabel(ax, f'{lab}, source {i}')
    wind(axa, x0+20, y1-95, L=70); scalebar(axb, x1-120, y0+15, 100, '100 m')
    cb = fig.colorbar(mm, cax=cax); cb.set_label('deposited tracers (m$^{-2}$)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5, width=.6, length=2); cb.outline.set_linewidth(.6)
    # (c) 소스별 비
    ts = list(csv.DictReader(open(f'{V2}/dep_timeseries.csv')))
    fin = {(r['run'], int(r['src'])): r for r in ts if int(r['step']) == 30000}
    rd = np.array([int(fin[('B', s)]['n_dep'])/int(fin[('A', s)]['n_dep']) for s in IDS])
    ra = np.array([float(fin[('B', s)]['area90_m2'])/float(fin[('A', s)]['area90_m2']) for s in IDS])
    o = np.argsort(rd); y = np.arange(len(IDS))
    axc.axvline(1.0, c=C['grey'], lw=.6, ls=(0, (3, 2)))
    axc.plot(rd[o], y, 'o', mfc=C['blue'], mec='white', mew=.5, ms=4.5, label='deposited count', zorder=3)
    axc.plot(ra[o], y, 's', mfc=C['orange'], mec='white', mew=.5, ms=4.0, label='footprint area', zorder=3)
    axc.set_yticks(y); axc.set_yticklabels([str(IDS[k]) for k in o], fontsize=6.5)
    axc.set_ylabel('source'); axc.set_xlabel('time-averaged / transient')
    axc.set_xlim(0.2, 1.5); axc.set_ylim(-0.7, len(IDS)-0.3)
    axc.legend(loc='lower right', fontsize=6.3, frameon=True, handletextpad=.3, borderpad=.35)
    plabel(axc, '(c)', x=0.03, y=0.99)
    save(fig, 'Fig13_transient_vs_mean')


FIT = None
if __name__ == '__main__':
    todo = sys.argv[1:] or ['fig07', 'fig08', 'fig09', 'fig10', 'fig11', 'fig12', 'fig13']
    for n in todo:
        r = globals()[n]()
        if n == 'fig11':
            Mx, stats = r
            for k, (rr, pp) in stats.items(): print(f'   {k:6s} rho_s={rr:+.3f}  p={pp:.4f}')
            np.save(f'{D10}/metrics13.npy', Mx, allow_pickle=True)
