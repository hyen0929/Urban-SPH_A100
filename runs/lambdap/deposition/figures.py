#!/usr/bin/env python3
"""낱장 그림 생성 — paper/figs/ 에 그래프 하나당 파일 하나 (PDF + 600 dpi PNG).

  python3 figures.py            전부
  python3 figures.py f05 f12    일부
"""
import csv, json, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
import figstyle as st
# 같은 그림을 다른 런에 그대로 다시 그리기 위한 경로 오버라이드 (기본값 = 중립 기준런)
PED_CSV    = os.environ.get('PED_CSV',    'ped_exposure.csv')
PED_BANDS  = os.environ.get('PED_BANDS',  'ped_exposure_bands.csv')
PED_FINE   = os.environ.get('PED_FINE',   'ped_grid_fine.npz')
PED_FINE_B = os.environ.get('PED_FINE_B', 'ped_grid_fine_B.npz')
PED_GRID   = os.environ.get('PED_GRID',   'ped_grid.npz')
from figstyle import C, MM, newfig, finish, legend_inside, note, log_minor, label_points, sources_scatter, expand

HERE = os.path.dirname(os.path.abspath(__file__))
if os.environ.get('FIG_OUT'): st.OUT = os.environ['FIG_OUT']
MORPH = os.path.normpath(os.path.join(HERE, '..'))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
COURT = {3, 6, 10}
R_MAIN = 150
XL = rf'$\lambda_\mathrm{{p}}$   (R = {R_MAIN} m)'

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']

lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
IDS = sorted(S); KEEP = np.array([i not in COURT for i in IDS])
X = np.array([lam[R_MAIN][i] for i in IDS])
met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
corr = list(csv.DictReader(open(os.path.join(HERE, 'dep_lambda_corr.csv'))))
ped = {int(r['src']): {k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(HERE, PED_CSV)))}
bands = list(csv.DictReader(open(os.path.join(HERE, PED_BANDS))))
flow = {int(r['src']): {k: float(v) for k, v in r.items()} for r in csv.DictReader(open(os.path.join(HERE, 'flow_predictors.csv')))}


def sp(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


# 지도용: 흰 바탕 → 밝은 노랑 → 주황 → 자주 (검은 바탕 대신)
MAPCMAP = LinearSegmentedColormap.from_list(
    'expo', ['#FFF9E8', '#FDD98A', '#F3924B', '#D94E2F', '#8C2150', '#3F1246'])
BLDG = '#DCE1E6'


def scatter_fig(name, y, ylab, log=False, loc='upper right', xvals=None, xlab=None, legend=True):
    fig, ax = newfig()
    xv = X if xvals is None else np.asarray(xvals, float)
    y = np.asarray(y, float)
    if log: ax.set_yscale('log')
    ax.set_xlim(xv.min(), xv.max()); ax.set_ylim(y.min(), y.max())
    expand(ax, logy=log)
    b, a = np.polyfit(xv[KEEP], np.log10(y[KEEP]) if log else y[KEEP], 1)
    xs = np.linspace(*ax.get_xlim(), 50)
    ax.plot(xs, 10**(a+b*xs) if log else a+b*xs, '-', c=C['grey'], lw=0.9, zorder=3)
    sources_scatter(ax, xv, y, IDS)
    rs = sp(xv[KEEP], np.log10(y[KEEP]) if log else y[KEEP])
    ax.set_xlabel(xlab or XL); ax.set_ylabel(ylab)
    if log: log_minor(ax)
    # 범례부터 구석에 앉히고, ρs 는 남은 빈 공간 한가운데에 자연스럽게 놓는다
    if legend: legend_inside(ax, ncol=2)
    note(ax, rf'$\rho_\mathrm{{s}} = {rs:+.2f}$   (n = 10)')
    finish(fig, ax, name)


# ---------------------------------------------------------------- 침적 형상
def f01(): scatter_fig('f01_footprint_vs_lambda', [met['A'][i]['area90_m2'] for i in IDS],
                       r'Footprint area $A_{90}$   (m$^2$)', log=True, loc='upper right')


def f02(): scatter_fig('f02_trapped_vs_lambda', [met['A'][i]['f_r150'] for i in IDS],
                       r'Trapped fraction   ($r \leq 150$ m)', loc='lower right')


def f03(): scatter_fig('f03_deposited_vs_lambda', [met['A'][i]['f_dep'] for i in IDS],
                       'Deposited fraction', loc='upper right')


def f04():
    fig, ax = newfig(st.W_ONE*1.12, st.H_ONE)
    style = [('area90_m2', 'Footprint area', C['blue'], 'o'), ('f_r150', 'Trapped fraction', C['green'], 's'),
             ('shift_m', 'Centroid shift', C['purple'], '^'), ('f_dep', 'Deposited fraction', C['orange'], 'D'),
             ('elong', 'Elongation', C['sky'], 'v')]
    for mt, lab, c, m in style:
        rr = sorted((int(float(r['R'])), float(r['spearman_nocourt'])) for r in corr
                    if r['metric'] == mt and r['run'] == 'A')
        ax.plot([a for a, _ in rr], [b for _, b in rr], marker=m, c=c, ms=4.2, lw=1.1,
                mec='white', mew=.5, label=lab)
    ax.axhspan(-1.05, 1.05, xmin=0, xmax=0, color='none')
    ax.axvspan(75, 100, color='#EFF2F4', zorder=0, lw=0)
    ax.axhline(0, c=C['grey'], lw=0.7, zorder=1)
    ax.set_xlim(0, 320); ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel('Window radius  R   (m)')
    ax.set_ylabel(r'Spearman $\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
    legend_inside(ax, ncol=3)
    finish(fig, ax, 'f04_radius_scan')


def f05():
    ts = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'dep_timeseries.csv'))):
        if r['run'] == 'A':
            ts.setdefault(int(r['src']), []).append((float(r['time']), int(r['n_released']), int(r['n_dep'])))
    for k in ts: ts[k].sort()
    hi = {8: ('S8  open', C['blue']), 13: ('S13  semi-open', C['green']),
          4: ('S4  dense', C['orange']), 10: ('S10  courtyard', C['red'])}
    fig, ax = newfig(st.W_ONE*1.06, st.H_ONE)
    for i in IDS:
        a = np.array([(t, n/max(rel, 1)) for t, rel, n in ts[i]])
        if i in hi: continue
        gl, = ax.plot(a[:, 0], a[:, 1], lw=0.8, c='#C2C9D0', zorder=2)
    for i, (lab, col) in hi.items():
        a = np.array([(t, n/max(rel, 1)) for t, rel, n in ts[i]])
        ax.plot(a[:, 0], a[:, 1], lw=1.4, c=col, zorder=4, label=lab)
    ax.set_xlim(150, 302); ax.set_ylim(0, .43)
    ax.set_xlabel('Time   (s)'); ax.set_ylabel('Deposited fraction')
    # 회색 다발도 범례 항목으로. ncol=2 는 열 우선이라 5번째가 1열 맨 아랫줄에 온다.
    h, l = ax.get_legend_handles_labels()
    gl.set_label('other nine sources')
    h = h[:2] + [gl] + h[2:]; l = l[:2] + ['other nine sources'] + l[2:]
    legend_inside(ax, ncol=2, fontsize=7, handles=h, labels=l)
    finish(fig, ax, 'f05_deposition_timeseries')



# ---------------------------------------------------------------- 지도 공통 규격
# 네 장의 지도(침적·노출 × 세기·소스별)가 같은 도메인을 그리므로 축 박스와 우측 띠의
# 크기를 인치 단위로 고정한다. subplots + colorbar(fraction=..) 는 제목·라벨 길이에 따라
# 박스가 흔들려서 쓰지 않는다.
MAP_FIG = (5.90, 4.50)                       # 그림 전체 [in]
MAP_AX = (0.62, 0.52, 3.92, 3.80)            # 축 박스 [in]  (도메인 종횡비 960/930)
MAP_CB = (4.70, 0.52, 0.17, 3.80)            # 우측 띠 [in]
MAP_XLIM, MAP_YLIM = (200, 1160), (-465, 465)


def _mapfig():
    fig = plt.figure(figsize=MAP_FIG)
    fw, fh = MAP_FIG
    ax = fig.add_axes([MAP_AX[0]/fw, MAP_AX[1]/fh, MAP_AX[2]/fw, MAP_AX[3]/fh])
    cax = fig.add_axes([MAP_CB[0]/fw, MAP_CB[1]/fh, MAP_CB[2]/fw, MAP_CB[3]/fh])
    return fig, ax, cax


def _buildings(ax):
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    bx = x0 + np.arange(H.shape[1]+1)*hc; by = y0 + np.arange(H.shape[0]+1)*hc
    ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                  cmap=LinearSegmentedColormap.from_list('b', [BLDG, BLDG]),
                  shading='flat', rasterized=True, zorder=1)


def _mapdeco(ax, scale=True):
    ax.set_aspect('equal'); ax.set_xlim(*MAP_XLIM); ax.set_ylim(*MAP_YLIM)
    ax.annotate('', xy=(360, 420), xytext=(240, 420),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.1), zorder=7)
    ax.text(240, 432, 'wind', fontsize=7.5, zorder=7)
    if scale:
        ax.plot([950, 1150], [-440, -440], '-', c=C['ink'], lw=1.8, zorder=7)
        ax.text(1050, -432, '200 m', fontsize=7.5, ha='center', zorder=7)
    ax.set_xlabel('x   (m)'); ax.set_ylabel('y   (m)')


# 소스 13개용 파스텔 고리. 이웃한 번호끼리 구분되도록 채도·명도를 번갈아 준다.
def _src_colors(n=13):
    import colorsys
    out = []
    for k in range(n):
        h = (k*7 % n)/n                        # 7 과 13 은 서로소 → 인접 번호끼리 색상이 멀어진다
        sat, val = (0.38, 0.99) if k % 2 == 0 else (0.52, 0.90)
        out.append(colorsys.hsv_to_rgb(h, sat, val))
    return out


def _src_map(grids, name, cblabel, gamma=0.55):
    """셀마다 가장 많이 기여한 소스의 색으로 칠하고, 세기는 흰색→파스텔→진한색 램프로."""
    import colorsys
    from matplotlib.colors import to_rgb
    cols = _src_colors(len(IDS))
    lim, cell = grids['lim'], float(grids['cell'])
    Gs = np.stack([grids[f'S{i}'] for i in IDS])          # (13, nx, ny)
    tot = Gs.sum(0); dom = Gs.argmax(0)
    hi = np.percentile(tot[tot > 0], 99.5)
    t = np.clip(np.log10(tot+1)/np.log10(hi+1), 0, 1)**gamma
    rgb = np.ones(tot.shape+(3,), float)
    for k, base in enumerate(cols):
        h, sv, vv = colorsys.rgb_to_hsv(*base)
        deep = colorsys.hsv_to_rgb(h, min(1, sv*2.0), vv*0.62)
        m = (dom == k) & (tot > 0)
        if not m.any(): continue
        tt = t[m][:, None]
        pale = np.array(base)[None, :]; dp = np.array(deep)[None, :]
        # 0→흰색, 0.5→파스텔, 1→진한색
        low = 1.0 + (pale-1.0)*np.clip(tt/0.5, 0, 1)
        rgb[m] = np.where(tt < 0.5, low, pale + (dp-pale)*np.clip((tt-0.5)/0.5, 0, 1))
    alpha = (tot > 0).astype(float)
    img = np.dstack([rgb, alpha])
    fig, ax, cax = _mapfig()
    _buildings(ax)
    ax.imshow(np.transpose(img, (1, 0, 2)), origin='lower',
              extent=[lim[0], lim[0]+Gs.shape[1]*cell, lim[2], lim[2]+Gs.shape[2]*cell],
              interpolation='nearest', zorder=2)
    for k, i in enumerate(IDS):
        sx, sy = S[i]
        ax.plot(sx, sy, 'o', mfc=cols[k], mec=C['ink'], ms=5.0, mew=.9, zorder=6)
    _mapdeco(ax)
    label_points(ax, [S[i][0] for i in IDS], [S[i][1] for i in IDS],
                 [str(i) for i in IDS], fontsize=7)
    # 우측 띠 = 소스 색 칩 (세기 지도의 컬러바와 같은 자리·같은 크기)
    cax.set_xlim(0, 1); cax.set_ylim(0, len(IDS)); cax.set_xticks([])
    for k, i in enumerate(IDS):
        cax.add_patch(plt.Rectangle((0, len(IDS)-1-k), 1, 1, facecolor=cols[k],
                                    edgecolor='white', lw=.6))
    cax.set_yticks(np.arange(len(IDS))+0.5)
    cax.set_yticklabels([f'S{i}' for i in IDS[::-1]], fontsize=6.5)
    cax.yaxis.tick_right(); cax.tick_params(length=0)
    for sp_ in cax.spines.values(): sp_.set_linewidth(.7)
    cax.set_ylabel(cblabel, fontsize=7.5)
    cax.yaxis.set_label_position('right'); cax.yaxis.labelpad = 16
    finish(fig, ax, name, minor=False, tight=False)


def f29():
    _src_map(np.load(os.path.join(HERE, 'grid_dep_src.npz')),
             'f29_depmap_by_source', 'dominant source   ·   ground deposition')


def f30():
    fp = os.path.join(HERE, PED_FINE)
    _src_map(np.load(fp), 'f30_exposure_by_source',
             'dominant source   ·   exposure below 2 m')


# ---------------------------------------------------------------- 침적 지도
DEPGRID = 'grid_dep5.npz'          # 5 m 격자 (노출 지도와 같은 규격). 없으면 10 m 로 되돌아간다


def _depgrid():
    fp = os.path.join(HERE, DEPGRID)
    return np.load(fp if os.path.exists(fp) else os.path.join(HERE, 'grid_dep.npz'))


def _depmap(G, name, lab=None, vmax=None, diff=False):
    """침적 지도. 제목은 쓰지 않는다 — 네 장의 지도가 같은 축 박스를 갖도록."""
    d = _depgrid()
    lim = d['lim']
    gx = np.linspace(lim[0], lim[1], G.shape[0]+1); gy = np.linspace(lim[2], lim[3], G.shape[1]+1)
    fig, ax, cax = _mapfig()
    _buildings(ax)
    if diff:
        m = ax.pcolormesh(gx, gy, G.T, cmap='RdBu_r', vmin=-1.2, vmax=1.2,
                          shading='flat', rasterized=True, zorder=2)
        cl = 'transient − mean flow   ($\\log_{10}$)'
    else:
        Gm = np.ma.masked_where(G <= 0, G)
        m = ax.pcolormesh(gx, gy, np.log10(Gm.T+1), cmap=MAPCMAP, vmin=0, vmax=vmax,
                          shading='flat', rasterized=True, zorder=2)
        cell = float(d['cell'])
        cl = rf'deposits per {cell*cell:.0f} m$^2$   ($\log_{{10}}(N+1)$)'
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(cl, fontsize=7.5); cax.tick_params(labelsize=7, width=.7, length=2.4)
    cb.outline.set_linewidth(.7)
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.4, mew=.8, zorder=6)
    _mapdeco(ax)
    finish(fig, ax, name, minor=False, tight=False)


def f06():
    d = _depgrid()
    vmax = np.log10(max(d['GA'].max(), d['GB'].max())+1)
    _depmap(d['GA'], 'f06_depmap_transient', 'Transient LES', vmax)


def f07():
    d = _depgrid()
    vmax = np.log10(max(d['GA'].max(), d['GB'].max())+1)
    _depmap(d['GB'], 'f07_depmap_meanflow', 'Frozen mean flow', vmax)


def f08():
    d = _depgrid()
    D = np.log10(d['GA']+1)-np.log10(d['GB']+1)
    D = np.ma.masked_where((d['GA'] <= 0) & (d['GB'] <= 0), D)   # 양쪽 다 침적 없는 셀은 비운다
    _depmap(D, 'f08_depmap_difference', 'Difference', diff=True)


def _ab(name, key, lab, log=False):
    a = np.array([met['A'][i][key] for i in IDS]); b = np.array([met['B'][i][key] for i in IDS])
    fig, ax = newfig()
    if log: ax.set_xscale('log'); ax.set_yscale('log')
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    pad = (1.35 if log else 1.12)
    lo, hi = (lo/pad, hi*pad) if log else (lo-(hi-lo)*.12, hi+(hi-lo)*.12)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.plot([lo, hi], [lo, hi], '-', c=C['grey'], lw=0.8, zorder=2, dashes=(6, 3))
    sources_scatter(ax, a, b, IDS)
    ax.set_xlabel(f'{lab},  transient'); ax.set_ylabel(f'{lab},  mean flow')
    ax.set_aspect('equal', adjustable='box')
    if log: log_minor(ax); log_minor(ax, 'x')
    legend_inside(ax, ncol=1)
    note(ax, f'median ratio  {np.median(b/a):.2f}', fontsize=7.5)
    finish(fig, ax, name)


def f09(): _ab('f09_AB_deposited_fraction', 'f_dep', 'Deposited fraction')
def f10(): _ab('f10_AB_footprint_area', 'area90_m2', r'$A_{90}$  (m$^2$)', log=True)


# ---------------------------------------------------------------- 캐노피
# f11(캐노피 풍속 프로파일)은 2026-09-17 사용자 요청으로 제외.
# 번호는 그대로 두었다 — 캡션·보고서와 대응을 유지하기 위함.


def f12():
    prof = list(csv.DictReader(open(os.path.join(HERE, 'canopy_profile.csv'))))
    z = np.array([float(r['z_m']) for r in prof]); K = np.array([float(r['k_turb']) for r in prof])
    Hb = json.load(open(os.path.join(HERE, 'canopy_aero.json')))['morphology']['H_built']
    fig, ax = newfig()
    ax.set_ylim(0, 160); ax.set_xlim(0, K.max()*1.08)
    ax.plot(K, z, 's', ms=3.0, c=C['green'], mec='white', mew=.4, zorder=4)
    ax.axhline(Hb, ls=(0, (5, 2.5)), c=C['grey'], lw=0.8)
    ax.annotate(r'$H_\mathrm{built}$', (K.max()*1.03, Hb+4.5), fontsize=7.5, color=C['grey'], ha='right')
    ax.set_xlabel(r'$\overline{k}$   (m$^2$ s$^{-2}$)'); ax.set_ylabel('z   (m)')
    finish(fig, ax, 'f12_canopy_tke_profile')


def f13():
    scatter_fig('f13_deposited_vs_local_tke', [met['A'][i]['f_dep'] for i in IDS], 'Deposited fraction',
                xvals=[flow[i]['k50'] for i in IDS],
                xlab=r'$\overline{k}$ within 50 m, z < 10 m   (m$^2$ s$^{-2}$)', loc='upper left')


# ---------------------------------------------------------------- 보행자 노출
def _expo_map(i, name):
    fp = os.path.join(HERE, PED_FINE)
    G = np.load(fp if os.path.exists(fp) else os.path.join(HERE, PED_GRID))
    lim, cell = G['lim'], float(G['cell'])
    E = G[f'S{i}']/(cell*cell*2.0)
    vmax = max(np.log10(G[f'S{j}'].max()/(cell*cell*2.0)+1) for j in (8, 13, 4, 10))
    sx, sy = S[i]; half = 200
    gx = np.linspace(lim[0], lim[1], E.shape[0]+1); gy = np.linspace(lim[2], lim[3], E.shape[1]+1)
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    fig, ax = plt.subplots(figsize=(st.W_ONE*0.92, st.W_ONE*0.92))
    bx = x0 + np.arange(H.shape[1]+1)*hc; by = y0 + np.arange(H.shape[0]+1)*hc
    ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                  cmap=LinearSegmentedColormap.from_list('b', [BLDG, BLDG]), shading='flat', rasterized=True)
    Em = np.ma.masked_where(E <= 0, E)
    m = ax.pcolormesh(gx, gy, np.log10(Em.T+1), cmap=MAPCMAP, vmin=0, vmax=vmax, shading='flat', rasterized=True)
    cb = fig.colorbar(m, ax=ax, fraction=.042, pad=.02)
    cb.set_label(r'exposure   ($\log_{10}$, particle s m$^{-3}$)', fontsize=7.5)
    cb.ax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    ax.plot(sx, sy, 'o', mfc='none', mec=C['ink'], ms=5, mew=.9, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(sx-half, sx+half); ax.set_ylim(sy-half, sy+half)
    ax.set_xticks([sx-150, sx, sx+150]); ax.set_xticklabels(['−150', '0', '150'])
    ax.set_yticks([sy-150, sy, sy+150]); ax.set_yticklabels(['−150', '0', '150'])
    ax.set_xlabel(r'$x - x_\mathrm{src}$   (m)'); ax.set_ylabel(r'$y - y_\mathrm{src}$   (m)')
    # 제목 없음 — 네 장의 축 박스를 같게 두고, 분류 이름은 슬라이드 캡션이 맡는다
    finish(fig, ax, name, minor=False, tight=False)


def f14(): _expo_map(8, 'f14_exposure_S08_open')
def f15(): _expo_map(13, 'f15_exposure_S13_semiopen')
def f16(): _expo_map(4, 'f16_exposure_S04_dense')
def f17(): _expo_map(10, 'f17_exposure_S10_courtyard')


def f18(): scatter_fig('f18_exposure_area_vs_lambda', [ped[i]['area_1pct_m2'] for i in IDS],
                       r'Exposure area   (m$^2$)', log=True, loc='upper right')


def f19(): scatter_fig('f19_exposure_peak_vs_lambda', [ped[i]['peak_s_per_m3'] for i in IDS],
                       r'Peak exposure   (particle s m$^{-3}$)', loc='upper left')


def f20(): scatter_fig('f20_exposure_total_vs_lambda', [ped[i]['exposure_s_per_particle'] for i in IDS],
                       'Exposure per released particle   (s)', loc='upper left')


def f21():
    fig, ax = newfig(st.W_ONE*1.06, st.H_ONE)
    bl = ['0-1.5', '0-2', '0-4', '2-5', '5-10', '10-20']
    for key, lab, c, m in [('area_1pct_m2', 'Exposure area', C['blue'], 'o'),
                           ('frac_within_150m', 'Trapped fraction', C['green'], 's'),
                           ('exposure_s_per_particle', 'Total exposure', C['orange'], 'D')]:
        v = []
        for bb in bl:
            rows = sorted([r for r in bands if r['band'] == bb], key=lambda r: int(r['src']))
            y = np.array([float(r[key]) for r in rows])
            xs_ = np.array([lam[R_MAIN][int(r['src'])] for r in rows])
            k2 = np.array([int(r['src']) not in COURT for r in rows])
            v.append(sp(xs_[k2], y[k2]))
        ax.plot(range(len(bl)), v, marker=m, c=c, ms=4.2, lw=1.1, mec='white', mew=.5, label=lab)
    ax.axhline(0, c=C['grey'], lw=0.7)
    ax.axvspan(0.5, 2.5, color='#EFF2F4', zorder=0, lw=0)
    ax.set_xticks(range(len(bl))); ax.set_xticklabels(bl)
    ax.set_xlim(-.3, len(bl)-.7); ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel('Height band   (m)'); ax.set_ylabel(r'Spearman $\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
    legend_inside(ax, ncol=2, prefer='upper right')
    finish(fig, ax, 'f21_exposure_height_bands', minor=False)


# ---------------------------------------------------------------- 벽면 침적
WRUNS = [('2 m', '2026-09-17_leipzig-dt01-wall-meanflow', C['blue'], 'o'),
         ('3 m', '2026-09-17_leipzig-dt01-wall3m-meanflow', C['green'], 's'),
         ('5 m', '2026-09-17_leipzig-dt01-wall5m-meanflow', C['purple'], '^')]


def _wall():
    W = {}
    for lab, d, _, _ in WRUNS:
        p, F = read_vtk(os.path.join(RUNS, d, 'out', 'LDM_30000stp.vtk'))
        b = F['buffer_type'].astype(int); s = F['p_type'].astype(int)-100
        W[lab] = (np.array([int(((s == i) & (b == 1)).sum()) for i in IDS]),
                  np.array([int(((s == i) & (b == 2)).sum()) for i in IDS]), p[b == 2, 2])
    return W


def f22():
    W = _wall()
    fig, ax = newfig(st.W_ONE*1.06, st.H_ONE)
    ax.set_yscale('log'); ax.set_xlim(0.03, 0.55); ax.set_ylim(2.5e-4, 4.0)
    for lab, _, col, mk in WRUNS:
        g, w, _z = W[lab]
        ratio = np.maximum(w/np.maximum(g, 1), 6e-4)
        ax.plot(X[KEEP], ratio[KEEP], ls='', marker=mk, ms=4.6, mfc=col, mec='white', mew=.5,
                label=rf'$d_\mathrm{{w}}$ = {lab}')
        ax.plot(X[~KEEP], ratio[~KEEP], ls='', marker=mk, ms=5.2, mfc='white', mec=C['red'], mew=1.0)
    ax.set_xlabel(XL); ax.set_ylabel('wall + roof  /  ground deposition')
    log_minor(ax)
    legend_inside(ax, ncol=3, prefer='lower right')     # 범례는 아래, ρs 는 그 위
    note(ax, r'$\rho_\mathrm{s}$ = +0.84 to +0.89')
    finish(fig, ax, 'f22_wall_ratio_vs_lambda')


def f23():
    W = _wall(); _g, _w, z = W['2 m']
    fig, ax = newfig()
    ax.hist(z, bins=np.arange(0, 60, 2.5), color=C['blue'], edgecolor='white', lw=.5)
    ax.axvline(18.8, ls=(0, (5, 2.5)), c=C['grey'], lw=.8)
    ax.annotate(r'$H_\mathrm{built}$', (19.8, ax.get_ylim()[1]*.88), fontsize=7.5, color=C['grey'])
    ax.set_xlabel('Deposit height  z   (m)'); ax.set_ylabel('Wall + roof deposits')
    finish(fig, ax, 'f23_wall_deposit_height')


def f24():
    W = _wall(); N = 99988
    vals = [sp(X[KEEP], (W['2 m'][0]/N)[KEEP])] + [sp(X[KEEP], ((W[l][0]+W[l][1])/N)[KEEP]) for l, _, _, _ in WRUNS]
    labs = ['ground\nonly', 'incl. wall\n$d_w$ = 2 m', 'incl. wall\n$d_w$ = 3 m', 'incl. wall\n$d_w$ = 5 m']
    fig, ax = newfig()
    ax.bar(range(4), vals, color=[C['orange']]+[C['blue']]*3, width=.6, edgecolor='white', lw=.6)
    ax.axhline(0, c=C['ink'], lw=.7)
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=7.5)
    ax.set_ylim(min(vals)-0.09, 0.015); ax.set_ylabel(r'$\rho_\mathrm{s}$ (deposited fraction vs $\lambda_\mathrm{p}$)')
    for i_, v in enumerate(vals):
        ax.annotate(f'{v:+.2f}', (i_, v), xytext=(0, -9), textcoords='offset points', ha='center',
                    fontsize=7.5, color=C['ink'])
    finish(fig, ax, 'f24_wall_effect_on_correlation', minor=False)


# ---------------------------------------------------------------- 도메인 노출
def _domain(fine=True):
    # 지도는 고해상도 격자(ped_grid_fine.py), 셀 단위 통계는 표준 10 m 격자를 쓴다 —
    # 면적·상관 같은 지표는 셀 크기에 의존하므로 보고서 표와 같은 격자여야 한다.
    fp = os.path.join(HERE, PED_FINE)
    G = np.load(fp if (fine and os.path.exists(fp)) else os.path.join(HERE, PED_GRID))
    lim, cell = G['lim'], float(G['cell'])
    E = sum(G[f'S{i}'] for i in S)/(cell*cell*2.0)
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    return E, lim, cell, H, x0, y0, hc


def f25():
    E, lim, cell, H, x0, y0, hc = _domain()
    gx = lim[0]+(np.arange(E.shape[0])+.5)*cell; gy = lim[2]+(np.arange(E.shape[1])+.5)*cell
    fig, ax, cax = _mapfig()
    _buildings(ax)
    one = 1.0/(cell*cell*2.0)                        # 1 입자·s 에 해당하는 밀도
    Em = np.ma.masked_where(E < 0.5*one, E)
    m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=MAPCMAP,
                      vmin=np.log10(one), vmax=np.log10(E.max()),
                      shading='gouraud', rasterized=True, zorder=2)
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r'exposure below 2 m   ($\log_{10}$, particle s m$^{-3}$)', fontsize=7.5)
    cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.4, mew=.8, zorder=6)
    _mapdeco(ax)
    label_points(ax, [S[i][0] for i in IDS], [S[i][1] for i in IDS], [str(i) for i in IDS], fontsize=7)
    finish(fig, ax, 'f25_domain_exposure_map', minor=False, tight=False)


def f26():
    E, lim, cell, H, x0, y0, hc = _domain(fine=False)
    B = (H > 1.0).astype(np.float64)
    gx = lim[0]+(np.arange(E.shape[0])+.5)*cell; gy = lim[2]+(np.arange(E.shape[1])+.5)*cell
    GX, GY = np.meshgrid(gx, gy, indexing='ij')
    ix = np.floor((GX-x0)/hc).astype(int); iy = np.floor((GY-y0)/hc).astype(int)
    ok = (ix >= 0) & (ix < H.shape[1]) & (iy >= 0) & (iy < H.shape[0])
    bld = np.zeros(E.shape, bool); bld[ok] = H[iy[ok], ix[ok]] > 1.0
    dmin = np.full(E.shape, 1e9)
    for i, (sx, sy) in S.items():
        dmin = np.minimum(dmin, np.hypot(GX-sx, GY-sy))
    rad = int(R_MAIN/hc); yy, xx = np.mgrid[-rad:rad+1, -rad:rad+1]
    K = ((xx*xx+yy*yy)*hc*hc <= R_MAIN**2).astype(np.float64); K /= K.sum()
    Kp = np.zeros_like(B); Kp[:K.shape[0], :K.shape[1]] = K
    Kp = np.roll(np.roll(Kp, -(K.shape[0]//2), 0), -(K.shape[1]//2), 1)
    L = np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp), B.shape)
    LAM = np.full(E.shape, np.nan); LAM[ok] = L[iy[ok], ix[ok]]
    edges = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.60]); mids = (edges[:-1]+edges[1:])/2
    fig, ax = newfig(st.W_ONE*1.12, st.H_ONE*1.02)
    ax.set_yscale('log')
    for (lo, hi), col, mk in [((25, 75), C['orange'], 'o'), ((50, 150), C['blue'], 's'),
                              ((100, 200), C['green'], '^'), ((150, 300), C['purple'], 'D')]:
        sel = (dmin > lo) & (dmin <= hi) & (~bld) & (E > 0) & np.isfinite(LAM)
        lv, ev = LAM[sel], E[sel]
        med = [np.median(ev[(lv >= a) & (lv < b)]) if ((lv >= a) & (lv < b)).sum() >= 20 else np.nan
               for a, b in zip(edges[:-1], edges[1:])]
        ax.plot(mids, med, marker=mk, c=col, ms=4.4, lw=1.1, mec='white', mew=.5,
                label=rf'{lo}–{hi} m  ($\rho_\mathrm{{s}}$ {sp(lv, np.log10(ev)):+.2f})')
    ax.set_xlim(0.06, 0.55); ax.set_ylim(0.02, 1.4)
    ax.set_xlabel(rf'local $\lambda_\mathrm{{p}}$   (R = {R_MAIN} m)')
    ax.set_ylabel(r'Median exposure   (particle s m$^{-3}$)')
    log_minor(ax)
    legend_inside(ax, ncol=2, fontsize=7, title='Distance from nearest source',
                  title_fontsize=7.5)
    finish(fig, ax, 'f26_exposure_vs_lambda_by_distance')


def f27():
    """방출고 민감도 — 지표별 λp 상관이 방출고에 견디는가."""
    fp = os.path.join(HERE, 'release_height.csv')
    if not os.path.exists(fp):
        print('  release_height.csv 없음 — f27 생략'); return
    rows = list(csv.DictReader(open(fp)))
    for r in rows:
        for k in r:
            if k != 'run': r[k] = int(r[k]) if k == 'src' else float(r[k])
    zs = sorted({r['z_rel'] for r in rows})
    fig, ax = newfig(st.W_ONE*1.12, st.H_ONE)
    style = [('f_r150', '150 m 갇힘', 'Trapped fraction', C['green'], 's'),
             ('area90_m2', 'A90', r'Footprint area $A_{90}$', C['blue'], 'o'),
             ('f_wall', '벽면', 'Wall + roof deposition', C['purple'], '^'),
             ('f_dep', '침적률', 'Deposited fraction', C['orange'], 'D')]
    for key, _, lab, col, mk in style:
        v = []
        for z in zs:
            K = [r for r in rows if r['z_rel'] == z and r['src'] not in COURT]
            v.append(sp([r['lam'] for r in K], [r[key] for r in K]))
        ax.plot(zs, v, marker=mk, c=col, ms=4.4, lw=1.1, mec='white', mew=.5, label=lab)
    ax.axhline(0, c=C['grey'], lw=0.7, zorder=1)
    ax.axhspan(-0.63, 0.63, color='#EFF2F4', zorder=0, lw=0)   # n=10, p>0.05 인 띠
    ax.set_xscale('log'); ax.set_xticks(zs); ax.set_xticklabels([f'{z:g}' for z in zs])
    ax.set_xlim(2.0, 26); ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel('Release height   (m)')
    ax.set_ylabel(r'Spearman $\rho_\mathrm{s}$ vs $\lambda_\mathrm{p}$')
    legend_inside(ax, ncol=2, fontsize=7)
    note(ax, 'shading:  not significant\n(n = 10, p > 0.05)', fontsize=7)
    finish(fig, ax, 'f27_release_height_sensitivity', minor=False)


def f28():
    """평균장이 어느 규모까지 전이 유동을 재현하는가 (meanflow_scale.py)."""
    fp = os.path.join(HERE, 'meanflow_scale.csv')
    if not os.path.exists(fp):
        print('  meanflow_scale.csv 없음 — f28 생략'); return
    rows = list(csv.DictReader(open(fp)))
    for r in rows:
        for k in r:
            if k != 'case': r[k] = float(r[k])
    fig, ax = newfig(st.W_ONE*1.12, st.H_ONE)
    sty = {'5x5 array': ('5×5 array  (idealised)', C['orange'], 'o'),
           'Leipzig': ('Leipzig  (real city)', C['blue'], 's')}
    for case, (lab, col, mk) in sty.items():
        R = sorted([r for r in rows if r['case'] == case], key=lambda r: r['cell_per_building'])
        x = [r['cell_per_building'] for r in R]
        ax.plot(x, [r['spearman'] for r in R], marker=mk, c=col, ms=4.4, lw=1.2,
                mec='white', mew=.5, label=lab)
        ax.plot(x, [r['pearson_linear'] for r in R], ls=(0, (4, 2)), c=col, lw=1.0, alpha=.55)
    ax.set_xscale('log'); ax.set_xlim(0.055, 3.2); ax.set_ylim(0.3, 1.02)
    ax.set_xlabel('Grid cell / building size')
    ax.set_ylabel('Deposition map correlation')
    log_minor(ax, 'x')
    legend_inside(ax, ncol=1)
    note(ax, 'solid: rank correlation\ndashed: Pearson on counts', fontsize=7.5)
    finish(fig, ax, 'f28_meanflow_scale', minor=False)


def _lamfield(R=R_MAIN):
    """높이장에서 λp(R) 장을 FFT 합성곱으로 만든다 (morph 방법론과 동일)."""
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    B = (H > 1.0).astype(np.float64)                 # 높이 1 m 초과 = 건물
    rad = int(R/hc); yy, xx = np.mgrid[-rad:rad+1, -rad:rad+1]
    K = ((xx*xx+yy*yy)*hc*hc <= R*R).astype(np.float64); K /= K.sum()
    Kp = np.zeros_like(B); Kp[:K.shape[0], :K.shape[1]] = K
    Kp = np.roll(np.roll(Kp, -(K.shape[0]//2), 0), -(K.shape[1]//2), 1)
    L = np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp), B.shape)
    return np.clip(L, 0, 1), H, x0, y0, hc


def f31():
    """λp(R) 장 위의 13개 방출원.

    다른 지도 넷과 나란히 놓지 않는 단독 그림이라 공통 규격(_mapfig)을 쓰지 않는다.
    λp 장이 배경이므로 축을 **높이장 실제 범위**에 딱 맞춰야 한다 — 공통 축범위
    (x 200~1160, y ±465)는 높이장(x 200~1120, y −452.5~455)보다 넓어서 오른쪽과
    위아래에 칠해지지 않은 흰 띠가 남았다.
    """
    L, H, x0, y0, hc = _lamfield()
    x1, y1 = x0 + H.shape[1]*hc, y0 + H.shape[0]*hc
    fig, ax = plt.subplots(figsize=(5.2, 4.35))
    bx = x0 + np.arange(H.shape[1]+1)*hc; by = y0 + np.arange(H.shape[0]+1)*hc
    m = ax.pcolormesh(bx, by, L, cmap='BuPu', vmin=0, vmax=0.6,
                      shading='flat', rasterized=True, zorder=1)
    ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                  cmap=LinearSegmentedColormap.from_list('b', ['#5A6670', '#5A6670']),
                  shading='flat', rasterized=True, zorder=2, alpha=.30)
    cb = fig.colorbar(m, ax=ax, fraction=.040, pad=.015)
    cb.set_label(rf'$\lambda_\mathrm{{p}}$   (R = {R_MAIN} m)', fontsize=9)
    cb.ax.tick_params(labelsize=8, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    for i in IDS:
        sx, sy = S[i]
        if i in COURT:
            ax.plot(sx, sy, 'o', mfc='none', mec=C['red'], ms=7.4, mew=1.6, zorder=6)
        else:
            ax.plot(sx, sy, 'o', mfc=C['blue'], mec='white', ms=6.4, mew=1.0, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.annotate('', xy=(x0+160, y1-40), xytext=(x0+40, y1-40),
                arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.1), zorder=7)
    ax.text(x0+40, y1-32, 'wind', fontsize=8, zorder=7)
    ax.plot([x1-220, x1-20], [y0+25, y0+25], '-', c=C['ink'], lw=1.8, zorder=7)
    ax.text(x1-120, y0+33, '200 m', fontsize=8, ha='center', zorder=7)
    ax.set_xlabel('x   (m)'); ax.set_ylabel('y   (m)')
    label_points(ax, [S[i][0] for i in IDS], [S[i][1] for i in IDS],
                 [str(i) for i in IDS], fontsize=8.5, color='#31383F')
    finish(fig, ax, 'f31_lambda_field', minor=False)


def _component(mask, seed_ij):
    """seed 가 속한 4이웃 연결 성분 (scipy 없이 반복 팽창)."""
    cur = np.zeros_like(mask, bool); cur[seed_ij] = True
    while True:
        nxt = cur.copy()
        nxt[1:, :] |= cur[:-1, :]; nxt[:-1, :] |= cur[1:, :]
        nxt[:, 1:] |= cur[:, :-1]; nxt[:, :-1] |= cur[:, 1:]
        nxt &= mask
        if nxt.sum() == cur.sum(): return cur
        cur = nxt


def f32(src=10, half=130.0):
    """중정 판정 — 개방 셀의 연결 성분이 가로망과 끊겨 있는가."""
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    sx, sy = S[src]
    j0 = int((sx-half-x0)/hc); j1 = int((sx+half-x0)/hc)
    i0 = int((sy-half-y0)/hc); i1 = int((sy+half-y0)/hc)
    Hc = H[i0:i1, j0:j1]
    openm = Hc <= 1.0
    si, sj = int((sy-y0)/hc)-i0, int((sx-x0)/hc)-j0
    court = _component(openm, (si, sj))
    # 창 테두리의 개방 셀에서 출발한 성분 = 가로망(외부와 이어진 쪽)
    seed = np.zeros_like(openm); seed[[0, -1], :] = True; seed[:, [0, -1]] = True
    street = _component(openm, np.where(seed & openm))
    other = openm & ~street & ~court                 # 그 밖의 고립 성분 = 다른 중정들
    BLD = (0.741, 0.784, 0.820); STR = (0.973, 0.984, 0.992)
    OTH = (0.937, 0.898, 0.812); CRT = (0.949, 0.749, 0.690)
    rgb = np.ones(Hc.shape+(3,))
    rgb[~openm] = BLD; rgb[street] = STR; rgb[other] = OTH; rgb[court] = CRT
    fig, ax = plt.subplots(figsize=(MAP_AX[3]*1.02, MAP_AX[3]))
    ex = [x0+j0*hc, x0+j1*hc, y0+i0*hc, y0+i1*hc]
    ax.imshow(rgb, origin='lower', extent=ex, interpolation='nearest')
    ax.contour(np.linspace(ex[0], ex[1], Hc.shape[1]), np.linspace(ex[2], ex[3], Hc.shape[0]),
               court.astype(float), levels=[0.5], colors=[C['red']], linewidths=1.2)
    ax.plot(sx, sy, 'o', mfc='none', mec=C['red'], ms=7, mew=1.6, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(ex[0], ex[1]); ax.set_ylim(ex[2], ex[3])
    ax.set_xlabel('x   (m)'); ax.set_ylabel('y   (m)')
    from matplotlib.patches import Patch
    leg = ax.legend(handles=[Patch(fc=BLD, ec='none', label='building'),
                             Patch(fc=STR, ec='#9AA6B0', lw=.7, label='open, reaches the street network'),
                             Patch(fc=OTH, ec='#C9B98E', lw=.7, label='open, isolated'),
                             Patch(fc=CRT, ec=C['red'], lw=1.1,
                                   label=f'S{src} release point — isolated, '
                                         f'{court.sum()*hc*hc:,.0f} m$^2$')],
                    loc='upper left', fontsize=7.5,
                    title='escape only over the roofline', title_fontsize=7.5)
    leg._legend_box.align = 'left'
    finish(fig, ax, f'f32_courtyard_S{src:02d}', minor=False)


def f33():
    """소스별 λp(R) — 반경을 바꾸면 값도 순위도 달라진다. 색은 소스별 지도(f29·f30)와 동일."""
    import colorsys
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    cols = _src_colors(len(IDS))
    Rs = sorted(lam)
    fig, ax = newfig(st.W_ONE*1.22, st.H_ONE*1.05)
    band = ax.axvspan(75, 100, color='#EFF2F4', zorder=0, lw=0)
    prim = ax.axvline(R_MAIN, c=C['grey'], lw=0.9, ls=(0, (4, 2.5)), zorder=1)
    deeps, ends = [], []
    for k, i in enumerate(IDS):
        h, sv, vv = colorsys.rgb_to_hsv(*cols[k])
        deep = colorsys.hsv_to_rgb(h, min(1, sv*1.9), vv*0.72)     # 색상은 지도와 같게, 선만 진하게
        y = [lam[R][i] for R in Rs]
        ax.plot(Rs, y, '-', c=deep, lw=1.3, zorder=4,
                marker='o' if i in COURT else None, ms=3.6, mec='white', mew=.5)
        deeps.append(deep); ends.append(y[-1])
    ax.set_xlim(15, 345); ax.set_ylim(0, 0.68)
    ax.set_xlabel('Window radius  R   (m)')
    ax.set_ylabel(r'$\lambda_\mathrm{p}$')
    # 선 끝 번호: 세로로만 밀어내 겹치지 않게 (같은 x 에 13개라 일반 라벨 배치로는 안 된다)
    fig.canvas.draw()
    yd = np.array([ax.transData.transform((Rs[-1], v))[1] for v in ends], float)
    gap = 8.0*fig.dpi/72.0
    order = np.argsort(yd)                                       # 정렬은 한 번만
    for n in range(1, len(order)):
        prev, cur = order[n-1], order[n]
        yd[cur] = max(yd[cur], yd[prev] + gap)
    top = ax.get_window_extent().y1 - 2
    if yd.max() > top: yd -= (yd.max()-top)
    inv = ax.transData.inverted()
    for k, i in enumerate(IDS):
        ax.annotate(f'{i}', (Rs[-1], ends[k]), xytext=(Rs[-1]+14, inv.transform((0, yd[k]))[1]),
                    textcoords='data', fontsize=7, color=deeps[k], va='center', zorder=5,
                    arrowprops=dict(arrowstyle='-', color=deeps[k], lw=.4, shrinkA=1, shrinkB=1))
    ax.legend(handles=[Patch(fc='#EFF2F4', ec='none', label='75–100 m:  peak correlation'),
                       Line2D([], [], c=C['grey'], lw=0.9, ls=(0, (4, 2.5)),
                              label=f'R = {R_MAIN} m:  primary radius'),
                       Line2D([], [], c=C['grey'], lw=0, marker='o', ms=3.6,
                              mfc='none', mec=C['grey'], label='enclosed courtyard')],
              loc='upper right', fontsize=7)
    finish(fig, ax, 'f33_lambda_vs_R')


def f34():
    """분류별 대표 소스 4곳의 보행자 노출 클로즈업 — 슬라이드용 2×2 타일.

    축·눈금을 없애고 분류 이름을 타일 안에 넣는다. 작게 넣어도 플룸 '모양'이 읽히도록.
    네 패널은 색 범위를 공유한다.
    """
    fp = os.path.join(HERE, PED_FINE)
    G = np.load(fp if os.path.exists(fp) else os.path.join(HERE, PED_GRID))
    lim, cell = G['lim'], float(G['cell'])
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    panels = [(8, 'open'), (13, 'semi-open'), (4, 'street canyon'), (10, 'enclosed courtyard')]
    half = 190.0
    E = {i: G[f'S{i}']/(cell*cell*2.0) for i, _ in panels}
    one = 1.0/(cell*cell*2.0)
    vmax = max(np.log10(v.max()) for v in E.values())
    fig = plt.figure(figsize=(4.95, 4.15))
    W, Hf = 4.95, 4.15
    tw = 1.92/W; th = 1.92/Hf
    pos = [(0.015, 0.505), (0.015+tw+0.012, 0.505), (0.015, 0.505-th-0.012),
           (0.015+tw+0.012, 0.505-th-0.012)]
    gx = lim[0]+(np.arange(E[8].shape[0])+.5)*cell; gy = lim[2]+(np.arange(E[8].shape[1])+.5)*cell
    bx = x0 + np.arange(H.shape[1]+1)*hc; by = y0 + np.arange(H.shape[0]+1)*hc
    for (i, lab), (px, py) in zip(panels, pos):
        ax = fig.add_axes([px, py, tw, th])
        sx, sy = S[i]
        ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                      cmap=LinearSegmentedColormap.from_list('b', [BLDG, BLDG]),
                      shading='flat', rasterized=True, zorder=1)
        Em = np.ma.masked_where(E[i] < 0.5*one, E[i])
        m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=MAPCMAP, vmin=np.log10(one), vmax=vmax,
                          shading='gouraud', rasterized=True, zorder=2)
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.0, mew=.8, zorder=6)
        ax.set_aspect('equal'); ax.set_xlim(sx-half, sx+half); ax.set_ylim(sy-half, sy+half)
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_linewidth(.7)
        ax.text(0.035, 0.965, f'{lab}\nS{i}   $\\lambda_\\mathrm{{p}}$ = {lam[R_MAIN][i]:.2f}',
                transform=ax.transAxes, ha='left', va='top', fontsize=7.5, color=C['ink'],
                linespacing=1.35, zorder=7)
        if i == 10:                                   # 축척은 한 패널에만
            ax.plot([sx+40, sx+140], [sy-160, sy-160], '-', c=C['ink'], lw=1.6, zorder=7)
            ax.text(sx+90, sy-152, '100 m', fontsize=7, ha='center', color=C['ink'], zorder=7)
    cax = fig.add_axes([0.015+2*tw+0.030, 0.505-th-0.012, 0.030, 2*th+0.012])
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r'exposure below 2 m   ($\log_{10}$, particle s m$^{-3}$)', fontsize=7.5)
    cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    os.makedirs(st.OUT, exist_ok=True)
    from matplotlib.transforms import Bbox
    kw = dict(bbox_inches=Bbox([[0, 0], fig.get_size_inches()]), pad_inches=0)
    fig.savefig(os.path.join(st.OUT, 'f34_closeups_by_type.pdf'), **kw)
    fig.savefig(os.path.join(st.OUT, 'f34_closeups_by_type.png'), dpi=600, **kw)
    plt.close(fig)


def f35():
    """전이 대 동결 평균장 — 분류별 대표 4곳의 보행자 노출 클로즈업 2×4.

    침적은 시간적분이라 차이가 뭉개진다(소스별 순위상관 0.72~0.88). 공중 노출장이
    비정상성에 직접 반응하므로 차이가 드러난다. 색 범위는 8개 패널 공통.
    """
    A = np.load(os.path.join(HERE, PED_FINE))
    B = np.load(os.path.join(HERE, PED_FINE_B))
    lim, cell = A['lim'], float(A['cell'])
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, hc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    # 전이/평균장 차이가 가장 큰 네 곳(순위상관 최저). 공교롭게 분류가 하나씩이다.
    panels = [(11, 'open'), (13, 'semi-open'), (7, 'street canyon'), (6, 'enclosed courtyard')]
    half, one = 190.0, 1.0/(cell*cell*2.0)
    E = {(g, i): (A if g == 0 else B)[f'S{i}']/(cell*cell*2.0) for g in (0, 1) for i, _ in panels}
    vmax = max(np.log10(v.max()) for v in E.values())
    k0 = panels[0][0]
    gx = lim[0]+(np.arange(E[(0, k0)].shape[0])+.5)*cell
    gy = lim[2]+(np.arange(E[(0, k0)].shape[1])+.5)*cell
    bx = x0 + np.arange(H.shape[1]+1)*hc; by = y0 + np.arange(H.shape[0]+1)*hc

    def sp_(a, b):
        return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])

    W, Hf = 8.35, 3.72
    tw, th = 1.62, 1.62
    x_lab, x0t, gap = 0.34, 0.60, 0.07
    y_rows = [(Hf-0.30-th), (Hf-0.30-2*th-0.06)]                 # 위=전이, 아래=평균장
    fig = plt.figure(figsize=(W, Hf))
    for gi, rowlab in enumerate(['Transient LES', 'Mean flow']):
        fig.text(x_lab/W, (y_rows[gi]+th/2)/Hf, rowlab, rotation=90, ha='left', va='center',
                 fontsize=9, color=C['ink'], fontweight='bold')
        for k, (i, lab) in enumerate(panels):
            px = x0t + k*(tw+gap)
            ax = fig.add_axes([px/W, y_rows[gi]/Hf, tw/W, th/Hf])
            sx, sy = S[i]
            ax.pcolormesh(bx, by, np.where(H > 1, 1.0, np.nan),
                          cmap=LinearSegmentedColormap.from_list('b', [BLDG, BLDG]),
                          shading='flat', rasterized=True, zorder=1)
            V = E[(gi, i)]
            Em = np.ma.masked_where(V < 0.5*one, V)
            m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=MAPCMAP, vmin=np.log10(one), vmax=vmax,
                              shading='gouraud', rasterized=True, zorder=2)
            ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=3.4, mew=.7, zorder=6)
            ax.set_aspect('equal'); ax.set_xlim(sx-half, sx+half); ax.set_ylim(sy-half, sy+half)
            ax.set_xticks([]); ax.set_yticks([])
            for sp2 in ax.spines.values(): sp2.set_linewidth(.7)
            if gi == 0:
                fig.text(px/W+tw/W/2, (y_rows[0]+th+0.06)/Hf, f'{lab}   S{i}',
                         ha='center', va='bottom', fontsize=8.5, color=C['ink'])
            else:
                a, b = E[(0, i)].ravel(), E[(1, i)].ravel()
                ok = (a+b) > 0
                ax.text(0.04, 0.05, f'rank ρ = {sp_(a[ok], b[ok]):.2f}', transform=ax.transAxes,
                        fontsize=8, color=C['red'], ha='left', va='bottom', zorder=7)
    cax = fig.add_axes([(x0t+4*tw+3*gap+0.13)/W, y_rows[1]/Hf, 0.16/W, (2*th+0.06)/Hf])
    cb = fig.colorbar(m, cax=cax)
    cb.set_label(r'exposure below 2 m   ($\log_{10}$, particle s m$^{-3}$)', fontsize=8)
    cax.tick_params(labelsize=7.5, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    os.makedirs(st.OUT, exist_ok=True)
    from matplotlib.transforms import Bbox
    kw = dict(bbox_inches=Bbox([[0, 0], fig.get_size_inches()]), pad_inches=0)
    fig.savefig(os.path.join(st.OUT, 'f35_transient_vs_mean_closeups.pdf'), **kw)
    fig.savefig(os.path.join(st.OUT, 'f35_transient_vs_mean_closeups.png'), dpi=600, **kw)
    plt.close(fig)


def f36():
    """도메인 전체 보행자 노출 차이 — 전이 − 평균장 (log 비). f08 과 같은 지도 규격."""
    A = np.load(os.path.join(HERE, PED_FINE))
    B = np.load(os.path.join(HERE, PED_FINE_B))
    lim, cell = A['lim'], float(A['cell'])
    EA = sum(A[f'S{i}'] for i in S)/(cell*cell*2.0)
    EB = sum(B[f'S{i}'] for i in S)/(cell*cell*2.0)
    one = 1.0/(cell*cell*2.0)
    # 한쪽만 0 이면 로그 차가 발산하므로 1 입자·s 밀도를 더해 비를 유계로 만든다
    D = np.log10((EA + 0.5*one)/(EB + 0.5*one))
    D = np.ma.masked_where((EA < 0.5*one) & (EB < 0.5*one), D)
    gx = lim[0]+(np.arange(EA.shape[0])+.5)*cell; gy = lim[2]+(np.arange(EA.shape[1])+.5)*cell
    fig, ax, cax = _mapfig()
    _buildings(ax)
    m = ax.pcolormesh(gx, gy, D.T, cmap='RdBu_r', vmin=-1.0, vmax=1.0,
                      shading='gouraud', rasterized=True, zorder=2)
    cb = fig.colorbar(m, cax=cax)
    cb.set_label('transient − mean flow   ($\\log_{10}$ exposure ratio)', fontsize=7.5)
    cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
    for i, (sx, sy) in S.items():
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.4, mew=.8, zorder=6)
    _mapdeco(ax)
    finish(fig, ax, 'f36_exposure_difference', minor=False, tight=False)


def f37():
    """침적 질량중심 이동의 풍하·횡풍 분해 — λp 가 억제하는 것은 풍하 이류다."""
    alo, cro = [], []
    for i in IDS:
        d = np.radians(met['A'][i]['shift_deg']); sm = met['A'][i]['shift_m']
        alo.append(sm*np.cos(d)); cro.append(abs(sm*np.sin(d)))
    alo, cro = np.array(alo), np.array(cro)
    fig, ax = newfig(st.W_ONE*1.12, st.H_ONE)
    ax.axhline(0, c=C['grey'], lw=0.7, zorder=1)
    b, a = np.polyfit(X[KEEP], alo[KEEP], 1)
    xs = np.linspace(0.05, 0.53, 50)
    ax.plot(xs, a+b*xs, '-', c=C['blue'], lw=0.9, zorder=3)
    for v, col, mk, lab in [(alo, C['blue'], 'o', 'along-wind'),
                            (cro, C['orange'], 's', 'cross-wind')]:
        ax.plot(X[KEEP], v[KEEP], ls='', marker=mk, ms=4.6, mfc=col, mec='white', mew=.5,
                label=lab, zorder=5)
        ax.plot(X[~KEEP], v[~KEEP], ls='', marker=mk, ms=5.2, mfc='white', mec=C['red'],
                mew=1.0, zorder=5)
    ax.set_xlim(0.04, 0.54); ax.set_ylim(-18, 92)
    ax.set_xlabel(XL); ax.set_ylabel('Centroid displacement   (m)')
    legend_inside(ax, ncol=2)
    note(ax, r'along-wind  $\rho_\mathrm{s} = -0.81$  (p = 0.006)' + '\n'
             r'cross-wind  $\rho_\mathrm{s} = -0.52$  (p = 0.13)')
    finish(fig, ax, 'f37_centroid_components')


if __name__ == '__main__':
    names = [n for n in sorted(globals()) if len(n) == 3 and n[0] == 'f' and n[1:].isdigit()]
    todo = sys.argv[1:] or names
    for n in todo:
        globals()[n]()
        print(f'{n} 완료', flush=True)
