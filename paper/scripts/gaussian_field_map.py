#!/usr/bin/env python3
"""f44 — 보행자층 시간적분 농도장: 건물해상 LES 대 운영 가우시안 (같은 색눈금).

양쪽 모두 **단위 방출량당 시간적분 농도** [s m^-3] 로 맞춘다.
  LES      : ped_grid_fine.npz 의 체류(입자·s/셀)를 소스별 방출입자수로 나누고
             셀 부피(cell^2 x 2 m)로 나눈다.
  가우시안 : 지면농도 Cg(단위 방출량당) x 그 지점의 노출시간 texp.
             Cg·texp 는 gaussian_reference.py 의 침적식에서 Vd 를 뺀 것과 같다.
             (0~2 m 층평균 대신 지면값을 쓴다 — 관심 범위에서 sigma_z >> 2 m 라 차이가 수 %다.)

소스 13개, 중립, 풍향 0°, 전이 런 `2026-09-15_leipzig-dt01-sponge-uy`.
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
import figures as F
from figstyle import C, finish
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DEP  = '/home/kdh/Urban-SPH/runs/lambdap/deposition'
st.OUT = os.environ.get('FIG_OUT', os.path.join(HERE, '..', 'figs', 'main'))

# 가우시안 설정·식을 그대로 가져온다 (U, VD, sigmas, depletion, H_SRC, T0/T1/TEND)
g = {'__file__': os.path.join(HERE, 'gaussian_reference.py')}
exec(open(g['__file__']).read().split('# ---- LES 침적 지도')[0], g)
sigmas, depletion = g['sigmas'], g['depletion']
U, H_SRC, T0, T1, TEND = g['U'], g['H_SRC'], g['T0'], g['T1'], g['TEND']

# ── LES 농도장 ────────────────────────────────────────────────────────────
G = np.load(os.path.join(DEP, 'ped_grid_fine.npz'))
cell = float(G['cell']); lim = G['lim']; nrel = G['nrel']
ids = sorted(int(k[1:]) for k in G.files if k.startswith('S'))
vol = cell*cell*2.0
E_les = sum(G[f'S{i}']/nrel[k] for k, i in enumerate(ids))/vol
gx = lim[0]+(np.arange(E_les.shape[0])+.5)*cell
gy = lim[2]+(np.arange(E_les.shape[1])+.5)*cell

# ── 가우시안 농도장 (같은 격자, 같은 소스) ────────────────────────────────
S13 = {int(r['id']): (float(r['x']), float(r['y'])) for r in
       csv.DictReader(open(os.path.join(DEP, '..', 'sources13.csv')))}
assert sorted(S13) == ids, (sorted(S13), ids)
X, Y = np.meshgrid(gx, gy, indexing='ij')

def gauss_field(sx, sy, kind='urban'):
    """소스 하나의 보행자층 시간적분 농도 [s m^-3], 단위 방출량당."""
    out = np.zeros_like(E_les)
    dx, dy = X-sx, Y-sy
    ok = dx > 1.0
    sy_, sz_ = sigmas(dx[ok], kind)
    fd = depletion(np.maximum(gx-sx, 1.0), kind)               # x 에만 의존
    FD = np.broadcast_to(fd[:, None], X.shape)[ok]
    texp = np.clip(TEND-T0-dx[ok]/U, 0.0, T1-T0)
    Cg = FD/((T1-T0)*np.pi*U*sy_*sz_)*np.exp(-dy[ok]**2/(2*sy_**2))*np.exp(-H_SRC**2/(2*sz_**2))
    out[ok] = Cg*texp
    return out

E_gau = sum(gauss_field(sx, sy) for sx, sy in S13.values())

floor = 1.0/(float(nrel.min())*vol)        # LES 검출한계 = 입자 1개가 1초 머문 셀
vmax = np.log10(max(E_les.max(), E_gau.max()))
print(f"LES  총 {E_les.sum()*vol:.3f}  최대 {E_les.max():.3e} s/m3")
print(f"가우 총 {E_gau.sum()*vol:.3f}  최대 {E_gau.max():.3e} s/m3   (비 {E_gau.sum()/E_les.sum():.3f})")
print(f"검출한계 {floor:.2e} s/m3,  색눈금 10^{np.log10(floor):.1f} ~ 10^{vmax:.1f}")

# ── 그림 (f39 와 같은 두 패널 규격) ───────────────────────────────────────
PW, PH = 3.30, 3.20
L, GAP, CBG, CBW, R, B, T = 0.58, 0.26, 0.13, 0.16, 0.60, 0.50, 0.26
FW, FH = L+PW+GAP+PW+CBG+CBW+R, B+PH+T
fig = plt.figure(figsize=(FW, FH))
axes = [fig.add_axes([(L+k*(PW+GAP))/FW, B/FH, PW/FW, PH/FH]) for k in (0, 1)]
cax = fig.add_axes([(L+2*PW+GAP+CBG)/FW, B/FH, CBW/FW, PH/FH])
for ax, E, lab in zip(axes, (E_les, E_gau),
                      ('LES, building-resolved', 'Operational Gaussian (urban Briggs)')):
    F._buildings(ax)
    Em = np.ma.masked_where(E < floor, E)
    m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=F.MAPCMAP, vmin=np.log10(floor),
                      vmax=vmax, shading='gouraud', rasterized=True, zorder=2)
    for sx, sy in S13.values():
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=3.6, mew=.7, zorder=6)
    ax.set_aspect('equal'); ax.set_xlim(*F.MAP_XLIM); ax.set_ylim(*F.MAP_YLIM)
    ax.set_xlabel('x   (m)'); ax.set_title(lab, fontsize=8.5, pad=4)
axes[0].annotate('', xy=(360, 420), xytext=(240, 420),
                 arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.1), zorder=7)
axes[0].text(240, 432, 'wind', fontsize=7.5, zorder=7)
axes[1].plot([950, 1150], [-440, -440], '-', c=C['ink'], lw=1.8, zorder=7)
axes[1].text(1050, -432, '200 m', fontsize=7.5, ha='center', zorder=7)
axes[0].text(240, -452, '13 release points · neutral · 150 s release', fontsize=7.5, zorder=7)
axes[0].set_ylabel('y   (m)'); axes[1].tick_params(labelleft=False)
cb = fig.colorbar(m, cax=cax)
cb.set_label(r'time-integrated concentration below 2 m'
             '\n' r'per unit release   ($\log_{10}$,  s m$^{-3}$)', fontsize=7.5)
cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
finish(fig, axes[0], 'f44_concentration_field_les_vs_gaussian', minor=False, tight=False)


# ══════════════════════════════════════════════════════════════════════════
# f45 — 유형별 대표 소스 4곳, LES(위) 대 가우시안(아래). f34 의 분류·창 규격을 따른다.
# ══════════════════════════════════════════════════════════════════════════
lam75 = {}
for r in csv.DictReader(open(os.path.join(DEP, '..', 'lambda_p_by_R.csv'))):
    if int(float(r['R_m'])) == 75:
        lam75 = {int(k[1:]): float(v) for k, v in r.items() if k.startswith('S')}
PANELS = [(8, 'open'), (13, 'semi-open'), (4, 'street canyon'), (10, 'enclosed courtyard')]
HALF = 190.0

les_i = {i: G[f'S{i}']/nrel[ids.index(i)]/vol for i, _ in PANELS}
gau_i = {i: gauss_field(*S13[i]) for i, _ in PANELS}
vmax5 = np.log10(max(max(v.max() for v in les_i.values()),
                     max(v.max() for v in gau_i.values())))

TW, GAPX, GAPY = 1.50, 0.030, 0.030
LFT, TOP, BOT, CBG, CBW, RGT = 0.30, 0.26, 0.06, 0.10, 0.17, 0.52
FW = LFT + 4*TW + 3*GAPX + CBG + CBW + RGT
FH = TOP + 2*TW + GAPY + BOT
fig = plt.figure(figsize=(FW, FH))
hf = np.load(os.path.join(DEP, '..', 'heightfield.npz'))
Hb, bx0, by0, bc = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
bx = bx0 + np.arange(Hb.shape[1]+1)*bc; by = by0 + np.arange(Hb.shape[0]+1)*bc
from matplotlib.colors import LinearSegmentedColormap
BCMAP = LinearSegmentedColormap.from_list('b', [F.BLDG, F.BLDG])

for row, (fields, rlab) in enumerate(((les_i, 'LES'), (gau_i, 'Gaussian'))):
    for col, (i, tlab) in enumerate(PANELS):
        ax = fig.add_axes([(LFT+col*(TW+GAPX))/FW, (BOT+(1-row)*(TW+GAPY))/FH, TW/FW, TW/FH])
        sx, sy = S13[i]
        ax.pcolormesh(bx, by, np.where(Hb > 1, 1.0, np.nan), cmap=BCMAP,
                      shading='flat', rasterized=True, zorder=1)
        E = fields[i]
        Em = np.ma.masked_where(E < floor, E)
        m = ax.pcolormesh(gx, gy, np.log10(Em.T), cmap=F.MAPCMAP, vmin=np.log10(floor),
                          vmax=vmax5, shading='gouraud', rasterized=True, zorder=2)
        ax.plot(sx, sy, 'o', mfc='white', mec=C['ink'], ms=4.0, mew=.8, zorder=6)
        ax.set_aspect('equal'); ax.set_xlim(sx-HALF, sx+HALF); ax.set_ylim(sy-HALF, sy+HALF)
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_linewidth(.7)
        if row == 0:
            ax.set_title(f'{tlab}\nS{i}   $\\lambda_\\mathrm{{p}}$ = {lam75[i]:.2f}',
                         fontsize=7.5, pad=3, linespacing=1.3)
        if col == 0:
            ax.text(-0.045, 0.5, rlab, transform=ax.transAxes, rotation=90,
                    ha='right', va='center', fontsize=8.5)
        if row == 1 and col == 0:
            ax.plot([sx-150, sx-50], [sy-160, sy-160], '-', c=C['ink'], lw=1.6, zorder=7)
            ax.text(sx-100, sy-152, '100 m', fontsize=7, ha='center', color=C['ink'], zorder=7)
        if row == 0 and col == 0:
            ax.annotate('', xy=(sx-70, sy+155), xytext=(sx-160, sy+155),
                        arrowprops=dict(arrowstyle='-|>', color=C['ink'], lw=1.0), zorder=7)
            ax.text(sx-160, sy+163, 'wind', fontsize=7, zorder=7)
cax = fig.add_axes([(LFT+4*TW+3*GAPX+CBG)/FW, BOT/FH, CBW/FW, (2*TW+GAPY)/FH])
cb = fig.colorbar(m, cax=cax)
cb.set_label(r'time-integrated concentration below 2 m per unit release'
             '\n' r'($\log_{10}$,  s m$^{-3}$)', fontsize=7.5)
cax.tick_params(labelsize=7, width=.7, length=2.4); cb.outline.set_linewidth(.7)
finish(fig, fig.axes[0], 'f45_closeups_les_vs_gaussian_by_type', minor=False, tight=False)
for i, lab in PANELS:
    print(f"  S{i:<3}{lab:<20} LES 최대 {les_i[i].max():.2e}   가우시안 최대 {gau_i[i].max():.2e}"
          f"   비 {les_i[i].max()/gau_i[i].max():5.1f}배")
