#!/usr/bin/env python3
"""논문 그림 5·6 — 캐노피 연직 프로파일(z0·d 적합)과 유동 기반 침적 예측자."""
import csv, json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'dep_paper_figs.py')).read().split('# ---------------- 데이터')[0].split('# ----------------')[0]
     .replace("import csv, os", "import csv, os"), globals())   # rcParams·팔레트·panel·save 재사용

prof = list(csv.DictReader(open(os.path.join(HERE, 'canopy_profile.csv'))))
z = np.array([float(r['z_m']) for r in prof]); U = np.array([float(r['U_mean']) for r in prof])
K = np.array([float(r['k_turb']) for r in prof])
J = json.load(open(os.path.join(HERE, 'canopy_aero.json')))
fit, mac, mo = J['fit'], J['macdonald'], J['morphology']
Hb = mo['H_built']

fig, axs = plt.subplots(1, 2, figsize=(W2*0.66, 62*MM))
ax = axs[0]
ax.plot(U, z, 'o', ms=2.4, c=OI['blue'], mec='none', label='LES 150-s mean')
zz = np.linspace(fit['fit_lo'], fit['fit_hi'], 50)
ax.plot(fit['ustar']/0.4*np.log((zz-fit['d'])/fit['z0']), zz, '-', c=OI['black'], lw=1,
        label=f"log law: $d$={fit['d']:.1f} m, $z_0$={fit['z0']:.2f} m")
ax.axhline(Hb, ls=':', c='#888', lw=.7)
ax.text(7.4, Hb*1.10, r'$H_\mathrm{built}$', fontsize=6, color='#555')
ax.axhspan(fit['fit_lo'], fit['fit_hi'], color='#f2f2f2', zorder=0, lw=0)
ax.set_xlabel(r'$\overline{U}$ (m s$^{-1}$)'); ax.set_ylabel('z (m)'); ax.set_ylim(0, 160)
ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.62))
ax2 = axs[1]
ax2.plot(K, z, 's', ms=2.4, c=OI['green'], mec='none')
ax2.axhline(Hb, ls=':', c='#888', lw=.7)
ax2.set_xlabel(r'$\overline{k}$ (m$^2$ s$^{-2}$)'); ax2.set_ylabel('z (m)'); ax2.set_ylim(0, 160)
for a, l in zip(axs, 'ab'): panel(a, l)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig5_canopy_profile')

# --- fig 6: 유동 기반 예측자 ---
rows = list(csv.DictReader(open(os.path.join(HERE, 'flow_predictors.csv'))))
ids = [int(r['src']) for r in rows]
COURT_ = {3, 6, 10}
y = np.array([float(r['f_dep']) for r in rows])
fig, axs = plt.subplots(1, 2, figsize=(W2*0.66, 55*MM))
for ax, key, xlab in [(axs[0], 'lambda_p75', r'$\lambda_\mathrm{p}$ (R = 75 m)'),
                      (axs[1], 'k50', r'$\overline{k}$ within 50 m, z < 10 m (m$^2$ s$^{-2}$)')]:
    x = np.array([float(r[key]) for r in rows])
    keep = np.array([i not in COURT_ for i in ids])
    ax.scatter(x[keep], y[keep], s=14, c=OI['blue'], lw=0, zorder=3)
    ax.scatter(x[~keep], y[~keep], s=16, facecolors='none', edgecolors=OI['red'], lw=.9, zorder=3)
    b, a = np.polyfit(x[keep], y[keep], 1)
    xs = np.linspace(x.min()*.95, x.max()*1.03, 20); ax.plot(xs, a + b*xs, '-', c='k', lw=.7)
    rs = np.corrcoef(np.argsort(np.argsort(x[keep])), np.argsort(np.argsort(y[keep])))[0, 1]
    ax.text(.04, .95, rf'$\rho_\mathrm{{s}} = {rs:+.2f}$', transform=ax.transAxes, va='top', fontsize=6.5)
    for i, xi, yi in zip(ids, x, y):
        ax.annotate(f'{i}', (xi, yi), fontsize=5.5, xytext=(3, 2.5), textcoords='offset points', color='#444')
    ax.set_xlabel(xlab); ax.set_ylabel('Deposited fraction')
for a, l in zip(axs, 'ab'): panel(a, l)
fig.tight_layout(w_pad=2.0)
save(fig, 'fig6_flow_predictor')
print('paper/fig5, fig6 생성')
