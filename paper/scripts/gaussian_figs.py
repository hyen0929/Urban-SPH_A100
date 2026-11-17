#!/usr/bin/env python3
"""운영 가우시안 대비 그림 둘 (RESULTS.md §2).

  f42  소스 간 편차 — LES 대 가우시안. 각 계열을 자기 평균으로 규격화해 '퍼짐'만 본다.
       총량 자체는 설정으로 조정 가능하므로(§2.2) 비교 대상이 아니다.
  f43  36개 설정 민감도 — 침적속도는 총량비(x)만 움직이고 소스 간 CV(y)는 못 올린다.

입력: ../data/08_gaussian/gaussian_vs_les.csv, gaussian_sensitivity.csv
양식: runs/lambdap/deposition/figstyle.py (88 mm 단일단, 그림 하나에 그래프 하나)
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
from figstyle import C, newfig, finish, legend_inside, note, log_minor, expand

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', '08_gaussian')
st.OUT = os.environ.get('FIG_OUT', os.path.join(HERE, '..', 'figs', 'main'))

cv = lambda v: float(np.std(v, ddof=1)/np.mean(v))

# ───────────────────────────────────────────────── f42  소스 간 편차
R = [r for r in csv.DictReader(open(os.path.join(DATA, 'gaussian_vs_les.csv')))]
d0 = [r for r in R if r['deg'] == '0']
les = np.array([float(r['les']) for r in d0])
gu  = np.array([float(r['gauss_urban']) for r in d0])
gr  = np.array([float(r['gauss_rural']) for r in d0])
o   = np.argsort(les)                      # LES 값 순으로 소스를 세운다
x   = np.arange(1, len(o)+1)

# 점 50개가 세 줄이라 기본 규격으로는 축이 꽉 찬다. 이 그림만 마커·글씨를 한 단계 줄이고
# 높이는 집안 기본(66 mm)으로 되돌려 여백을 준다.
fig, ax = newfig()
SER = [(les, 'LES, building-resolved', C['blue'],   'o'),
       (gu,  'Gaussian, urban',        C['orange'], 's'),
       (gr,  'Gaussian, rural',        C['green'],  '^')]
for v, lab, col, mk in SER:
    ax.plot(x, v[o]/v.mean(), mk, color=col, mec='white', mew=0.4, ms=3.0,
            ls='-', lw=0.6, label=f'{lab}   CV {cv(v):.3f}', zorder=3)
ax.axhline(1.0, color=C['grey'], lw=0.5, ls=(0, (4, 3)), zorder=1)
# 긴 라벨은 글씨를 줄여도 축 폭을 다 먹는다 — 문구 자체를 짧게 한다.
ax.set_xlabel('release point, ranked by LES', fontsize=6.5)
ax.set_ylabel(r'$f_\mathrm{dep}\ /\ \overline{f_\mathrm{dep}}$'  '\n'  r'(each series by its own mean)', fontsize=6.5, linespacing=1.5)
ax.tick_params(labelsize=6)
ax.set_xlim(0, 50); ax.set_ylim(0.25, 2.25)
ax.set_xticks(np.arange(0, 51, 10)); ax.set_yticks(np.arange(0.5, 2.1, 0.5))
# 축 범위를 못 박았으므로 legend_inside(자리가 모자라면 축을 넓힌다)를 쓰지 않는다.
ax.legend(loc='upper left', bbox_to_anchor=(0.015, 0.985), ncol=1, fontsize=6,
          handlelength=1.1, handletextpad=0.35, labelspacing=0.3, borderpad=0.35)
finish(fig, ax, 'f42_source_spread_les_vs_gaussian')
print(f"f42: LES CV {cv(les):.3f} ({les.max()/les.min():.2f}×), "
      f"urban {cv(gu):.4f} ({gu.max()/gu.min():.3f}×), rural {cv(gr):.4f}")

# ───────────────────────────────────────────────── f43  36개 설정 민감도
S = [r for r in csv.DictReader(open(os.path.join(DATA, 'gaussian_sensitivity.csv')))]
LES_CV = 0.408                                   # 50소스 × 4풍향 (§2.3)
vds = sorted({r['vd'] for r in S}, key=lambda s: -float(s.split()[-1]))
RAMP = ['#20618D', '#5C9BC7', '#A9C7DD']         # 순차 1색조: Vd 큰→작은
fig, ax = newfig(h=62/25.4)
ax.axhline(LES_CV, color=C['red'], lw=1.0, zorder=2)
ax.axvline(1.0, color=C['grey'], lw=0.6, ls=(0, (4, 3)), zorder=1)
for k, vd in enumerate(vds):
    for fam, mk in (('urban', 's'), ('rural', '^')):
        g = [r for r in S if r['vd'] == vd and r['sigma'].startswith(fam)]
        if not g: continue
        ax.plot([float(r['ratio']) for r in g], [float(r['cv']) for r in g], mk,
                color=RAMP[k], mec='white', mew=0.5, ms=4.4, ls='none', zorder=3,
                label=f"$V_d$ = {vd.split()[-1]} m s$^{{-1}}$" if fam == 'urban' else None)
# 모양도 정보를 지므로 범례에 넣는다 (색 = 침적속도, 모양 = 분산계수 계열)
for mk, lab in (('s', 'urban $\\sigma$ (Briggs)'), ('^', 'rural $\\sigma$')):
    ax.plot([], [], mk, color=C['grey'], mec='white', mew=0.5, ms=4.4, ls='none', label=lab)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('total deposition, Gaussian / LES')
ax.set_ylabel('CV between release points')
expand(ax, fx=0.10, fy=0.14, logx=True, logy=True)
log_minor(ax, 'both') if False else (log_minor(ax, 'y'), log_minor(ax, 'x'))
legend_inside(ax, ncol=2, prefer='upper left')
rat = np.array([float(r['ratio']) for r in S]); cvs = np.array([float(r['cv']) for r in S])
note(ax, f"36 settings\ntotal spans {rat.max()/rat.min():.0f}×\n"
         f"CV stays {cvs.min():.3f}–{cvs.max():.3f}", fontsize=7, loc='lower right')
ax.annotate('LES, 50 sources × 4 wind directions', xy=(rat.min()*1.1, LES_CV),
            xytext=(0, 3), textcoords='offset points', fontsize=7, color=C['red'], va='bottom')
finish(fig, ax, 'f43_gaussian_sensitivity_ratio_vs_cv', minor=False)
print(f"f43: ratio {rat.min():.3f}~{rat.max():.3f} ({rat.max()/rat.min():.0f}×), "
      f"CV {cvs.min():.4f}~{cvs.max():.4f}, 서로 다른 CV 값 {len(set(cvs))}개")
