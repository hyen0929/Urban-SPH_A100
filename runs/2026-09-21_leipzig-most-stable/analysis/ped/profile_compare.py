#!/usr/bin/env python3
"""중립 대 안정 성층 — 공중 LDM 입자의 연직 체류 프로파일.

노출 정의(z<=2 m 적분)와 같은 표본(t=150-300 s, 1 s 간격)을 쓰되 고도로 분해한다.
  R(z) dz = 방출 입자 1개가 고도 [z, z+dz) 에 머문 시간 [s/m]
  R 을 0-2 m 에 적분하면 ped_exposure.csv 의 exposure_s_per_particle 총합과 같다.
"""
import os, sys, re
import numpy as np
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
if os.environ.get('FIG_OUT'): st.OUT = os.environ['FIG_OUT']
from figstyle import C, newfig, finish, legend_inside, note, log_minor

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = '/home/kdh/Urban-SPH/runs'
ns = {'__file__': '/home/kdh/Urban-SPH/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def load(')[0], ns)
read = ns['read_vtk']

EVERY, DT, NREL = 100, 1.0, 13*99988
EDGES = np.arange(0.0, 60.001, 0.25)
CEN = 0.5*(EDGES[:-1]+EDGES[1:])
DZ  = EDGES[1]-EDGES[0]

def profile(run):
    d = os.path.join(RUNS, run, 'out')
    steps = sorted(int(re.search(r'_(\d+)stp', f).group(1))
                   for f in os.listdir(d) if f.startswith('LDM_'))
    steps = [s for s in steps if s >= 15000 and s % EVERY == 0]
    h = np.zeros(len(CEN))
    for s in steps:
        p, F = read(os.path.join(d, f'LDM_{s}stp.vtk'))
        air = np.linalg.norm(F['velocity'], axis=1) != 0
        h += np.histogram(p[air, 2], bins=EDGES)[0]*DT
    return h/NREL/np.diff(EDGES), len(steps)

def smooth(y, win=2.5):
    """격자 간격(2.5 m) 이동평균 — 입자가 2.5 m 격자 유동장으로 이류되면서 생기는
    격자 규모 계단을 지운다. 창이 격자 간격과 같아 2Δz 성분이 정확히 상쇄된다.
    경계에서는 창을 줄여 가장자리 왜곡을 막는다."""
    out = np.empty_like(y)
    for i, z in enumerate(CEN):
        m = np.abs(CEN - z) <= win/2
        out[i] = y[m].mean()
    return out

cache = os.path.join(HERE, 'profile_compare.npz')
if os.path.exists(cache):
    Z = np.load(cache); N, S = Z['N'], Z['S']
else:
    N, nN = profile('2026-09-15_leipzig-dt01-sponge-uy')
    S, nS = profile('2026-09-21_leipzig-most-stable')
    np.savez(cache, N=N, S=S, cen=CEN)
    print(f'스냅샷 중립 {nN} / 안정 {nS}')

fig, ax = newfig(st.W_ONE*1.10, st.H_ONE*1.25)
ax.axhspan(0, 2, color='#EFF2F4', zorder=0, lw=0)
ax.plot(smooth(N), CEN, '-', c=C['blue'],   lw=1.3, label='Neutral', zorder=3)
ax.plot(smooth(S), CEN, '-', c=C['orange'], lw=1.3, label=r'Stable ($T_\mathrm{s}$ = 285 K)', zorder=3)
ax.set_xlim(0, 5.2); ax.set_ylim(0, 50)
ax.set_xlabel(r'Residence time density   (s m$^{-1}$ per particle)')
ax.set_ylabel(r'Height  $z$   (m)')
i2 = CEN <= 2.0
# 교차 고도 — 여기서 위아래가 뒤집힌다
Ns, Ss = smooth(N), smooth(S)
xo = CEN[np.argmax((Ss < Ns) & (CEN > 10))]
ax.axhline(xo, c=C['grey'], lw=0.7, ls=(0, (4, 2.5)), zorder=2)
ax.annotate(f'crossover  {xo:.0f} m', xy=(5.05, xo-0.7), fontsize=7, color=C['grey'],
            ha='right', va='top')
# 캐노피 앵커 — 교차고가 지붕 높이 바로 아래라는 것이 이 그림의 요점이다
HR = np.median(np.load('/home/kdh/Urban-SPH/runs/lambdap/heightfield.npz')['H'][
               np.load('/home/kdh/Urban-SPH/runs/lambdap/heightfield.npz')['H'] > 1])
ax.axhline(HR, c=C['grey'], lw=0.6, ls=(0, (1.2, 2.0)), zorder=2)
ax.annotate(f'median roof  {HR:.0f} m', xy=(5.05, HR+0.7), fontsize=7, color=C['grey'],
            ha='right', va='bottom')
ax.annotate('pedestrian layer', xy=(0.14, 2.35), fontsize=7, color=C['grey'], ha='left', va='bottom')
legend_inside(ax, ncol=1, prefer='upper right')
finish(fig, ax, 'f38_vertical_profile_stable_vs_neutral', minor=False)
print('0-2 m  중립 %.3f s  안정 %.3f s  배율 %.2f' % (N[i2].sum()*DZ, S[i2].sum()*DZ, S[i2].sum()/N[i2].sum()))
for zlo, zhi in [(0,2),(2,5),(5,10),(10,20),(20,50)]:
    m = (CEN>zlo)&(CEN<=zhi)
    print(f'  {zlo:2d}-{zhi:2d} m : 중립 {N[m].mean():.4f}  안정 {S[m].mean():.4f}  ×{S[m].mean()/N[m].mean():.2f}')
