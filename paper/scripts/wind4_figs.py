#!/usr/bin/env python3
"""50 방출점 × 4 풍향에서 나오는 두 가지 주장의 근거 그림.

  f46  건물이 빽빽할수록 오염물이 땅보다 벽·지붕에 붙는다. 네 풍향이 같은 관계를 그린다.
       → 형태로 설명되는 부분은 풍향에 흔들리지 않는다.
  f47  같은 방출점이라도 침적량은 풍향에 따라 2배 안팎으로 바뀐다.
       → 한 풍향 결과로 일반화할 수 없다.

f46 은 벽면 침적을 켠 계열(DD_WALL_DEP=1), f47 은 끈 계열이다. 둘 다 50소스 × 4풍향 평균장.
"""
import csv, os, sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '/home/kdh/Urban-SPH/runs/lambdap/deposition')
import figstyle as st
from figstyle import C, newfig, finish

HERE = os.path.dirname(os.path.abspath(__file__))
st.OUT = os.environ.get('FIG_OUT', os.path.join(HERE, '..', 'figs', 'main'))
DEGC = [('0', C['blue'], 'o'), ('90', C['orange'], 's'),
        ('180', C['green'], '^'), ('270', C['purple'], 'D')]

def sp(a, b):
    ra = np.argsort(np.argsort(a)).astype(float); rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])

# ── f46  침적물이 어디로 가는가 — 땅 대 벽·지붕 ─────────────────────────
#  방출점에 따라 제염 대상이 갈린다는 것을 보이는 그림. 벽면 침적을 켠 계열
#  (DD_WALL_DEP=1), 50소스 × 4풍향 평균. 벽면 모델이 단순하므로(벽에서 2 m 안의
#  입자를 후보로 삼고 재부유 없음) 절대값보다 방출점 간 상대 비교로 읽는다.
W = defaultdict(dict)
for r in csv.DictReader(open(os.path.join(HERE, '..', 'data/02_wall/wall4_n50_meanflow_metrics.csv'))):
    W[int(r['src'])][r['deg']] = (float(r['f_ground']), float(r['f_wall']))
SRC = sorted(W)
share = np.array([np.mean([w/(g+w) for g, w in W[s].values()]) for s in SRC])*100
o = np.argsort(share)
x = np.arange(1, len(SRC)+1)
fig, ax = newfig()
ax.bar(x, 100-share[o], width=0.78, color=C['light'], ec='none', zorder=2,
       label='ground')
ax.bar(x, share[o], bottom=100-share[o], width=0.78, color=C['orange'], ec='none',
       zorder=2, label='walls and roofs')
ax.set_xlabel('release point, ranked by share on walls and roofs', fontsize=6.5)
ax.set_ylabel('where the deposited mass went  (%)', fontsize=6.5)
ax.tick_params(labelsize=6)
ax.set_xlim(0, 51); ax.set_ylim(0, 100)
ax.set_xticks(np.arange(0, 51, 10)); ax.set_yticks(np.arange(0, 101, 25))
ax.legend(loc='lower left', bbox_to_anchor=(0.015, 0.03), ncol=2, fontsize=6,
          handlelength=1.1, handletextpad=0.35, borderpad=0.35, columnspacing=0.9)
finish(fig, ax, 'f46_deposition_path_by_source', minor_y=False)
print(f"f46  벽·지붕 비중: 최소 {share.min():.1f} %  중앙값 {np.median(share):.1f} %  최대 {share.max():.1f} %")

# ── f47  같은 방출점, 풍향에 따른 변동 ────────────────────────────────────
rows = [r for r in csv.DictReader(open(os.path.join(HERE, '..', 'data/04_nsources/per_source_metrics.csv')))
        if r['n_sources'] == '50']
D = {d: {int(r['src']): float(r['f_dep'])*100 for r in rows if r['deg'] == d} for d, _, _ in DEGC}
S = sorted(D['0'])
M = np.array([[D[d][s] for d, _, _ in DEGC] for s in S])          # (50, 4)
o = np.argsort(M.mean(1))
x = np.arange(1, len(S)+1)
fig, ax = newfig()
ax.vlines(x, M[o].min(1), M[o].max(1), color=C['light'], lw=1.6, zorder=1)
for k, (d, col, mk) in enumerate(DEGC):
    ax.plot(x, M[o, k], mk, color=col, mec='white', mew=0.35, ms=2.8, ls='none',
            label=f'wind {d}°', zorder=3)
ax.set_xlabel('release point, ranked by 4-wind mean', fontsize=6.5)
ax.set_ylabel('deposited fraction  (%)', fontsize=6.5)
ax.tick_params(labelsize=6)
ax.set_xlim(0, 51); ax.set_ylim(0, 65)
ax.set_xticks(np.arange(0, 51, 10))
ax.legend(loc='upper left', bbox_to_anchor=(0.015, 0.985), ncol=2, fontsize=6,
          handlelength=1.1, handletextpad=0.35, labelspacing=0.3, columnspacing=0.9, borderpad=0.35)
finish(fig, ax, 'f47_per_source_wind_spread')
rng = M.max(1)/M.min(1)
print(f"f47  방출점별 4풍향 최대/최소 비: 중앙값 {np.median(rng):.2f}배 (범위 {rng.min():.2f}~{rng.max():.2f})")
print(f"     방출점 간 (4풍향 평균) 최대/최소 {M.mean(1).max()/M.mean(1).min():.2f}배,"
      f"  단일 풍향 최대 {max((M[:,k].max()/M[:,k].min()) for k in range(4)):.2f}배")
