#!/usr/bin/env python3
"""시간 간격 수렴 — dt=0.005 런과 dt=0.01 런의 같은 물리시각 침적 비교.

dt=0.005: runs/2026-09-10_leipzig-full-bcfix_2 (정전으로 t=237.2 s 에서 중단, LDM 방출 t=150~)
dt=0.01 : runs/2026-09-15_leipzig-dt01-sponge-uy

주의 — 두 런은 (1) dt, (2) 천장 스펀지(없음 / uz+uy) 가 다르고 (3) 난류라 궤적이 카오스적이다.
따라서 입자 대 입자가 아니라 **소스별 통계량과 침적 지도의 통계적 일치**를 본다.
스펀지는 z>90 m 에서만 uz·uy 를 감쇠하므로 z<10 m 침적에는 직접 영향이 없다.

출력: dt_convergence.csv (소스별), 화면 요약
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
OLD = os.path.join(RUNS, '2026-09-10_leipzig-full-bcfix_2')   # dt=0.005
NEW = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy')  # dt=0.01
TIMES = [200.0, 210.0, 220.0, 230.0, 236.0]
CELL = 10.0

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(HERE, '..', 'sources13.csv')))}


def snap(run, t, dt):
    step = int(round(t/dt/25)*25)
    p, F = read_vtk(os.path.join(run, 'out', f'LDM_{step}stp.vtk'))
    src = F['p_type'].astype(int) - 100
    dep = np.linalg.norm(F['velocity'], axis=1) == 0.0
    return step, p, src, dep


def area90(xy, cell=CELL):
    if len(xy) < 5: return 0.0
    k = (np.floor(xy[:, 0]/cell).astype(np.int64)*1_000_000 + np.floor(xy[:, 1]/cell).astype(np.int64))
    _, c = np.unique(k, return_counts=True); c = np.sort(c)[::-1]
    return float((np.searchsorted(np.cumsum(c)/c.sum(), .90) + 1) * cell * cell)


def grid(xy, lim, cell=CELL):
    nx = int((lim[1]-lim[0])/cell)+1; ny = int((lim[3]-lim[2])/cell)+1
    g = np.zeros((nx, ny))
    ix = np.clip(((xy[:, 0]-lim[0])/cell).astype(int), 0, nx-1)
    iy = np.clip(((xy[:, 1]-lim[2])/cell).astype(int), 0, ny-1)
    np.add.at(g, (ix, iy), 1.0); return g


rows = []
print(f"{'t[s]':>6} {'src':>4} {'dt.005 n_dep':>12} {'dt.01 n_dep':>11} {'비':>5} {'f_dep .005':>10} {'.01':>7} "
      f"{'A90 .005':>9} {'.01':>8} {'지도 r':>6}")
for t in TIMES:
    s_o, p_o, src_o, dep_o = snap(OLD, t, 0.005)
    s_n, p_n, src_n, dep_n = snap(NEW, t, 0.01)
    lim = [min(p_o[dep_o][:, 0].min(), p_n[dep_n][:, 0].min()), max(p_o[dep_o][:, 0].max(), p_n[dep_n][:, 0].max()),
           min(p_o[dep_o][:, 1].min(), p_n[dep_n][:, 1].min()), max(p_o[dep_o][:, 1].max(), p_n[dep_n][:, 1].max())]
    for i in sorted(S):
        mo, mn = src_o == i, src_n == i
        do, dn = mo & dep_o, mn & dep_n
        fo = do.sum()/max(mo.sum(), 1); fn = dn.sum()/max(mn.sum(), 1)
        ao, an = area90(p_o[do][:, :2]), area90(p_n[dn][:, :2])
        go, gn = grid(p_o[do][:, :2], lim), grid(p_n[dn][:, :2], lim)
        b = (go > 0) | (gn > 0)
        r = float(np.corrcoef(go[b], gn[b])[0, 1]) if b.sum() > 2 else np.nan
        rows.append(dict(time=t, src=i, step_dt005=s_o, step_dt01=s_n, n_rel_dt005=int(mo.sum()), n_rel_dt01=int(mn.sum()),
                         n_dep_dt005=int(do.sum()), n_dep_dt01=int(dn.sum()), f_dep_dt005=fo, f_dep_dt01=fn,
                         area90_dt005=ao, area90_dt01=an, map_pearson=r))
        print(f"{t:6.0f} {i:4d} {int(do.sum()):12,} {int(dn.sum()):11,} {dn.sum()/max(do.sum(),1):5.2f} "
              f"{fo:10.3f} {fn:7.3f} {ao:9.0f} {an:8.0f} {r:6.2f}")
    # 전역
    go, gn = grid(p_o[dep_o][:, :2], lim), grid(p_n[dep_n][:, :2], lim)
    b = (go > 0) | (gn > 0)
    rows.append(dict(time=t, src=0, step_dt005=s_o, step_dt01=s_n, n_rel_dt005=int(len(p_o)), n_rel_dt01=int(len(p_n)),
                     n_dep_dt005=int(dep_o.sum()), n_dep_dt01=int(dep_n.sum()),
                     f_dep_dt005=dep_o.mean(), f_dep_dt01=dep_n.mean(), area90_dt005=area90(p_o[dep_o][:, :2]),
                     area90_dt01=area90(p_n[dep_n][:, :2]), map_pearson=float(np.corrcoef(go[b], gn[b])[0, 1])))
    print(f"{t:6.0f}  ALL {int(dep_o.sum()):12,} {int(dep_n.sum()):11,} {dep_n.sum()/dep_o.sum():5.2f} "
          f"{dep_o.mean():10.3f} {dep_n.mean():7.3f} {'':9} {'':8} {rows[-1]['map_pearson']:6.2f}", flush=True)

with open(os.path.join(HERE, 'dt_convergence.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\ndt_convergence.csv 저장')
