#!/usr/bin/env python3
"""소스별 침적 시계열 — LDM 덤프를 1 s 간격으로 훑어 침적 입자 수와 footprint 를 기록.

침적 판정은 dep_metrics.py 와 같다(속도 == 0). 출력: dep_timeseries.csv
  step, time, run, src, n_released, n_dep, area90_m2, r_med_m
"""
import csv, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A_DIR = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy')
B_DIR = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy-meanflow')
EVERY = 100          # step 간격 (1 s). LDM 덤프는 25 step 마다 있다.
CELL = 10.0

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']
S = {int(r['id']): (float(r['x']), float(r['y'])) for r in
     csv.DictReader(open(os.path.join(HERE, '..', 'sources13.csv')))}


def area90(xy, cell=CELL):
    if len(xy) < 5: return 0.0
    key = (np.floor(xy[:, 0]/cell).astype(np.int64)*1_000_000 + np.floor(xy[:, 1]/cell).astype(np.int64))
    _, c = np.unique(key, return_counts=True)
    c = np.sort(c)[::-1]; cum = np.cumsum(c)/c.sum()
    return float((np.searchsorted(cum, 0.90) + 1) * cell * cell)


def main():
    out = open(os.path.join(HERE, 'dep_timeseries.csv'), 'w', newline='')
    w = csv.writer(out); w.writerow(['step', 'time', 'run', 'src', 'n_released', 'n_dep', 'area90_m2', 'r_med_m'])
    for tag, d in (('A', A_DIR), ('B', B_DIR)):
        steps = sorted(int(re.search(r'_(\d+)stp', f).group(1))
                       for f in os.listdir(os.path.join(d, 'out')) if f.startswith('LDM_'))
        steps = [s for s in steps if s >= 15000 and s % EVERY == 0]
        for s in steps:
            p, F = read_vtk(os.path.join(d, 'out', f'LDM_{s}stp.vtk'))
            src = F['p_type'].astype(int) - 100
            dep = np.linalg.norm(F['velocity'], axis=1) == 0.0
            for i in sorted(S):
                m = src == i; md = m & dep
                xy = p[md][:, :2]
                r = np.linalg.norm(xy - np.array(S[i]), axis=1) if len(xy) else np.array([0.0])
                w.writerow([s, round(s*0.01, 2), tag, i, int(m.sum()), int(md.sum()),
                            round(area90(xy), 1), round(float(np.median(r)), 2)])
            out.flush()
        print(f'{tag}: {len(steps)} dumps 완료', flush=True)
    out.close()


if __name__ == '__main__':
    main()
