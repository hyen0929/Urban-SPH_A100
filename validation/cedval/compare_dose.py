#!/usr/bin/env python3
"""모델 vs 실험 평균 dosage (모양 비교 + 단일 상수 맞춤 후 FAC2/FB/NMSE).

단위 정의(T*, C*)가 미확정이라 모델 dosage 에 상수 하나를 곱해 맞춘다:
  상수 = 기하평균(실험/모델), 지면·연직 전체.  FAC2 는 이 상수에 무관하지 않으므로 '모양 기준'으로만 읽을 것.
사용: python3 -I validation/cedval/compare_dose.py runs/<run> <exp_dts_dir>
"""
import glob, re, sys, os
import numpy as np

RUN, EXP = sys.argv[1], sys.argv[2]


def exp_stats():
    res = {}
    for f in sorted(glob.glob(EXP + '/*.dts')):
        t = open(f, 'rb').read().decode('latin1')
        g = lambda k: float(re.search(r'([-\d.]+) \[' + re.escape(k), t).group(1))
        key = (g('measurement location X - full'), g('measurement location Y - full'),
               round(g('measurement location Z - full'), 2))
        for z in t.split('ZONE')[1:]:
            a = np.array([[float(v) for v in l.split()[:2]] for l in z.split('\n')[1:] if l.strip()])
            if len(a) < 3:
                continue
            dt = np.median(np.diff(a[:, 0]))
            res.setdefault(key, []).append(a[:, 1].sum() * dt)
    return res


d = np.load(RUN + '/analysis/receptor_ts.npz')
t, R, C = d['t'], d['rec'], d['C']
dt = np.median(np.diff(t))
D = C.sum(axis=2) * dt                      # [nr, npuff]
E = exp_stats()

rows = []
for r in range(len(R)):
    key = (R[r, 0], R[r, 1], round(R[r, 2], 2))
    if key not in E:
        continue
    e = np.array(E[key])
    rows.append((R[r], D[r].mean(), np.percentile(D[r], [25, 75]), e.mean(), np.percentile(e, [25, 75]),
                 D[r].std() / max(D[r].mean(), 1e-30), e.std() / e.mean()))

m = np.array([x[1] for x in rows]); o = np.array([x[3] for x in rows])
ok = m > 0
k = np.exp(np.mean(np.log(o[ok] / m[ok])))
p = m * k
fac2 = np.mean((p / o >= 0.5) & (p / o <= 2))
fb = 2 * (o.mean() - p.mean()) / (o.mean() + p.mean())
nmse = np.mean((o - p) ** 2) / (o.mean() * p.mean())
print(f'측정점 {len(rows)}  상수 k={k:.3e}  FAC2={fac2:.2f}  FB={fb:.2f}  NMSE={nmse:.2f}')
print('   x      y     z |   실험   모델(×k) 비 | CV실험 CV모델')
for (rec, mm, _, oo, _, cvm, cve), pp in zip(rows, p):
    print('%6.0f %6.0f %5.1f | %7.1f %7.1f %5.2f | %4.2f %4.2f' % (*rec, oo, pp, pp / oo, cve, cvm))
np.savez(RUN + '/analysis/compare_dose.npz', rec=np.array([x[0] for x in rows]), model=m, k=k, exp=o,
         model_iqr=np.array([x[2] for x in rows]), exp_iqr=np.array([x[4] for x in rows]),
         cv_model=np.array([x[5] for x in rows]), cv_exp=np.array([x[6] for x in rows]))
