#!/usr/bin/env python3
"""침적량을 형태가 아니라 **국소 유동**으로 설명할 수 있는가 (REPORT §8 의 후속).

λp 로는 소스별 침적률이 설명되지 않았다(Spearman −0.6, 2번째 형태지표를 넣어도 LOO 악화).
평균장(mean_field.bin, t=150~299 s)에서 각 소스 주변의 유동량을 뽑아 다시 건다.

또한 대수법칙 적합의 d–z0 트레이드오프를 RMSE 지형으로 확인한다(canopy_aero.py 보완).

출력: flow_predictors.csv, 화면 요약
"""
import csv, json, os, struct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.normpath(os.path.join(HERE, '..', '..'))
A = os.path.join(RUNS, '2026-09-15_leipzig-dt01-sponge-uy')
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}
KAPPA = 0.4

ns = {'__file__': os.path.join(HERE, 'dep_metrics.py')}
exec(open(os.path.join(HERE, 'dep_metrics.py')).read().split('def load(')[0], ns)
read_vtk = ns['read_vtk']


def spearman(a, b):
    return float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1])


def loo_rmse(X, y):
    A_ = np.column_stack([np.ones(len(y))] + [X[:, j] for j in range(X.shape[1])])
    e = []
    for i in range(len(y)):
        k = np.arange(len(y)) != i
        b, *_ = np.linalg.lstsq(A_[k], y[k], rcond=None)
        e.append((A_[i] @ b - y[i])**2)
    b, *_ = np.linalg.lstsq(A_, y, rcond=None)
    r2 = 1 - ((y - A_ @ b)**2).sum()/((y - y.mean())**2).sum()
    return float(np.sqrt(np.mean(e))), float(r2)


def main():
    with open(os.path.join(A, 'out', 'mean_field.bin'), 'rb') as f:
        hdr = f.read(48); n, ac = struct.unpack_from('<Ii', hdr, 16)
        rec = np.fromfile(f, dtype=np.dtype([('ux', '<f4'), ('uy', '<f4'), ('uz', '<f4'),
                                             ('k', '<f4'), ('e', '<f4')]), count=n)
    p, F = read_vtk(os.path.join(A, 'out', 'fluid_30000stp.vtk'))
    fluid = (F['p_type'] == 1) & (rec['ux'] != 0)
    hf = np.load(os.path.join(MORPH, 'heightfield.npz'))
    H, x0, y0, cell = hf['H'], float(hf['x0']), float(hf['y0']), float(hf['cell'])
    ix = np.floor((p[:, 0]-x0)/cell).astype(int); iy = np.floor((p[:, 1]-y0)/cell).astype(int)
    ok = (iy >= 0) & (iy < H.shape[0]) & (ix >= 0) & (ix < H.shape[1])
    hz = np.zeros(len(p)); hz[ok] = H[iy[ok], ix[ok]]
    use = fluid & (p[:, 2] > hz + 0.1)

    S = {int(r['id']): (float(r['x']), float(r['y']), float(r['z'])) for r in
         csv.DictReader(open(os.path.join(MORPH, 'sources13.csv')))}
    met = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
        met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
    lam = {}
    for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
        lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}

    rows = []
    for i, (sx, sy, sz) in S.items():
        rr = np.hypot(p[:, 0]-sx, p[:, 1]-sy)
        rec_i = {}
        for R in (25, 50, 100):
            m = use & (rr <= R) & (p[:, 2] <= 10.0)
            if m.sum() < 50: continue
            U = np.hypot(rec['ux'][m], rec['uy'][m])
            rec_i.update({f'U{R}': float(U.mean()), f'k{R}': float(rec['k'][m].mean()),
                          f'w{R}': float(rec['uz'][m].mean()), f'wabs{R}': float(np.abs(rec['uz'][m]).mean()),
                          f'n{R}': int(m.sum())})
        # 캐노피 높이(주변 평균 건물높이) 에서의 연직속도 — 환기 프록시
        mm = use & (rr <= 50)
        hloc = float(H[(np.abs(np.arange(H.shape[0])*cell + y0 - sy) <= 50)][:, :][0].mean()) if False else 0.0
        zc = 18.8
        mcan = use & (rr <= 50) & (np.abs(p[:, 2] - zc) <= 2.5)
        rec_i['w_canopy'] = float(rec['uz'][mcan].mean()) if mcan.sum() > 50 else np.nan
        rec_i['k_canopy'] = float(rec['k'][mcan].mean()) if mcan.sum() > 50 else np.nan
        rows.append(dict(src=i, lambda_p75=lam[75][i], f_dep=met['A'][i]['f_dep'],
                         area90=met['A'][i]['area90_m2'], f_r150=met['A'][i]['f_r150'], **rec_i))
    with open(os.path.join(HERE, 'flow_predictors.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    ids = [r['src'] for r in rows]
    keep = np.array([i not in COURT for i in ids])
    y = np.array([r['f_dep'] for r in rows])
    print(f"{'예측자':>12} {'Spearman(13)':>13} {'Spearman(10)':>13} {'R²(10)':>8} {'LOO(10)':>9}")
    cand = ['lambda_p75', 'U25', 'U50', 'U100', 'k25', 'k50', 'k100', 'wabs25', 'wabs50', 'w_canopy', 'k_canopy']
    base = None
    for c in cand:
        x = np.array([r.get(c, np.nan) for r in rows])
        if np.isnan(x).any(): continue
        loo, r2 = loo_rmse(x[keep][:, None], y[keep])
        if c == 'lambda_p75': base = loo
        print(f"{c:>12} {spearman(x, y):13.2f} {spearman(x[keep], y[keep]):13.2f} {r2:8.3f} {loo:9.4f}")
    print(f"\n2변수 (중정 제외 10점, 기준 λp 단독 LOO={base:.4f})")
    lamv = np.array([r['lambda_p75'] for r in rows])[keep]
    for c in cand[1:]:
        x = np.array([r.get(c, np.nan) for r in rows])
        if np.isnan(x).any(): continue
        loo, r2 = loo_rmse(np.column_stack([lamv, x[keep]]), y[keep])
        flag = '  ← 개선' if loo < base else ''
        print(f"  λp + {c:<10} R²={r2:.3f}  LOO={loo:.4f}{flag}")

    # --- 대수법칙 d–z0 트레이드오프 ---
    prof = list(csv.DictReader(open(os.path.join(HERE, 'canopy_profile.csv'))))
    zc_ = np.array([float(r['z_m']) for r in prof]); U_ = np.array([float(r['U_mean']) for r in prof])
    Hb = json.load(open(os.path.join(HERE, 'canopy_aero.json')))['morphology']['H_built']
    res = []
    for d in np.arange(0, 0.95*Hb, 0.25):
        m = (zc_ > max(2*Hb, d+2)) & (zc_ < 8*Hb)
        X = np.log(zc_[m] - d); b, a = np.polyfit(X, U_[m], 1)
        rmse = float(np.sqrt(((U_[m] - (b*X+a))**2).mean()))
        res.append((d, rmse, b*KAPPA, float(np.exp(-a/b))))
    res = np.array(res); best = res[res[:, 1].argmin()]
    ok5 = res[res[:, 1] < 1.05*best[1]]
    print(f"\n[대수법칙 적합 강건성] 최적 d={best[0]:.1f} m, z0={best[3]:.3f} m, u*={best[2]:.3f}")
    print(f"  RMSE 가 최적의 5% 이내인 범위: d {ok5[:,0].min():.1f}~{ok5[:,0].max():.1f} m, "
          f"z0 {ok5[:,3].min():.3f}~{ok5[:,3].max():.3f} m, u* {ok5[:,2].min():.3f}~{ok5[:,2].max():.3f}")


if __name__ == '__main__':
    main()
