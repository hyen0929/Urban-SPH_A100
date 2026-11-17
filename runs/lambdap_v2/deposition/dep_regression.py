#!/usr/bin/env python3
"""침적률의 λp 잔차를 2번째 형태지표가 설명하는가 (handoff §4.2).

λp 하나로는 침적률 설명력이 약하므로(Spearman −0.6), λp 와 상관이 낮은 지표
(H_built r=0.14, λf/λp 바람정렬)를 두 번째 예측자로 넣어 본다.

  - 단순회귀  y ~ λp
  - 중회귀    y ~ λp + z   (z = H_built, lf_ratio=λf/λp, encl_50, W_char)
  - 편상관    corr(resid(y~λp), resid(z~λp))
  - LOO 교차검증 RMSE 로 과적합 확인 (10점에 2예측자)

출력: dep_regression.csv, 화면 요약
"""
import csv, itertools, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}


def load(R=150):
    morph = {}
    for r in csv.DictReader(open(os.path.join(MORPH, 'morph_sources.csv'))):
        if abs(float(r['R']) - R) < 1e-6:
            morph[int(r['id'])] = {k: float(v) for k, v in r.items() if k not in ('label',)}
    met = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
        met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}
    return morph, met


def fit(X, y):
    """최소제곱 + R², adj R², LOO RMSE"""
    A = np.column_stack([np.ones(len(y))] + [X[:, j] for j in range(X.shape[1])])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    ss = ((y - pred)**2).sum(); st = ((y - y.mean())**2).sum()
    r2 = 1 - ss/st
    n, p = len(y), X.shape[1]
    adj = 1 - (1-r2)*(n-1)/(n-p-1)
    loo = []
    for i in range(n):
        k = np.arange(n) != i
        b, *_ = np.linalg.lstsq(A[k], y[k], rcond=None)
        loo.append((A[i] @ b - y[i])**2)
    return beta, r2, adj, float(np.sqrt(np.mean(loo)))


def partial(y, x, z):
    ry = y - np.polyval(np.polyfit(x, y, 1), x)
    rz = z - np.polyval(np.polyfit(x, z, 1), x)
    return float(np.corrcoef(ry, rz)[0, 1])


def main():
    out = []
    for R in (75, 150, 250):
        morph, met = load(R)
        ids = [i for i in sorted(morph) if i not in COURT]
        lam = np.array([morph[i]['lambda_p'] for i in ids])
        extras = {'H_built': np.array([morph[i]['H_built'] for i in ids]),
                  'lf_ratio': np.array([morph[i]['lambda_f']/morph[i]['lambda_p'] for i in ids]),
                  'lambda_f': np.array([morph[i]['lambda_f'] for i in ids]),
                  'encl_50': np.array([morph[i]['encl_50'] for i in ids]),
                  'W_char': np.array([morph[i]['W_char'] for i in ids])}
        for run in ('A', 'B'):
            for tgt in ('f_dep', 'area90_m2'):
                y = np.array([met[run][i][tgt] for i in ids])
                if tgt == 'area90_m2': y = np.log10(y)
                _, r2_1, adj_1, loo_1 = fit(lam[:, None], y)
                out.append(dict(R=R, run=run, target=tgt, model='lambda_p', r2=r2_1, adj_r2=adj_1, loo_rmse=loo_1,
                                partial_r=float('nan'), corr_with_lambda=1.0))
                for nm, z in extras.items():
                    _, r2_2, adj_2, loo_2 = fit(np.column_stack([lam, z]), y)
                    out.append(dict(R=R, run=run, target=tgt, model=f'lambda_p+{nm}', r2=r2_2, adj_r2=adj_2,
                                    loo_rmse=loo_2, partial_r=partial(y, lam, z),
                                    corr_with_lambda=float(np.corrcoef(lam, z)[0, 1])))
    with open(os.path.join(HERE, 'dep_regression.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print(f"{'R':>4} {'run':>3} {'target':>10} {'model':>22} {'R²':>6} {'adj':>6} {'LOO':>8} {'편상관':>7} {'z~λp':>6}")
    for r in out:
        if r['run'] != 'A': continue
        print(f"{r['R']:>4} {r['run']:>3} {r['target']:>10} {r['model']:>22} {r['r2']:6.3f} {r['adj_r2']:6.3f} "
              f"{r['loo_rmse']:8.4f} {r['partial_r']:+7.2f} {r['corr_with_lambda']:+6.2f}")


if __name__ == '__main__':
    main()
