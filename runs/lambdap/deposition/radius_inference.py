#!/usr/bin/env python3
"""창 반경 R 선택을 통계적으로 정직하게 다루는 두 가지 방법.

문제: λp 는 창 반경 R 의 함수라 9개 반경 × 여러 지표를 훑으면 "가장 좋은 R" 이 늘 나온다.
n=10 에서 그 최대값을 그대로 보고하면 다중비교로 부풀려진 값이다.

방법 1 — **최대통계 치환검정**: 관측된 max_R |ρs(R)| 를, 결과를 무작위로 섞어 같은
          max_R 연산을 반복한 분포와 비교한다. 반경 탐색 자체가 귀무분포에 들어가므로
          "R 을 골랐다" 는 비판이 자동으로 해소된다.
방법 2 — **중첩 LOO**: 폴드 안에서 R 을 고른 모형과 R 을 150 으로 고정한 모형의
          leave-one-out 예측오차를 비교한다. 개선이 없으면 R 선택은 정보가 아니라 잡음이다.

  python3 radius_inference.py
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MORPH = os.path.normpath(os.path.join(HERE, '..'))
COURT = {3, 6, 10}
NPERM = 20000
rng = np.random.default_rng(20260920)

lam = {}
for r in csv.DictReader(open(os.path.join(MORPH, 'lambda_p_by_R.csv'))):
    lam[int(float(r['R_m']))] = {i: float(r[f'S{i}']) for i in range(1, 14)}
met = {}
for r in csv.DictReader(open(os.path.join(HERE, 'dep_metrics.csv'))):
    met.setdefault(r['run'], {})[int(r['src'])] = {k: float(v) for k, v in r.items() if k != 'run'}

KEEP = [i for i in sorted(met['A']) if i not in COURT]
RS = sorted(lam)
X = np.array([[lam[R][i] for i in KEEP] for R in RS])      # (nR, n)


def rank(a):
    return np.argsort(np.argsort(a)).astype(float)


def sp(a, b):
    return float(np.corrcoef(rank(a), rank(b))[0, 1])


METRICS = [('area90_m2', 'A90 (log)', True), ('f_r150', '150 m 갇힘', False),
           ('shift_m', '질량중심 이동', False), ('f_dep', '침적률', False)]

print(f'n = {len(KEEP)} (중정 제외),  반경 {len(RS)}개 {RS}')
print(f'치환 {NPERM:,}회\n')
print(f'{"지표":<14}{"R*":>5}{"ρs(R*)":>9}{"p 순진":>9}{"p 최대통계":>11}'
      f'{"ρs(150)":>10}{"p(150)":>9}')

for key, lab, uselog in METRICS:
    y0 = np.array([met['A'][i][key] for i in KEEP])
    y = np.log10(y0) if uselog else y0
    rs = np.array([sp(X[k], y) for k in range(len(RS))])
    kbest = int(np.argmax(np.abs(rs)))
    Tobs = abs(rs[kbest])
    k150 = RS.index(150)

    cnt_max = cnt_best = cnt_150 = 0
    for _ in range(NPERM):
        yp = rng.permutation(y)
        rp = np.array([sp(X[k], yp) for k in range(len(RS))])
        if np.abs(rp).max() >= Tobs - 1e-12: cnt_max += 1          # 반경 탐색 보정
        if abs(rp[kbest]) >= Tobs - 1e-12: cnt_best += 1           # 보정 없음
        if abs(rp[k150]) >= abs(rs[k150]) - 1e-12: cnt_150 += 1    # 사전 고정 R
    p_max = (cnt_max+1)/(NPERM+1); p_best = (cnt_best+1)/(NPERM+1); p150 = (cnt_150+1)/(NPERM+1)
    print(f'{lab:<14}{RS[kbest]:>5}{rs[kbest]:>9.2f}{p_best:>9.3f}{p_max:>11.3f}'
          f'{rs[k150]:>10.2f}{p150:>9.3f}')

# ---------------------------------------------------------------- 중첩 LOO
print('\n중첩 LOO — 폴드 안에서 R 을 고르면 예측이 나아지는가')
print(f'{"지표":<14}{"R 고정 150":>12}{"R 폴드내 선택":>14}{"개선":>8}{"선택된 R":>22}')
for key, lab, uselog in METRICS:
    y0 = np.array([met['A'][i][key] for i in KEEP])
    y = np.log10(y0) if uselog else y0
    n = len(y)
    err_fix, err_sel, picks = [], [], []
    for t in range(n):
        tr = np.arange(n) != t
        # (a) R = 150 고정
        b, a = np.polyfit(X[RS.index(150)][tr], y[tr], 1)
        err_fix.append(y[t] - (a + b*X[RS.index(150)][t]))
        # (b) 훈련 폴드에서만 |ρs| 최대인 R 선택
        k = int(np.argmax([abs(sp(X[kk][tr], y[tr])) for kk in range(len(RS))]))
        b, a = np.polyfit(X[k][tr], y[tr], 1)
        err_sel.append(y[t] - (a + b*X[k][t]))
        picks.append(RS[k])
    rf, rsl = np.sqrt(np.mean(np.square(err_fix))), np.sqrt(np.mean(np.square(err_sel)))
    print(f'{lab:<14}{rf:>12.4f}{rsl:>14.4f}{(rf-rsl)/rf*100:>7.1f}%   {sorted(set(picks))}')
