#!/usr/bin/env python3
"""소스 입자수 수렴성 — 지표가 N 에 대해 언제 안정되는가.

다섯 런은 같은 동결 평균장 위에서 입자 개수만 다르다. 지표 정의는
deposition/dep_metrics.py 를 그대로 쓴다. 기준은 최대 N(200,000).
"""
import csv, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
DEP=os.path.normpath(os.path.join(HERE,'..','deposition'))
MORPH=os.path.normpath(os.path.join(HERE,'..'))
RUNS=os.path.normpath(os.path.join(HERE,'..','..'))
ns={'__file__':os.path.join(DEP,'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
load,footprint,sp=ns['load'],ns['footprint'],ns['spearman']

NS=[10000,25000,50000,100000,200000]
COURT={3,6,10}
S={int(r['id']):(float(r['x']),float(r['y'])) for r in csv.DictReader(open(f'{MORPH}/sources13.csv'))}
lam={int(float(r['R_m'])):{i:float(r[f'S{i}']) for i in range(1,14)}
     for r in csv.DictReader(open(f'{MORPH}/lambda_p_by_R.csv'))}
ids=sorted(S); keep=np.array([i not in COURT for i in ids])
x150=np.array([lam[150][i] for i in ids]); x75=np.array([lam[75][i] for i in ids])
KEYS=['f_dep','area90_m2','f_r150','shift_m','elong']

M={}; rows=[]
for N in NS:
    d=f'{RUNS}/2026-09-28_leipzig-pconv-N{N}'
    p,src,dep=load(d)
    per={}
    for i in ids:
        m=src==i
        f=footprint(p[m&dep][:,:2], np.array(S[i]))
        f['f_dep']=float((m&dep).sum()/max(m.sum(),1))
        per[i]=f
        rows.append(dict(N=N,src=i,**{k:f[k] for k in KEYS}))
    M[N]=dict(per=per, dep_tot=float(dep.mean()), n=len(p))

print('전체 침적률 과 λp 상관 (중정 제외 n=10)')
print(f"{'N/소스':>8} {'총입자':>10} {'침적%':>7} | " + ' '.join(f'{k:>10}' for k in KEYS))
for N in NS:
    c=[sp(x150[keep], np.array([M[N]['per'][i][k] for i in ids])[keep]) for k in KEYS]
    print(f'{N:8,} {M[N]["n"]:10,} {100*M[N]["dep_tot"]:6.2f}  | ' + ' '.join(f'{v:+10.2f}' for v in c))

REF=NS[-1]
print(f'\n소스별 지표의 N={REF:,} 대비 상대오차 [%] (13개 소스 중앙값 / 최대)')
print(f"{'N/소스':>8} | " + ' '.join(f'{k:>16}' for k in KEYS))
for N in NS[:-1]:
    out=[]
    for k in KEYS:
        a=np.array([M[N]['per'][i][k] for i in ids]); b=np.array([M[REF]['per'][i][k] for i in ids])
        e=100*np.abs(a-b)/np.abs(b)
        out.append(f'{np.median(e):6.1f} /{np.max(e):6.1f}')
    print(f'{N:8,} | ' + ' '.join(f'{o:>16}' for o in out))

with open(os.path.join(HERE,'conv_metrics.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\n저장: conv_metrics.csv')
