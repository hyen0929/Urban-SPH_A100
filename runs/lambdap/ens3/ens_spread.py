#!/usr/bin/env python3
"""난류 실현 변동의 크기 — 같은 소스를 세 시간창에 방출해 그 흩어짐을 잰다.

창 w 는 t0_w 부터 150 s 동안 방출한다(t0 = 150/300/450 s). 같은 '나이'(t - t0)에서
비교해야 공정하다. LDM 덤프가 t=200~700 s 뿐이라 공통 최대 나이는 250 s.
"""
import csv, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
DEP=os.path.normpath(os.path.join(HERE,'..','deposition')); MORPH=os.path.normpath(os.path.join(HERE,'..'))
D='/home/kdh/Urban-SPH/runs/2026-09-30_leipzig-ens3-src29'
ns={'__file__':os.path.join(DEP,'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,footprint,sp=ns['read_vtk'],ns['footprint'],ns['spearman']

S=[dict(id=int(r['id']),x=float(r['x']),y=float(r['y']),lam75=float(r['lam75']),enc=float(r['enclosure']))
   for r in csv.DictReader(open(os.path.join(MORPH,'sources_v2','sources_v2.csv')))]
NS=len(S); ids=[s['id'] for s in S]
court=np.array([s['enc']>0.8 for s in S]); keep=~court
X=np.array([s['lam75'] for s in S])
T0={1:150.0,2:300.0,3:450.0}
KEYS=['f_dep','area90_m2','f_r150','shift_m']
AGES=[150,250]

res={}
for age in AGES:
    for w in (1,2,3):
        step=int((T0[w]+age)*100)
        f=f'{D}/out/LDM_{step}stp.vtk'
        if not os.path.exists(f): print(f'  없음 {f}'); continue
        p,F=read(f); pt=F['p_type'].astype(int); v=np.linalg.norm(F['velocity'],axis=1)
        lo=101+(w-1)*NS
        for k,s in enumerate(S):
            m=pt==lo+k; dep=m&(v==0)
            if m.sum()==0: continue
            fp=footprint(p[dep][:,:2], np.array([s['x'],s['y']]))
            fp['f_dep']=float(dep.sum()/m.sum())
            res[(age,w,s['id'])]=fp

print('=== 창별 λp(75) 상관 — 같은 나이에서 ===')
print(f"{'나이':>5} {'창':>3} | " + '  '.join(f'{k:>10}' for k in KEYS))
spread={}
for age in AGES:
    vals={k:[] for k in KEYS}
    for w in (1,2,3):
        row=[]
        for k in KEYS:
            y=np.array([res[(age,w,i)][k] for i in ids]); r=sp(X[keep],y[keep])
            vals[k].append(r); row.append(r)
        print(f'{age:4d}s {w:3d} | ' + '  '.join(f'{r:+10.2f}' for r in row))
    print(f'{"":>5} {"폭":>3} | ' + '  '.join(f'{max(vals[k])-min(vals[k]):10.2f}' for k in KEYS))
    spread[age]=vals

print('\n=== 소스별 지표의 창 간 변동 (실현 잡음) ===')
print(f"{'나이':>5} | " + '  '.join(f'{k:>22}' for k in KEYS) + '   (상대표준편차 중앙값 / 최대)')
for age in AGES:
    out=[]
    for k in KEYS:
        A=np.array([[res[(age,w,i)][k] for w in (1,2,3)] for i in ids])
        cv=100*A.std(1,ddof=1)/np.abs(A.mean(1))
        out.append(f'{np.median(cv):6.1f}% / {np.max(cv):6.1f}%')
    print(f'{age:4d}s | ' + '  '.join(f'{o:>22}' for o in out))

print('\n=== 핵심 비교: 실현 변동 vs 풍향 변동 (ρs 폭) ===')
print(f"{'지표':>12} {'창 간(실현)':>12} {'풍향 간':>10}")
WIND={'f_dep':0.40,'area90_m2':0.48,'f_r150':0.56,'shift_m':None}
for k in KEYS:
    s=max(spread[250][k])-min(spread[250][k])
    w=WIND.get(k)
    print(f'{k:>12} {s:12.2f} {("%.2f"%w) if w else "-":>10}')
import json; json.dump({f'{a}_{w}_{i}':{k:res[(a,w,i)][k] for k in KEYS}
    for (a,w,i) in res}, open(os.path.join(HERE,'ens_spread.json'),'w'), indent=1)
print('\n저장: ens_spread.json')
