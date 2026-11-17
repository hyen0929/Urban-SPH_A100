#!/usr/bin/env python3
"""50소스 × 4풍향, 벽면침적 ON — 벽면/지면 침적의 형태 의존성.
buffer_type: 0=공중, 1=지면, 2=벽면·지붕. 대조군은 같은 소스·유동장의 벽면 OFF 런.
"""
import csv, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
DEP=os.path.normpath(os.path.join(HERE,'..','deposition')); MORPH=os.path.normpath(os.path.join(HERE,'..'))
RUNS=os.path.normpath(os.path.join(HERE,'..','..'))
ns={'__file__':os.path.join(DEP,'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,sp=ns['read_vtk'],ns['spearman']
ANG=[0,90,180,270]
S=[(int(r['id']),float(r['x']),float(r['y']),float(r['lam75']),float(r['enclosure']))
   for r in csv.DictReader(open(os.path.join(MORPH,'sources_v2','sources_v3_50.csv')))]
ids=[s[0] for s in S]; X=np.array([s[3] for s in S])
R={}
for a in ANG:
    p,F=read(f'{RUNS}/2026-10-07_leipzig-src50w-rot{a}/out/LDM_30000stp.vtk')
    pt=F['p_type'].astype(int)-100; b=F['buffer_type'].astype(int)
    g=np.array([float(((pt==i)&(b==1)).sum()/max((pt==i).sum(),1)) for i in ids])
    w=np.array([float(((pt==i)&(b==2)).sum()/max((pt==i).sum(),1)) for i in ids])
    R[a]=(g,w)
print('풍향별 ρs(λp75, ·)   n=50')
print(f"{'각도':>5} {'지면 f_g':>9} {'벽면 f_w':>9} {'총 f_g+f_w':>11} {'벽/지면 비':>11}")
acc={k:[] for k in ('g','w','t','r')}
for a in ANG:
    g,w=R[a]; t=g+w; rr=w/np.maximum(g,1e-9)
    v=[sp(X,g),sp(X,w),sp(X,t),sp(X,rr)]
    for k,q in zip(('g','w','t','r'),v): acc[k].append(q)
    print(f'{a:4d}° {v[0]:+9.2f} {v[1]:+9.2f} {v[2]:+11.2f} {v[3]:+11.2f}')
print(f"{'평균':>5} " + ' '.join(f'{np.mean(acc[k]):+9.2f}' if k!='t' else f'{np.mean(acc[k]):+11.2f}' for k in ('g','w','t','r')))
print(f"{'폭':>5} " + ' '.join(f'{max(acc[k])-min(acc[k]):9.2f}' if k!='t' else f'{max(acc[k])-min(acc[k]):11.2f}' for k in ('g','w','t','r')))

print('\n벽면 포함이 λp 상관을 바꾸는가 (지면만 → 지면+벽면)')
for a in ANG:
    g,w=R[a]
    print(f'  {a:4d}°  {sp(X,g):+.2f} → {sp(X,g+w):+.2f}   (Δ {sp(X,g+w)-sp(X,g):+.2f})')

print('\n벽/지면 비의 분포 (소스 50개)')
for a in ANG:
    g,w=R[a]; rr=w/np.maximum(g,1e-9)
    print(f'  {a:4d}°  중앙값 {np.median(rr):.3f}  범위 {rr.min():.3f}~{rr.max():.3f}')
rows=[]
for a in ANG:
    g,w=R[a]
    for k,i in enumerate(ids):
        rows.append(dict(deg=a,src=i,lam75=X[k],f_ground=g[k],f_wall=w[k],enc=S[k][4]))
with open(os.path.join(HERE,'wall50_metrics.csv'),'w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print('\n저장: wall50_metrics.csv')
