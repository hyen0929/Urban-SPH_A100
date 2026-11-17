#!/usr/bin/env python3
"""29 소스 × 4 풍향 — 형태-확산 상관이 풍향을 넘어 살아남는가.

지표 정의는 deposition/dep_metrics.py 재사용. 소스는 도시와 함께 회전한다.
주 예측변수는 λp(75): 소스 최소간격 120 m 이라 R=75 창은 10% 만 겹친다.
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
CX,CY=657.5,0.0; ANG=[0,90,180,270]; KEYS=['f_dep','area90_m2','f_r150','shift_m','elong']

S=[]
for r in csv.DictReader(open(os.path.join(MORPH,'sources_v2','sources_v2.csv'))):
    S.append(dict(id=int(r['id']), x=float(r['x']), y=float(r['y']),
                  lam75=float(r['lam75']), lam150=float(r['lam150']), enc=float(r['enclosure'])))
ids=[s['id'] for s in S]
court=np.array([s['enc']>0.8 for s in S]); keep=~court
X75=np.array([s['lam75'] for s in S])
def rot(x,y,deg):
    t=np.radians(deg);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y-CY
    return CX+dx*c-dy*s, CY+dx*s+dy*c

M={}; rows=[]
for deg in ANG:
    d=f'{RUNS}/2026-09-30_leipzig-src29-rot{deg}'
    p,src,dep=load(d)                      # step 30000 (t=300 s)
    per={}
    for k,s in enumerate(S):
        pid=s['id']; m=src==pid
        f=footprint(p[m&dep][:,:2], np.array(rot(s['x'],s['y'],deg)))
        f['f_dep']=float((m&dep).sum()/max(m.sum(),1))
        per[pid]=f; rows.append(dict(deg=deg,src=pid,lam75=s['lam75'],court=int(court[k]),
                                     **{q:f[q] for q in KEYS}))
    M[deg]=per

def boot(x,y,n=4000,seed=0):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(np.unique(x[i]))<4: continue
        out.append(sp(x[i],y[i]))
    return np.percentile(out,[2.5,97.5])

print('ρs(λp75, 지표)  — 중정 제외 n=%d,  [ ] 는 부트스트랩 95%% 신뢰구간' % keep.sum())
print(f"{'각도':>5} | " + '  '.join(f'{k:^22}' for k in KEYS))
tab={k:[] for k in KEYS}
for deg in ANG:
    cells=[]
    for k in KEYS:
        y=np.array([M[deg][i][k] for i in ids])
        r=sp(X75[keep],y[keep]); lo,hi=boot(X75[keep],y[keep])
        tab[k].append(r); cells.append(f'{r:+.2f} [{lo:+.2f},{hi:+.2f}]')
    print(f'{deg:4d}° | ' + '  '.join(f'{c:^22}' for c in cells))
print(f"{'폭':>5} | " + '  '.join(f'{max(tab[k])-min(tab[k]):^22.2f}' for k in KEYS))
print(f"{'평균':>5} | " + '  '.join(f'{np.mean(tab[k]):^22.2f}' for k in KEYS))

print('\n[대조] 기존 13 소스만으로 같은 계산')
old13=[s['id'] for s in S if s['id']<=13]
k13=np.array([s['enc']<=0.8 for s in S if s['id']<=13])
x13=np.array([s['lam75'] for s in S if s['id']<=13])
for k in ['area90_m2','f_dep']:
    v=[sp(x13[k13], np.array([M[deg][i][k] for i in old13])[k13]) for deg in ANG]
    print(f'  {k:10s} ' + '  '.join(f'{deg}° {x:+.2f}' for deg,x in zip(ANG,v)) + f'   폭 {max(v)-min(v):.2f}')
with open(os.path.join(HERE,'rot29_metrics.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\n저장: rot29_metrics.csv')
