#!/usr/bin/env python3
"""소스 13 / 29 / 50 × 풍향 4방향 — 소스를 늘리면 검정력이 실제로 오르는가.
전부 평균장(Run B), LDM 1 s 출력, t=300 s. λp(75) 주 예측변수.
"""
import csv, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
DEP=os.path.normpath(os.path.join(HERE,'..','deposition')); MORPH=os.path.normpath(os.path.join(HERE,'..'))
RUNS=os.path.normpath(os.path.join(HERE,'..','..'))
ns={'__file__':os.path.join(DEP,'dep_metrics.py')}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
load,footprint,sp=ns['load'],ns['footprint'],ns['spearman']
CX=657.5; ANG=[0,90,180,270]; KEYS=['f_dep','area90_m2','f_r150']

# λp 장 (회전 불변) — 소스 좌표에서 직접 읽는다
hf=np.load(os.path.join(MORPH,'heightfield.npz'))
H,cell,x0,y0=hf['H'],float(hf['cell']),float(hf['x0']),float(hf['y0'])
B=(H>1.0).astype(np.float64); ny,nx=H.shape
def lamfield(R):
    r=int(round(R/cell)); yy,xx=np.mgrid[-r:r+1,-r:r+1]
    K=((xx*xx+yy*yy)*cell*cell<=R*R).astype(np.float64); K/=K.sum()
    Kp=np.zeros_like(B); Kp[:K.shape[0],:K.shape[1]]=K
    Kp=np.roll(np.roll(Kp,-(K.shape[0]//2),0),-(K.shape[1]//2),1)
    return np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp),B.shape)
L75=lamfield(75.0)
def lam_at(x,y):
    j=int(np.floor((x-x0)/cell)); k=int(np.floor((y-y0)/cell))
    return float(L75[min(max(k,0),ny-1), min(max(j,0),nx-1)])
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c

SETS=[('13', '../sources13.csv',          {0:'2026-09-30_leipzig-src29b-rot0'}),  # 13 은 29 런의 부분집합
      ('29', '../sources_v2/sources_v2_29.csv', None),
      ('50', '../sources_v2/sources_v3_50.csv', None)]
def runs_for(tag,deg):
    pre='src29b' if tag in ('13','29') else 'src50'
    return f'{RUNS}/2026-10-06_leipzig-{pre}-rot{deg}'

def boot(x,y,n=4000,seed=0):
    rng=np.random.default_rng(seed); o=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(np.unique(x[i]))<4: continue
        o.append(sp(x[i],y[i]))
    return np.percentile(o,[2.5,97.5])

print(f"{'n':>4} {'지표':>10} | " + '  '.join(f'{d:>6}°' for d in ANG) + f"  {'폭':>6} {'평균':>6} {'CI폭(평균)':>10}")
store={}
for tag,f,_ in SETS:
    S=[(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(os.path.join(HERE,f)))]
    if tag=='13': S=[s for s in S if s[0]<=13]
    X=np.array([lam_at(x,y) for _,x,y in S])
    per={}
    for deg in ANG:
        p,src,dep=load(runs_for(tag,deg))
        M={}
        for sid,x,y in S:
            m=src==sid
            fp=footprint(p[m&dep][:,:2], np.array(rot(x,y,deg)))
            fp['f_dep']=float((m&dep).sum()/max(m.sum(),1))
            M[sid]=fp
        per[deg]=M
    for k in KEYS:
        v=[];ciw=[]
        for deg in ANG:
            y=np.array([per[deg][sid][k] for sid,_,_ in S])
            v.append(sp(X,y)); lo,hi=boot(X,y); ciw.append(hi-lo)
        store[(tag,k)]=v
        print(f'{tag:>4} {k:>10} | ' + '  '.join(f'{q:+7.2f}' for q in v) +
              f'  {max(v)-min(v):6.2f} {np.mean(v):+6.2f} {np.mean(ciw):10.2f}')

# 풀링 순열검정
print('\n4방향 풀링 순열검정 (20,000회)')
for tag,f,_ in SETS:
    S=[(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(os.path.join(HERE,f)))]
    if tag=='13': S=[s for s in S if s[0]<=13]
    X=np.array([lam_at(x,y) for _,x,y in S])
    Y={}
    for deg in ANG:
        p,src,dep=load(runs_for(tag,deg))
        Y[deg]=np.array([float(((src==sid)&dep).sum()/max((src==sid).sum(),1)) for sid,_,_ in S])
    obs=np.mean([sp(X,Y[d]) for d in ANG])
    rng=np.random.default_rng(0)
    null=np.array([np.mean([sp(rng.permutation(X),Y[d]) for d in ANG]) for _ in range(20000)])
    print(f'  n={tag:>3}  f_dep 4방향 평균 ρs = {obs:+.3f}   p = {(np.abs(null)>=abs(obs)).mean():.4f}   귀무 표준편차 {null.std():.3f}')
