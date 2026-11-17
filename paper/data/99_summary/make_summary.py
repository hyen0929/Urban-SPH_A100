#!/usr/bin/env python3
"""모든 실험의 ρs 를 한 표로 모은다. paper/data/99_summary/master_correlations.csv
실행: python3 data/99_summary/make_summary.py   (paper/ 에서)
"""
import csv, json, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.path.normpath(os.path.join(HERE,'..','..'))
REPO=os.path.normpath(os.path.join(D,'..'))
ns={'__file__':f'{REPO}/runs/lambdap/deposition/dep_metrics.py'}
exec(open(ns['__file__']).read().split('def main(')[0], ns)
read,load,footprint,sp=ns['read_vtk'],ns['load'],ns['footprint'],ns['spearman']
ANG=[0,90,180,270]; CX=657.5
def rot(x,y,d):
    t=np.radians(d);c,s=np.cos(t),np.sin(t);dx,dy=x-CX,y
    return CX+dx*c-dy*s, dx*s+dy*c
hf=np.load(f'{D}/data/00_setup/heightfield.npz')
H,cell,x0,y0=hf['H'],float(hf['cell']),float(hf['x0']),float(hf['y0'])
B=(H>1.0).astype(float); ny,nx=H.shape
def lamfield(R):
    r=int(round(R/cell)); yy,xx=np.mgrid[-r:r+1,-r:r+1]
    K=((xx*xx+yy*yy)*cell*cell<=R*R).astype(float); K/=K.sum()
    Kp=np.zeros_like(B); Kp[:K.shape[0],:K.shape[1]]=K
    Kp=np.roll(np.roll(Kp,-(K.shape[0]//2),0),-(K.shape[1]//2),1)
    return np.fft.irfft2(np.fft.rfft2(B)*np.fft.rfft2(Kp),B.shape)
L75=lamfield(75.0)
def lam(x,y):
    j=int(np.floor((x-x0)/cell)); k=int(np.floor((y-y0)/cell))
    return float(L75[min(max(k,0),ny-1),min(max(j,0),nx-1)])
def srcs(f):
    return [(int(r['id']),float(r['x']),float(r['y'])) for r in csv.DictReader(open(f))]
def boot(x,y,n=4000,seed=0):
    rng=np.random.default_rng(seed); o=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(np.unique(x[i]))<4: continue
        o.append(sp(x[i],y[i]))
    return np.percentile(o,[2.5,97.5])

rows=[]
def emit(exp,n,flow,wall,metric,vals,ci=None):
    v=np.array(vals)
    rows.append(dict(experiment=exp,n_sources=n,flow=flow,wall_dep=wall,metric=metric,
        rho_0=round(v[0],3),rho_90=round(v[1],3),rho_180=round(v[2],3),rho_270=round(v[3],3),
        mean=round(v.mean(),3),spread=round(v.max()-v.min(),3),
        ci_width_mean=round(ci,3) if ci else ''))

SETS=[('n13_transient',13,'transient',0,f'{D}/data/00_setup/sources_n13.csv',
       {a:f'{REPO}/runs/'+{0:'2026-09-15_leipzig-dt01-sponge-uy',90:'2026-09-23_leipzig-rot90',
        180:'2026-09-23_leipzig-rot180',270:'2026-09-26_leipzig-rot270'}[a] for a in ANG}),
      ('n29_meanflow',29,'mean',0,f'{D}/data/00_setup/sources_n29.csv',
       {a:f'{REPO}/runs/2026-10-06_leipzig-src29b-rot{a}' for a in ANG}),
      ('n50_meanflow',50,'mean',0,f'{D}/data/00_setup/sources_n50.csv',
       {a:f'{REPO}/runs/2026-10-06_leipzig-src50-rot{a}' for a in ANG}),
      ('n50_meanflow_wall',50,'mean',1,f'{D}/data/00_setup/sources_n50.csv',
       {a:f'{REPO}/runs/2026-10-07_leipzig-src50w-rot{a}' for a in ANG})]
for exp,n,flow,wall,sf,rd in SETS:
    S=srcs(sf); ids=[s[0] for s in S]; X=np.array([lam(x,y) for _,x,y in S])
    per={}; gw={}
    for a in ANG:
        p,src,dep=load(rd[a])
        M={}
        for sid,x,y in S:
            m=src==sid
            fp=footprint(p[m&dep][:,:2],np.array(rot(x,y,a)))
            fp['f_dep']=float((m&dep).sum()/max(m.sum(),1)); M[sid]=fp
        per[a]=M
        if wall:
            _,F=read(f'{rd[a]}/out/LDM_30000stp.vtk')
            pt=F['p_type'].astype(int)-100; b=F['buffer_type'].astype(int)
            g=np.array([float(((pt==i)&(b==1)).sum()/max((pt==i).sum(),1)) for i in ids])
            w=np.array([float(((pt==i)&(b==2)).sum()/max((pt==i).sum(),1)) for i in ids])
            gw[a]=(g,w)
    for k in ['f_dep','area90_m2','f_r150','shift_m','elong']:
        v=[];c=[]
        for a in ANG:
            y=np.array([per[a][i][k] for i in ids]); v.append(sp(X,y))
            lo,hi=boot(X,y); c.append(hi-lo)
        emit(exp,n,flow,wall,k,v,np.mean(c))
    if wall:
        for lab,fn in (('f_ground',lambda g,w:g),('f_wall',lambda g,w:w),
                       ('f_dep_total',lambda g,w:g+w),('wall_to_ground',lambda g,w:w/np.maximum(g,1e-9))):
            v=[sp(X,fn(*gw[a])) for a in ANG]
            c=[np.diff(boot(X,fn(*gw[a])))[0] for a in ANG]
            emit(exp,n,flow,wall,lab,v,np.mean(c))
f=f'{HERE}/master_correlations.csv'
with open(f,'w',newline='') as g:
    w=csv.DictWriter(g,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f'저장: {f}  ({len(rows)} 행)')
for r in rows: print(f"  {r['experiment']:<20} {r['metric']:<14} {r['mean']:+.2f}  폭 {r['spread']:.2f}  CI {r['ci_width_mean']}")
