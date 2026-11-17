import sys, numpy as np
sys.path.insert(0,'/home/kdh/Urban-SPH/python')
from read_transient import Transient
R={'c':'runs/2026-11-03_cedval-c2-puff-dx15/out','M':'runs/2026-11-07_cedval-c2-puff-dx15-muscl/out'}
zb=np.arange(0,51,1.5); out={'zb':zb}
for k,o in R.items():
    T=Transient(o); X=T.xyz(); s=np.zeros(X.shape); s2=np.zeros(X.shape); n=0
    for c,t,u in T.frames():
        if c<20000 or c>99800 or c%400: continue
        u=u.astype(np.float64); s+=u; s2+=u*u; n+=1
    m=s/n; sd=np.sqrt(np.maximum(s2/n-m*m,0))
    for reg,(x0,x1) in {'up':(-1e9,41.25),'arr':(51.25,591.25)}.items():
        base=(X[:,0]>=x0)&(X[:,0]<x1)&(np.abs(X[:,1])<150)&(np.abs(m).sum(1)>1e-6)
        P=[]
        for a,b in zip(zb[:-1],zb[1:]):
            sel=base&(X[:,2]>=a)&(X[:,2]<b)
            P.append([m[sel,0].mean(),*sd[sel].mean(0),0.5*(sd[sel]**2).sum(1).mean()])
        out[f'{k}_{reg}']=np.array(P)
np.savez(sys.argv[1],**out); print('ok')
