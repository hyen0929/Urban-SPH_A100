import numpy as np, sys
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/115d389f-4ea3-4674-a3fc-756b70aea151/scratchpad')
from vtkread import read_vtk

DX=2.5; R=1.5*DX            # 벽 첫 층 판정 반경 3.75 m
CS=R                         # 셀 크기 = 탐색 반경 -> 27 이웃으로 충분
STEP=sys.argv[1] if len(sys.argv)>1 else '31500'
D={t:read_vtk(f'runs/2026-10-16_ibm-ppe-wall-{t}/out/fluid_{STEP}stp.vtk') for t in ('off','on')}
pts=D['off']['points']; pt=D['off']['p_type']
mark=pt==1000; fl=pt==1
org=pts.min(0)

# 내부/외부: 열 최대 높이장
ix=np.floor((pts[:,0]-org[0])/DX).astype(np.int64)
iy=np.floor((pts[:,1]-org[1])/DX).astype(np.int64)
NY=iy.max()+1; col=ix*NY+iy
ztop=np.full((ix.max()+1)*NY,-np.inf,dtype=np.float32)
np.maximum.at(ztop,col[mark],pts[mark,2])
exterior=pts[:,2]>ztop[col]

# --- 마커 셀 해시 (고정 폭 패딩) ---
mp=pts[mark].astype(np.float64)
mc=np.floor((mp-org)/CS).astype(np.int64)
NXc,NYc,NZc=mc[:,0].max()+2,mc[:,1].max()+2,mc[:,2].max()+2
mkey=(mc[:,0]*NYc+mc[:,1])*NZc+mc[:,2]
order=np.argsort(mkey,kind='stable'); mkey_s=mkey[order]; mp_s=mp[order]
ukey,ustart,ucount=np.unique(mkey_s,return_index=True,return_counts=True)
CAP=int(ucount.max()); print(f"셀당 최대 마커 수 = {CAP}")
slots=np.full((len(ukey),CAP,3),np.inf)
rel=np.arange(len(mkey_s))-np.repeat(ustart,ucount)
slots[np.repeat(np.arange(len(ukey)),ucount),rel]=mp_s

cand=np.where(fl&exterior)[0]
cp=pts[cand].astype(np.float64)
cc=np.floor((cp-org)/CS).astype(np.int64)
best=np.full(len(cand),np.inf); bvec=np.zeros((len(cand),3))
for a in(-1,0,1):
  for b in(-1,0,1):
    for c in(-1,0,1):
      k=((cc[:,0]+a)*NYc+(cc[:,1]+b))*NZc+(cc[:,2]+c)
      p=np.searchsorted(ukey,k)
      hit=(p<len(ukey))&(ukey[np.minimum(p,len(ukey)-1)]==k)
      if not hit.any(): continue
      hi=np.where(hit)[0]
      S=slots[p[hi]]                                  # (m,CAP,3)
      dv=S-cp[hi][:,None,:]
      dd=np.sqrt((dv**2).sum(2))
      amin=np.argmin(dd,axis=1); dmin=dd[np.arange(len(hi)),amin]
      upd=dmin<best[hi]
      if upd.any():
          rows=hi[upd]
          best[rows]=dmin[upd]
          bvec[rows]=dv[np.arange(len(hi)),amin][upd]

sel=best<=R
idx=cand[sel]; dist=best[sel]; vec=bvec[sel]/best[sel][:,None]
print(f"\n=== step {STEP} ===")
print(f"벽 첫 층 외부유체 = {len(idx)}   거리 {dist.min():.2f}~{dist.max():.2f} m, 중앙값 {np.median(dist):.2f} m\n")
hdr=f"{'지표':<32}{'OFF':>12}{'ON':>12}{'변화':>10}"
print(hdr); print('-'*66)
res={}
for t in ('off','on'):
    v=D[t]['velocity'][idx].astype(np.float64)
    un=(v*vec).sum(1); sp=np.linalg.norm(v,axis=1)
    ut=np.sqrt(np.maximum(sp*sp-un*un,0))
    res[t]={'평균 u_n (+ = 벽 향함)':un.mean(),'평균 |u_n|':np.abs(un).mean(),
            'rms u_n':np.sqrt((un**2).mean()),'투과(u_n>0) 비율 %':(un>0).mean()*100,
            '투과입자 평균 u_n':un[un>0].mean(),'p95 |u_n|':np.percentile(np.abs(un),95),
            '평균 |u_t|':ut.mean(),'rms u_t / rms u_n':np.sqrt((ut**2).mean())/np.sqrt((un**2).mean())}
for k_ in res['off']:
    a,b=res['off'][k_],res['on'][k_]
    ch=f"{(b-a)/abs(a)*100:+9.1f}%" if a!=0 else "       —"
    print(f"{k_:<32}{a:>12.5g}{b:>12.5g}{ch:>10}")
