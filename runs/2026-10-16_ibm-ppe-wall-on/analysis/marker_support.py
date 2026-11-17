import numpy as np, sys
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/115d389f-4ea3-4674-a3fc-756b70aea151/scratchpad')
from vtkread import read_vtk
rng=np.random.default_rng(0)
DX=2.5; H_IBM=3.75; KAPPA=2.0; R=KAPPA*H_IBM; V=19.1406/1.225
A=1365.0/(64*np.pi*R**3)
def W(r):
    q=np.minimum(np.asarray(r)/R,1.0); t=1-q
    return A*t**8*(1+8*q+25*q*q+32*q**3)

d=read_vtk('runs/2026-10-16_ibm-ppe-wall-off/out/fluid_31500stp.vtk',want=('p_type',))
pts=d['points'].astype(np.float64); pt=d['p_type']
mark=pt==1000; fluid=pt==1
org=pts.min(0)
ix=np.floor((pts[:,0]-org[0])/DX).astype(np.int64); iy=np.floor((pts[:,1]-org[1])/DX).astype(np.int64)
NY=iy.max()+1; col=ix*NY+iy
ztop=np.full((ix.max()+1)*NY,-np.inf); np.maximum.at(ztop,col[mark],pts[mark,2])
interior=fluid&(pts[:,2]<=ztop[col])
print(f"마커 {mark.sum():,}   외부유체 {(fluid&~interior).sum():,}   내부유체 {interior.sum():,}\n")

fidx=np.where(fluid)[0]
CS=R
fc=np.floor((pts[fidx]-org)/CS).astype(np.int64)
NYc,NZc=fc[:,1].max()+2,fc[:,2].max()+2
fkey=(fc[:,0]*NYc+fc[:,1])*NZc+fc[:,2]
o=np.argsort(fkey,kind='stable'); fkey_s=fkey[o]; fpos=pts[fidx[o]]; fint=interior[fidx[o]]
uk,us,uc=np.unique(fkey_s,return_index=True,return_counts=True)

midx=np.where(mark)[0]
sample=rng.choice(midx,size=3000,replace=False)
mc=np.floor((pts[sample]-org)/CS).astype(np.int64)
OFF=[(a,b,c) for a in(-1,0,1) for b in(-1,0,1) for c in(-1,0,1)]
Si,Se,Oi,Oe,fr=[],[],[],[],[]
for s in range(len(sample)):
    xm=pts[sample[s]]; cc=mc[s]; parts=[]
    for a,b,c in OFF:
        k=((cc[0]+a)*NYc+(cc[1]+b))*NZc+(cc[2]+c)
        p=np.searchsorted(uk,k)
        if p<len(uk) and uk[p]==k: parts.append((us[p],us[p]+uc[p]))
    if not parts: continue
    sl=np.concatenate([np.arange(a_,b_) for a_,b_ in parts])
    dv=fpos[sl]-xm; r=np.sqrt((dv*dv).sum(1)); m=r<R
    if m.sum()==0: continue
    sl=sl[m]; dv=dv[m]; w=W(r[m])*V
    ii=fint[sl]
    wi=w.sum(); we=w[~ii].sum()
    if we<=0: continue
    Si.append(wi); Se.append(we)
    Oi.append(np.linalg.norm((dv*w[:,None]).sum(0)/wi))
    Oe.append(np.linalg.norm((dv[~ii]*w[~ii,None]).sum(0)/we))
    fr.append(ii.sum()/len(ii))
Si,Se,Oi,Oe,fr=map(np.array,(Si,Se,Oi,Oe,fr))
print(f"표본 마커 {len(Si):,}개,  내부입자를 하나라도 가진 마커 {(fr>0).mean()*100:.1f}%\n")
hdr=f"{'':<32}{'내부 포함':>14}{'내부 제외':>14}"
print(hdr); print('-'*60)
print(f"{'(AIJ 때)':<32}{'    ← 2026-03-09 이전':>14}{'  ← 현재':>14}")
print('-'*60)
print(f"{'커널합 ΣVW  평균':<32}{Si.mean():>14.4f}{Se.mean():>14.4f}")
print(f"{'커널합 ΣVW  중앙값':<32}{np.median(Si):>14.4f}{np.median(Se):>14.4f}")
print(f"{'커널합 < 0.9 인 마커 %':<32}{(Si<0.9).mean()*100:>14.1f}{(Se<0.9).mean()*100:>14.1f}")
print(f"{'중심 치우침 [m] 평균':<32}{Oi.mean():>14.4f}{Oe.mean():>14.4f}")
print(f"{'중심 치우침 [m] 중앙값':<32}{np.median(Oi):>14.4f}{np.median(Oe):>14.4f}")
print(f"{'중심 치우침 / Δx 중앙값':<32}{np.median(Oi)/DX:>14.3f}{np.median(Oe)/DX:>14.3f}")
