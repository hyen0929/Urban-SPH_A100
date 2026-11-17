import numpy as np, sys
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/115d389f-4ea3-4674-a3fc-756b70aea151/scratchpad')
from vtkread import read_vtk
DX=2.5
d=read_vtk('runs/2026-10-16_ibm-ppe-wall-off/out/fluid_31500stp.vtk',want=('p_type','velocity'))
pts=d['points'].astype(np.float64); pt=d['p_type']; v=d['velocity'].astype(np.float64)
mark=pt==1000; fluid=pt==1
org=pts.min(0)
ix=np.floor((pts[:,0]-org[0])/DX).astype(np.int64); iy=np.floor((pts[:,1]-org[1])/DX).astype(np.int64)
NY=iy.max()+1; col=ix*NY+iy
ztop=np.full((ix.max()+1)*NY,-np.inf); np.maximum.at(ztop,col[mark],pts[mark,2])
inte=fluid&(pts[:,2]<=ztop[col]); ext=fluid&~inte
sp_i=np.linalg.norm(v[inte],axis=1); sp_e=np.linalg.norm(v[ext],axis=1)
print(f"건물 내부 유체 {inte.sum():,}개  /  외부 유체 {ext.sum():,}개\n")
print(f"{'':<22}{'내부':>12}{'외부':>12}")
print('-'*46)
for lab,f in [('평균 |u| [m/s]',np.mean),('중앙값 |u|',np.median)]:
    print(f"{lab:<22}{f(sp_i):>12.4f}{f(sp_e):>12.4f}")
for q in (90,99):
    print(f"{'p'+str(q)+' |u|':<22}{np.percentile(sp_i,q):>12.4f}{np.percentile(sp_e,q):>12.4f}")
print(f"{'최대 |u|':<22}{sp_i.max():>12.4f}{sp_e.max():>12.4f}")
print(f"{'|u| < 0.1 비율 %':<22}{(sp_i<0.1).mean()*100:>12.1f}{(sp_e<0.1).mean()*100:>12.1f}")
print(f"{'|u| > 1.0 비율 %':<22}{(sp_i>1.0).mean()*100:>12.1f}{(sp_e>1.0).mean()*100:>12.1f}")
# 깊이별 (표면에서 얼마나 안쪽인가)
depth=ztop[col][inte]-pts[inte,2]
print("\n건물 표면으로부터 깊이별 내부유체 |u|")
print(f"{'깊이 [m]':<14}{'n':>9}{'평균 |u|':>12}")
for lo,hi in [(0,2.5),(2.5,5),(5,10),(10,20),(20,1e9)]:
    m=(depth>=lo)&(depth<hi)
    if m.sum()==0: continue
    print(f"{lo:>5.1f}–{hi if hi<1e9 else 999:<8.1f}{m.sum():>9}{sp_i[m].mean():>12.4f}")
