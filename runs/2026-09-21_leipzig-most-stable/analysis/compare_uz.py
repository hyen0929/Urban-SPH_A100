import sys, numpy as np
sys.path.insert(0,'/tmp/claude-1014/-home-kdh-Urban-SPH/90fe7215-f210-4eca-bb99-e4c2226c7eef/scratchpad')
from vtkread import read

RUNS = {
 '중립 (09-15 sponge-uy)': 'runs/2026-09-15_leipzig-dt01-sponge-uy/out/fluid_condition_15000stp.vtk',
 '안정 T_con=285 (09-21)': 'runs/2026-09-21_leipzig-most-stable/out/fluid_condition_15000stp.vtk',
}
BANDS = [(0,5),(5,10),(10,20),(20,30),(30,50),(50,100),(100,200)]

data={}
for lab,f in RUNS.items():
    n,pts,fl = read(f)
    x=np.array(pts[0::3]); y=np.array(pts[1::3]); z=np.array(pts[2::3])
    v=np.array(fl['velocity']); pt=np.array(fl['p_type'])
    ux,uy,uz = v[0::3],v[1::3],v[2::3]
    m = (pt==1) & (np.abs(y)<1e-3)          # y=0 평면, 유체만
    data[lab]=dict(x=x[m],z=z[m],ux=ux[m],uy=uy[m],uz=uz[m])
    print(f"{lab}: y=0 평면 유체 {m.sum():,}개")

print("\n" + "="*78)
print(f"{'높이 [m]':>12} | {'평균 uz [m/s]':>26} | {'uz 표준편차':>22} | 하강비율")
print(f"{'':>12} | {'중립':>12} {'안정':>12} | {'중립':>10} {'안정':>10} | 중립  안정")
print("="*78)
for lo,hi in BANDS:
    row=[f"{lo:4d}~{hi:<4d}  "]
    means=[];stds=[];fdn=[]
    for lab in RUNS:
        d=data[lab]; s=(d['z']>=lo)&(d['z']<hi)
        if s.sum()==0: means.append(np.nan);stds.append(np.nan);fdn.append(np.nan); continue
        u=d['uz'][s]; means.append(u.mean()); stds.append(u.std()); fdn.append((u<0).mean()*100)
    print(f"{lo:5d}~{hi:<5d} | {means[0]:+12.4f} {means[1]:+12.4f} | {stds[0]:10.4f} {stds[1]:10.4f} | {fdn[0]:4.1f}% {fdn[1]:4.1f}%")
print("="*78)

print("\n=== 지면 최근접 층 (z<5 m) 상세 ===")
for lab in RUNS:
    d=data[lab]; s=d['z']<5
    u=d['uz'][s]; uh=np.hypot(d['ux'][s],d['uy'][s])
    print(f"{lab}")
    print(f"   uz : 평균 {u.mean():+.4f}  중앙 {np.median(u):+.4f}  p10 {np.percentile(u,10):+.4f}  p90 {np.percentile(u,90):+.4f}  |uz|평균 {np.abs(u).mean():.4f}")
    print(f"   수평: 평균 {uh.mean():.3f} m/s")
