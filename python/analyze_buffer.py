#!/usr/bin/env python3
"""Analyze p_type=0 (buffer) particles layer by layer."""
import numpy as np, re, glob

DTYPE_MAP = {"float":(">f4",4),"double":(">f8",8),"int":(">i4",4),"unsigned_int":(">u4",4)}
def parse_vtk(filepath):
    with open(filepath,"rb") as f: content=f.read()
    m=re.search(rb"POINTS\s+(\d+)\s+(\w+)\n",content)
    n=int(m.group(1)); dt,sz=DTYPE_MAP[m.group(2).decode("ascii")]
    pts=np.frombuffer(content,dtype=dt,count=n*3,offset=m.end()).reshape(n,3).astype(np.float64)
    fields={}
    mf=re.search(rb"FIELD FieldData\s+(\d+)\n",content)
    if not mf: return pts,fields
    off=mf.end()
    for _ in range(int(mf.group(1))):
        nl=content.index(b"\n",off); h=content[off:nl].decode("ascii").split()
        fn,nc,nt,ds=h[0],int(h[1]),int(h[2]),h[3]; off=nl+1
        dt2,isz=DTYPE_MAP.get(ds,(">f4",4))
        a=np.frombuffer(content,dtype=dt2,count=nc*nt,offset=off).astype(np.float64)
        if nc>1: a=a.reshape(nt,nc)
        fields[fn]=a; off+=nc*nt*isz
    return pts,fields

pts,fields=parse_vtk("Leipzig_small_point/fluid_condition_3000stp.vtk")
vel=fields["velocity"]; ptype=fields["p_type"].astype(int)
speed=np.sqrt(vel[:,0]**2+vel[:,1]**2+vel[:,2]**2)

buf=ptype==0
bp=pts[buf]; bv=vel[buf]; bs=speed[buf]

print("p_type=0 (buffer) z-layer별 속도")
print(f"  {'z':>8} | {'N':>6} | {'|u|_mean':>10} | {'ux':>10} | {'uy':>10} | {'uz':>10}")
for z in sorted(np.unique(np.round(bp[:,2],2))):
    m=np.abs(bp[:,2]-z)<0.1
    if m.sum()==0: continue
    print(f"  {z:>8.2f} | {m.sum():>6} | {bs[m].mean():>10.3f} | {bv[m,0].mean():>10.3f} | "
          f"{bv[m,1].mean():>10.3f} | {bv[m,2].mean():>10.3f}")

# 최하단(z=-8.75) vs 최상단(z=-1.25) vs 유체 최하층(z=1.25) 비교
print(f"\n속도 비교:")
fl = (ptype==1) & (np.abs(pts[:,2]-1.25)<0.5)
print(f"  z=-8.75 buffer:  ux={bv[np.abs(bp[:,2]-(-8.75))<0.5,0].mean():.3f}")
print(f"  z=-1.25 buffer:  ux={bv[np.abs(bp[:,2]-(-1.25))<0.5,0].mean():.3f}")
print(f"  z= 1.25 fluid:   ux={vel[fl,0].mean():.3f}")

# step별 buffer 속도 변화
print(f"\nstep별 z=-8.75 buffer 속도:")
files=sorted(glob.glob("Leipzig_small_point/fluid_condition_*stp.vtk"),
             key=lambda f:int(re.search(r"_(\d+)stp",f).group(1)))
print(f"  {'step':>6} | {'buf_ux(z=-8.75)':>16} | {'buf_ux(z=-1.25)':>16} | {'fluid_ux(z=1.25)':>18}")
for fp in files:
    st=int(re.search(r"_(\d+)stp",fp).group(1))
    if st%500!=0 and st!=0: continue
    p,f=parse_vtk(fp)
    v=f["velocity"]; pt=f["p_type"].astype(int)
    b0=pt==0
    bz8=(b0)&(np.abs(p[:,2]-(-8.75))<0.5)
    bz1=(b0)&(np.abs(p[:,2]-(-1.25))<0.5)
    fz1=(pt==1)&(np.abs(p[:,2]-1.25)<0.5)
    print(f"  {st:>6} | {v[bz8,0].mean():>16.3f} | {v[bz1,0].mean():>16.3f} | {v[fz1,0].mean():>18.3f}")
