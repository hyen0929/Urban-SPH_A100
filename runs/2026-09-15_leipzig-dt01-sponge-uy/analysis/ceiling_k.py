#!/usr/bin/env python3
"""천장 불안정 추적 — 슬래브(fluid_condition, y∈[0,7.5))에서 z>=250 m Σk_turb, z>=355 m |uz| 최대,
z 120~250 m Σk, z>=340 m uy 표준편차를 새 슬래브마다 analysis/ceiling_k.csv 에 추가하고 최근 행을 출력한다.
스펀지 없는 분기 원런(2026-09-14_leipzig-dt01-full): z>=250 Σk = 36(t=90) → 314(110) → 1752(130), |uz|max 1.3~1.6.
dt=0.005 Run A: 8.6(120) → 11(150) → 31(200) → 125(236), |uz|max ≤0.2, z 340~372 uy std 0.055(t=130).
스펀지 uz 전용 분기(2026-09-15_leipzig-dt01-sponge): uy std 0.31(110) → 0.37(120) → 0.49(130), Σk z>=250 425(120) → 555(130).
"""
import re, os, glob, numpy as np
RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, CSV = os.path.join(RUN, 'out'), os.path.join(RUN, 'analysis', 'ceiling_k.csv')

def read(path):
    buf = open(path, 'rb').read()
    m = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", buf); n = int(m.group(1))
    dt = {"float": ">f4", "double": ">f8"}[m.group(2).decode()]
    pts = np.frombuffer(buf, dtype=dt, count=3*n, offset=m.end()).reshape(-1, 3).astype(np.float64)
    f = re.search(rb"FIELD\s+FieldData\s+(\d+)\n", buf); pos = f.end(); F = {}
    for _ in range(int(f.group(1))):
        hm = re.compile(rb"([^\s]+)\s+(\d+)\s+(\d+)\s+(\w+)\n").match(buf, pos)
        if not hm: break
        name, nc, nt, ty = hm.group(1).decode(), int(hm.group(2)), int(hm.group(3)), hm.group(4).decode()
        pos = hm.end(); npd = {"float": ">f4", "double": ">f8", "int": ">i4"}[ty]
        a = np.frombuffer(buf, dtype=npd, count=nc*nt, offset=pos).astype(np.float64)
        pos += nc*nt*np.dtype(npd).itemsize
        if buf[pos:pos+1] == b"\n": pos += 1
        F[name] = a.reshape(nt, nc) if nc > 1 else a
    return pts, F

done = set()
if os.path.exists(CSV):
    done = {int(l.split(',')[0]) for l in open(CSV).read().splitlines()[1:] if l}
else:
    open(CSV, 'w').write('step,time,k_sum_z250,uz_max_z355,k_sum_z120_250,uy_std_z340\n')
import time
for p in sorted(glob.glob(os.path.join(OUT, 'fluid_condition_*stp.vtk')), key=lambda s: int(re.search(r'_(\d+)stp', s).group(1))):
    step = int(re.search(r'_(\d+)stp', p).group(1))
    if step in done or os.path.getmtime(p) > time.time() - 30:
        continue
    pts, F = read(p)
    fl = F['p_type'] == 1; z = pts[:, 2]; k = F['k_turb']; v = F['velocity']
    row = (step, step*0.01, k[fl & (z >= 250)].sum(), np.abs(v[fl & (z >= 355), 2]).max(), k[fl & (z >= 120) & (z < 250)].sum(),
           v[fl & (z >= 340), 1].std())
    with open(CSV, 'a') as f:
        f.write('%d,%.2f,%.4e,%.4f,%.4e,%.4f\n' % row)
rows = open(CSV).read().splitlines()[1:]
for l in rows[-4:]:
    print(l)
