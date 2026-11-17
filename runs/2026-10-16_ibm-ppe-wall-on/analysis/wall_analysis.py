import numpy as np, sys

def read_vtk(path, want=('p_type','velocity','pressure')):
    f=open(path,'rb')
    def rl():
        b=bytearray()
        while True:
            c=f.read(1)
            if not c or c==b'\n': break
            b+=c
        return b.decode('latin1')
    for _ in range(4): rl()
    n=int(rl().split()[1])
    pts=np.frombuffer(f.read(n*3*4),dtype='>f4').reshape(n,3).astype(np.float32)
    rl(); k=int(rl().split()[-1])
    out={'points':pts,'n':n}
    for _ in range(k):
        p=rl().split(); name,nc,nt,typ=p[0],int(p[1]),int(p[2]),p[3]
        sz={'float':4,'int':4,'double':8}[typ]
        nbytes=nc*nt*sz
        if name in want:
            dt='>f4' if typ=='float' else '>i4'
            a=np.frombuffer(f.read(nbytes),dtype=dt)
            out[name]=a.reshape(nt,nc).squeeze().astype(np.float32 if typ=='float' else np.int32)
        else:
            f.seek(nbytes,1)
    f.close(); return out

DX=2.5
STEP=sys.argv[1] if len(sys.argv)>1 else '31500'
runs={t:f'runs/2026-10-16_ibm-ppe-wall-{t}/out/fluid_{STEP}stp.vtk' for t in ('off','on')}

print(f"=== step {STEP} ===")
D={}
for t,p in runs.items():
    D[t]=read_vtk(p); print(f"  {t}: N={D[t]['n']}")

# 입자 순서가 같은가 (Eulerian 이면 위치 고정 -> 정렬키 동일)
same = np.array_equal(D['off']['points'], D['on']['points'])
print(f"  위치 배열 동일: {same}")

pts=D['off']['points']; pt=D['off']['p_type']
mark = pt==1000
fl   = pt==1
print(f"  마커={mark.sum()}  유체(p_type=1)={fl.sum()}")

# --- 열 최대 높이장 (내부/외부 판별) ---
ix=np.floor((pts[:,0]-pts[:,0].min())/DX).astype(np.int64)
iy=np.floor((pts[:,1]-pts[:,1].min())/DX).astype(np.int64)
iz=np.floor((pts[:,2]-pts[:,2].min())/DX).astype(np.int64)
NX=ix.max()+1; NY=iy.max()+1
col=ix*NY+iy
ztop=np.full(NX*NY, -np.inf, dtype=np.float32)
np.maximum.at(ztop, col[mark], pts[mark,2])
exterior = pts[:,2] > ztop[col]        # 내부 <=> z <= z_top
print(f"  내부로 판정된 유체 = {(fl & ~exterior).sum()}  ({(fl&~exterior).sum()/fl.sum()*100:.2f}%)")

# --- 벽 인접 셀 (마커 셀의 26-이웃) ---
cell = (ix*NY+iy)*(iz.max()+1)+iz
NZ=iz.max()+1
mcells=np.unique(cell[mark])
# 마커 셀을 26방향으로 팽창
mi,mj,mk=np.unravel_index(mcells,(NX,NY,NZ))
offs=[(a,b,c) for a in(-1,0,1) for b in(-1,0,1) for c in(-1,0,1)]
dil=[]
for a,b,c in offs:
    ii,jj,kk=mi+a,mj+b,mk+c
    ok=(ii>=0)&(ii<NX)&(jj>=0)&(jj<NY)&(kk>=0)&(kk<NZ)
    dil.append((ii[ok]*NY+jj[ok])*NZ+kk[ok])
shell=np.unique(np.concatenate(dil))
near=np.isin(cell,shell,assume_unique=False)

band = fl & exterior & near
print(f"  벽 인접 외부유체 = {band.sum()}  ({band.sum()/ (fl&exterior).sum()*100:.2f}% of 외부유체)")
print()

hdr=f"{'지표':<30}{'OFF':>13}{'ON':>13}{'변화':>10}"
print(hdr); print('-'*len(hdr))
rows=[]
for t in ('off','on'):
    v=D[t]['velocity']; pr=D[t]['pressure']
    sp=np.linalg.norm(v[band].astype(np.float64),axis=1)
    spa=np.linalg.norm(v[fl&exterior].astype(np.float64),axis=1)
    rows.append(dict(
        mean_speed_wall=sp.mean(), p95_speed_wall=np.percentile(sp,95),
        max_speed_wall=sp.max(), clip_wall=(sp>=19.99).sum(),
        mean_speed_all=spa.mean(), clip_all=(spa>=19.99).sum(),
        mean_p_wall=pr[band].mean(), std_p_wall=pr[band].std(),
    ))
o,n=rows
for k_ in o:
    a,b=o[k_],n[k_]
    ch=f"{(b-a)/a*100:+9.1f}%" if a!=0 else "       —"
    print(f"{k_:<30}{a:>13.5g}{b:>13.5g}{ch:>10}")
