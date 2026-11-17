import numpy as np
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
        sz={'float':4,'int':4,'double':8}[typ]; nb=nc*nt*sz
        if name in want:
            dt='>f4' if typ=='float' else '>i4'
            out[name]=np.frombuffer(f.read(nb),dtype=dt).reshape(nt,nc).squeeze().astype(
                np.float32 if typ=='float' else np.int32)
        else: f.seek(nb,1)
    f.close(); return out
