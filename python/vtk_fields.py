#!/usr/bin/env python3
"""VTK legacy binary(POLYDATA + FIELD)에서 각 필드의 min/max/NaN 개수를 뽑는다."""
import sys, re, numpy as np

def scan(path):
    with open(path, "rb") as f:
        buf = f.read()
    # FIELD FieldData N 이후부터 파싱
    m = re.search(rb"FIELD\s+FieldData\s+(\d+)\n", buf)
    if not m:
        return []
    pos, n = m.end(), int(m.group(1))
    out = []
    for _ in range(n):
        hm = re.compile(rb"([^\s]+)\s+(\d+)\s+(\d+)\s+(\w+)\n").match(buf, pos)
        if not hm:
            break
        name, ncomp, ntup, dtype = (hm.group(1).decode(), int(hm.group(2)),
                                    int(hm.group(3)), hm.group(4).decode())
        pos = hm.end()
        np_dt = {"float": ">f4", "double": ">f8", "int": ">i4"}.get(dtype)
        if np_dt is None:
            break
        nbytes = ncomp * ntup * np.dtype(np_dt).itemsize
        arr = np.frombuffer(buf, dtype=np_dt, count=ncomp * ntup, offset=pos)
        pos += nbytes
        if buf[pos:pos+1] == b"\n":
            pos += 1
        out.append((name, arr))
    return out

for path in sys.argv[1:]:
    print(f"── {path}")
    for name, a in scan(path):
        if a.dtype.kind == "f":
            nan = int(np.isnan(a).sum()); inf = int(np.isinf(a).sum())
            fin = a[np.isfinite(a)]
            rng = f"[{fin.min():.4g}, {fin.max():.4g}]" if fin.size else "(전부 비정상)"
            flag = f"  NaN {nan} / Inf {inf}" if (nan or inf) else ""
            print(f"   {name:<16} {rng}{flag}")
        else:
            print(f"   {name:<16} [{a.min()}, {a.max()}]  (int)")
