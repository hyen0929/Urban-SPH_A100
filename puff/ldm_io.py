"""LDM_<step>stp.vtk (POLYDATA, big-endian) 읽기: 좌표 + p_type(소스 태그)."""
import re
import numpy as np


def read_ldm(path):
    with open(path, 'rb') as f:
        buf = f.read()
    m = re.search(rb'POINTS\s+(\d+)\s+float\n', buf)
    n = int(m.group(1))
    xyz = np.frombuffer(buf, '>f4', 3 * n, m.end()).reshape(n, 3).astype(np.float32)
    out = {'xyz': xyz}
    pos = re.search(rb'FIELD\s+FieldData\s+(\d+)\n', buf).end()
    hdr = re.compile(rb'(\S+)\s+(\d+)\s+(\d+)\s+(\w+)\n')
    while True:
        h = hdr.match(buf, pos)
        if not h:
            break
        name, nc, nt, ty = h.group(1).decode(), int(h.group(2)), int(h.group(3)), h.group(4).decode()
        dt = {'float': '>f4', 'int': '>i4', 'double': '>f8'}[ty]
        a = np.frombuffer(buf, dt, nc * nt, h.end())
        out[name] = a.reshape(nt, nc) if nc > 1 else a
        pos = h.end() + a.nbytes
        if buf[pos:pos + 1] == b'\n':
            pos += 1
    return out
