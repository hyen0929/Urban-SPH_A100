#!/usr/bin/env python3
"""가벼운 순간 속도장 읽기 (functions_OUTPUT.cuh save_transient_velocity, -D TRANSIENT_DUMP_EVERY=N).

  <out>/transient_xyz.bin          좌표 (처음 한 번; 순서가 바뀌면 transient_xyz_<count>stp.bin 추가)
  <out>/transient_u_<count>stp.bin  (ux,uy,uz)
  머리 32 B: magic[8], int32 count, float32 time, int32 n, int32 pad, float64 key(= Σ(k+1)·x_k)

사용 (모듈):
    from read_transient import Transient
    T = Transient('runs/<run>/out')
    xyz = T.xyz()                       # [n,3]
    for step, t, u in T.frames():       # u [n,3]
        ...
    U = T.stack(every=1)                # [nt,n,3] 메모리 주의
사용 (명령줄): python3 python/read_transient.py runs/<run>/out   → 프레임 목록·크기 요약
"""
import glob, os, re, sys
import numpy as np

HDR = np.dtype([('magic', 'S8'), ('count', '<i4'), ('time', '<f4'), ('n', '<i4'), ('pad', '<i4'), ('key', '<f8')])
# 2026-11-04 이전 빌드: 필드 count/time 이 매크로 vii[27]/vif[3] 로 치환된 144 B 헤더.
#   count@116, n@128, key@136, time 은 기록 안 됨 → count·DT_OLD
HDR_OLD = np.dtype({'names': ['magic', 'count', 'n', 'key'], 'formats': ['S8', '<i4', '<i4', '<f8'],
                    'offsets': [0, 116, 128, 136], 'itemsize': 144})
DT_OLD = 0.01


def _head(path):
    size = os.path.getsize(path)
    with open(path, 'rb') as f:
        b = f.read(HDR_OLD.itemsize)
    h = np.frombuffer(b[:HDR.itemsize], HDR)[0]
    if size == HDR.itemsize + 12 * int(h['n']):
        return h, HDR.itemsize
    o = np.frombuffer(b, HDR_OLD)[0]
    if size != HDR_OLD.itemsize + 12 * int(o['n']):
        raise ValueError(f'{path}: 헤더를 알 수 없음')
    h = np.zeros((), HDR)
    h['magic'], h['count'], h['n'], h['key'] = o['magic'], o['count'], o['n'], o['key']
    h['time'] = o['count'] * DT_OLD
    return h, HDR_OLD.itemsize


def _read(path):
    h, off = _head(path)
    a = np.fromfile(path, '<f4', count=3 * int(h['n']), offset=off).reshape(-1, 3)
    return h, a


class Transient:
    def __init__(self, out):
        self.out = out
        step = lambda f: int(re.search(r'_(\d+)stp', f).group(1))
        self.ufiles = sorted(glob.glob(f'{out}/transient_u_*stp.bin'), key=step)
        self.xfiles = glob.glob(f'{out}/transient_xyz*.bin')
        if not self.xfiles:
            raise FileNotFoundError(f'{out}/transient_xyz.bin 없음')
        self._xyz = {}                       # key → 좌표

    def _xyz_for(self, h):
        k = float(h['key'])
        for kk, v in self._xyz.items():
            if abs(kk - k) <= 1e-9 * abs(k):
                return v
        for f in self.xfiles:
            hx = _head(f)[0]
            if int(hx['n']) == int(h['n']) and abs(float(hx['key']) - k) <= 1e-9 * abs(k):
                self._xyz[float(hx['key'])] = _read(f)[1]
                return self._xyz[float(hx['key'])]
        raise ValueError(f"step {int(h['count'])}: 맞는 좌표 파일이 없음 (n={int(h['n'])})")

    def xyz(self, i=0):
        return self._xyz_for(_head(self.ufiles[i])[0])

    def frames(self, every=1):
        for f in self.ufiles[::every]:
            h, u = _read(f)
            yield int(h['count']), float(h['time']), u

    def stack(self, every=1):
        return np.stack([u for _, _, u in self.frames(every)])


if __name__ == '__main__':
    T = Transient(sys.argv[1])
    h = _head(T.ufiles[0])[0]; hl = _head(T.ufiles[-1])[0]
    n = int(h['n'])
    print(f"프레임 {len(T.ufiles)}  step {int(h['count'])}–{int(hl['count'])}  t {h['time']:.2f}–{hl['time']:.2f} s")
    print(f"입자 {n:,}  프레임당 {12 * n / 1e6:.1f} MB  좌표 파일 {len(T.xfiles)}")
    xyz = T.xyz()
    print('범위', xyz.min(0), xyz.max(0))
