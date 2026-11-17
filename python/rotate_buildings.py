#!/usr/bin/env python3
"""건물(IBM 마커)만 도메인 중심 기준으로 회전시켜 풍향 변경 실험 입력을 만든다.

풍향을 바꾸려면 유입·유출·측면 경계를 모두 옮겨야 해서 위험하다. 대신 **도시를 돌린다.**
이 케이스에서 안전한 이유:

  * 도메인이 정확히 정사각(x 0~1315, y ±657.5)이라 90° 회전이 도메인을 자기 자신으로 보낸다.
  * 격자 간격 2.5 m, 중심 (657.5, 0) 이 격자 위에 있어 유체 격자가 자기 자신으로 사상된다.
  * 지면(p_type 0)·천장(p_type -1)이 도메인 전체를 덮는 균일 층이고 지형 기복이 없다.
    즉 회전 대칭이므로 건드릴 필요가 없다.
  * 따라서 IBM 마커(p_type >= 1000)의 x, y 만 돌리면 경계조건 코드는 그대로다.

  python3 rotate_buildings.py <입력.txt> <출력.txt> [각도=90]
"""
import os
import sys
import math

SRC = sys.argv[1]
DST = sys.argv[2]
DEG = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
CX, CY = 657.5, 0.0                      # 도메인 중심 (x 0~1315, y ±657.5)

c, s = math.cos(math.radians(DEG)), math.sin(math.radians(DEG))
print(f'{SRC}\n  → {DST}\n  {DEG:g}° 회전, 중심 ({CX}, {CY})')

n_tot = n_rot = 0
xmin = ymin = 1e18
xmax = ymax = -1e18
with open(SRC) as fi, open(DST, 'w') as fo:
    fo.write(fi.readline())              # 열 인덱스 머리글 그대로
    for line in fi:
        n_tot += 1
        f = line.split('\t')
        if int(float(f[5])) >= 1000:     # 6번째 열 = p_type
            x, y = float(f[0]) - CX, float(f[1]) - CY
            nx, ny = CX + (x*c - y*s), CY + (x*s + y*c)
            f[0] = f'{nx:.6e}'; f[1] = f'{ny:.6e}'
            xmin = min(xmin, nx); xmax = max(xmax, nx)
            ymin = min(ymin, ny); ymax = max(ymax, ny)
            n_rot += 1
            fo.write('\t'.join(f))
        else:
            fo.write(line)
        if n_tot % 5_000_000 == 0:
            print(f'  {n_tot:,} 행', flush=True)

print(f'\n전체 {n_tot:,} 행,  회전한 IBM 마커 {n_rot:,} 개')
print(f'회전 후 건물 범위  x {xmin:.1f} ~ {xmax:.1f}   y {ymin:.1f} ~ {ymax:.1f}')
print(f'도메인            x 0.0 ~ 1315.0        y -657.5 ~ 657.5')
if xmin < 0 or xmax > 1315 or ymin < -657.5 or ymax > 657.5:
    print('!! 건물이 도메인을 벗어났다 — 각도를 줄이거나 배치를 조정할 것')
else:
    print(f'유입 fetch {xmin:.0f} m,  유출 여유 {1315-xmax:.0f} m,  '
          f'측면 여유 {ymin+657.5:.0f} / {657.5-ymax:.0f} m')
print(f'파일 크기 {os.path.getsize(DST)/1e9:.2f} GB')
