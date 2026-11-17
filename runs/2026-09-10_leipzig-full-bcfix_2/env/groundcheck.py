#!/usr/bin/env python3
"""지표층(z=1.25 m) 속도 분포 감시.

tu_max 는 t=6.5 s 부터 클립(20 m/s)에 붙어 포화돼 건전성 신호로 못 쓴다. 대신
fluid_condition 슬래브(y∈[0,7.5))의 지표층 1,607 개 입자 **분포**를 30 분마다 기록하고,
스핀업 피크(t=17.5 s: p99=16.5, >15 m/s 19개) 수준으로 되돌아가면 ALERT 를 낸다.

  setsid nohup runs/2026-09-10_leipzig-full-bcfix_2/env/groundcheck.py >/dev/null 2>&1 &
"""
import re, os, glob, time, datetime, numpy as np

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, CSV = os.path.join(RUN, "out"), os.path.join(RUN, "analysis", "ground_layer.csv")
WD  = os.path.join(RUN, "analysis", "watchdog.log")
INTERVAL, RESUPPRESS = 1800, 6*3600
# 임계 재조정 (2026-09-12 23:20, t=200 s 전체장 분석 후):
#   지표층 이상치는 진폭이 큰 간헐현상이다 — p99 가 6.5(t=193.5) ~ 15.2(t=182) 를 왕복한다.
#   그래서 첫 임계(p99>14, >15m/s 12개)는 정상 진동의 봉우리에서도 울렸다. 진짜 신호는
#   (1) 봉우리가 더 높아지는 것과 (2) **중앙값이 계속 올라가는 것**(t=150 1.535 → t=197 1.732).
P99_LIM, N15_LIM, P50_LIM = 18.0, 30, 2.2

def say(msg):
    with open(WD, "a") as f:
        f.write(f"[{datetime.datetime.now():%F %T}] {msg}\n")

def read(path):
    buf = open(path, "rb").read()
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

if not os.path.exists(CSV):
    open(CSV, "w").write("step,time,N,mean,p50,p99,max,n_gt15,n_gt19,slab_mean,slab_max,slab_max_uz\n")
say(f"INFO  groundcheck 시작 (pid={os.getpid()}) — 30분 주기, "
    f"ALERT 조건 중앙값>{P50_LIM} 또는 p99>{P99_LIM} 또는 >15m/s {N15_LIM}개↑")

last_step, last_alert = -1, 0.0
while True:
    try:
        files = glob.glob(os.path.join(OUT, "fluid_condition_*stp.vtk"))
        if files:
            step = max(int(re.search(r"_(\d+)stp", f).group(1)) for f in files)
            path = os.path.join(OUT, f"fluid_condition_{step}stp.vtk")
            if step > last_step and os.path.getmtime(path) < time.time() - 30:   # 쓰기 완료 대기
                pts, F = read(path)
                v = F["velocity"]; u = np.linalg.norm(v, axis=1)
                g = np.abs(pts[:, 2] - 1.25) < 0.1; ug = u[g]
                row = (step, step*0.005, ug.size, ug.mean(), np.median(ug), np.percentile(ug, 99),
                       ug.max(), int((ug > 15).sum()), int((ug > 19).sum()),
                       u.mean(), u.max(), np.abs(v[:, 2]).max())
                with open(CSV, "a") as f:
                    f.write("%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%d,%d,%.4f,%.4f,%.4f\n" % row)
                if (row[5] > P99_LIM or row[7] >= N15_LIM or row[4] > P50_LIM) \
                        and time.time() - last_alert > RESUPPRESS:
                    say(f"ALERT 지표층 악화 — step={step} t={row[1]:.1f} s: 중앙값={row[4]:.2f}, "
                        f"p99={row[5]:.2f} m/s, >15m/s {row[7]}개, >19m/s {row[8]}개 "
                        f"(기준: 중앙값>{P50_LIM} / p99>{P99_LIM} / >15m/s {N15_LIM}개↑)")
                    last_alert = time.time()
                last_step = step
    except Exception as e:
        say(f"INFO  groundcheck 예외(무시하고 계속): {type(e).__name__}: {e}")
    time.sleep(INTERVAL)
