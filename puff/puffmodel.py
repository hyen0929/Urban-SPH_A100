"""SPH 평균장 위의 가우시안 퍼프 모델 (Leipzig, 중립).

퍼프 하나 = 질량 m, 중심 X, 등방 크기 sigma.
  중심 이동 : dX = (u_mean(X) + u') dt,  u' 는 LDM 과 같은 OU 과정이지만 분산을 beta 만큼만 준다
              u' <- u' R + sqrt(beta) sigma_u sqrt(1-R^2) xi,  R = exp(-dt/T_L)
  크기 성장 : d(sigma^2)/dt = 2 (1-beta) sigma_u^2 T_L (1 - exp(-age/T_L))
              → 나머지 (1-beta) 몫의 난류 확산은 퍼프가 커지는 것으로 표현. 장시간 극한에서
                총 확산 = beta K + (1-beta) K = K (LDM 과 같음). beta 는 "출렁임 vs 매끈함" 비율.
  벽        : 건물 칸에 들어가면 그 축 변위·u' 를 뒤집는다(축별 반사). 지면 z=0 반사.
  침적      : dm/dt = -vd · m · 2 g(0),  g = 퍼프 연직 가우시안의 지면 밀도(거울상 포함)
  농도      : 퍼프 가우시안 합 + 지면 거울상. 건물 칸은 0.

sigma_u = sqrt(2k/3), T_L = 2 sigma_u^2/(C0 eps) — LDM(KERNEL_time_update_LDM) 정의 그대로.
"""
import numpy as np


class Field:
    def __init__(self, path=None, uniform=None):
        if uniform is not None:                        # 건물 없는 균일장 (비교용 '보통 퍼프')
            self.uniform = uniform                     # dict(u=, sig=, TL=)
            self.solid = None
            return
        self.uniform = None
        d = np.load(path)
        U, sig = d['U'], d['sig']
        self.F = np.ascontiguousarray(np.concatenate([U, sig[None]], 0).transpose(1, 2, 3, 0))  # nx,ny,nz,5
        self.solid = d['solid']; self.H = d['H']
        self.o = d['origin'].astype(np.float64); self.dx = float(d['dx'])
        self.n = np.array(self.solid.shape)
        self.us = float(d['us'])

    def sample(self, p):
        """p (N,3) → ux,uy,uz,TL,sig (N,5), 삼선형 보간."""
        if self.uniform is not None:
            u = self.uniform
            return np.tile([u['u'], 0.0, 0.0, u['TL'], u['sig']], (len(p), 1))
        f = (p - self.o) / self.dx
        i0 = np.clip(np.floor(f).astype(int), 0, self.n - 2)
        w = np.clip(f - i0, 0.0, 1.0)
        out = 0.0
        for a in (0, 1):
            wa = w[:, 0] if a else 1 - w[:, 0]
            for b in (0, 1):
                wb = w[:, 1] if b else 1 - w[:, 1]
                for c in (0, 1):
                    wc = w[:, 2] if c else 1 - w[:, 2]
                    out = out + (wa * wb * wc)[:, None] * self.F[i0[:, 0] + a, i0[:, 1] + b, i0[:, 2] + c]
        return out

    def is_solid(self, p):
        if self.solid is None:
            return np.zeros(len(p), bool)
        ii = np.rint((p - self.o) / self.dx).astype(int)
        ok = np.all((ii >= 0) & (ii < self.n), 1)
        s = np.zeros(len(p), bool)
        s[ok] = self.solid[ii[ok, 0], ii[ok, 1], ii[ok, 2]]
        return s

    def solid_slice(self, xc, yc, z):
        if self.solid is None:
            return np.zeros((len(xc), len(yc)), bool)
        i = np.rint((xc - self.o[0]) / self.dx).astype(int)
        j = np.rint((yc - self.o[1]) / self.dx).astype(int)
        k = int(np.clip(np.rint((z - self.o[2]) / self.dx), 0, self.n[2] - 1))
        return self.solid[np.ix_(i, j, [k])][:, :, 0]


DOMAIN = (0.0, 1315.0, -657.5, 657.5, 400.0)          # x0, x1, y0, y1, 천장


def simulate(field, src, t_start, t_stop, t_end, out_times, grid, zlev,
             beta=0.5, vd=0.0, dt=0.2, dt_rel=0.5, sigma0=2.5, seed=0, record=None):
    """src (Ns,3) 소스 좌표. [t_start, t_stop) 동안 연속 방출(소스당 방출률 1 → c 는 c/Q).
    out_times 의 각 시각에 grid=(xc, yc) 평면, zlev 높이들의 농도장을 돌려준다.
    반환: C (nt, nz, Ns, nx, ny) [s/m^3], dep (nt, Ns) 누적 침적 비율.
    record: 시각 → 콜백(t, state) (드론 학습용 순간장 저장 등)."""
    rng = np.random.default_rng(seed)
    Ns = len(src)
    x0, x1, y0, y1, ztop = DOMAIN
    xc, yc = grid
    n_cap = int(np.ceil((t_stop - t_start) / dt_rel)) * Ns
    P = np.zeros((n_cap, 3)); V = np.zeros((n_cap, 3)); S2 = np.zeros(n_cap); M = np.zeros(n_cap)
    AGE = np.zeros(n_cap); SRC = np.zeros(n_cap, int); alive = np.zeros(n_cap, bool)
    n_used = 0
    dep = np.zeros(Ns)
    out_times = list(out_times)
    C = np.zeros((len(out_times), len(zlev), Ns, len(xc), len(yc)), np.float32)
    D = np.zeros((len(out_times), Ns))
    solid_masks = [field.solid_slice(xc, yc, z) for z in zlev]
    next_rel = t_start
    t = t_start
    it_out = 0
    nsteps = int(round((t_end - t_start) / dt))
    for step in range(nsteps + 1):
        t = t_start + step * dt
        # --- 방출
        while next_rel < t_stop and next_rel <= t + 1e-9:
            sl = slice(n_used, n_used + Ns)
            P[sl] = src; V[sl] = 0; S2[sl] = sigma0 ** 2; M[sl] = dt_rel; AGE[sl] = 0
            SRC[sl] = np.arange(Ns); alive[sl] = True
            n_used += Ns; next_rel += dt_rel
        # --- 출력
        while it_out < len(out_times) and out_times[it_out] <= t + 1e-9:
            C[it_out] = concentration(P, S2, M, SRC, alive, Ns, xc, yc, zlev, solid_masks)
            D[it_out] = dep / ((t_stop - t_start))
            if record is not None:
                record(t, dict(P=P[alive], S2=S2[alive], M=M[alive], SRC=SRC[alive]))
            it_out += 1
        if step == nsteps:
            break
        # --- 이동
        a = np.where(alive)[0]
        if len(a) == 0:
            continue
        p = P[a]
        f = field.sample(p)
        um, TL, su = f[:, :3], np.maximum(f[:, 3], 1e-3), f[:, 4]
        R = np.exp(-dt / TL)[:, None]
        V[a] = V[a] * R + (np.sqrt(beta) * su)[:, None] * np.sqrt(1 - R ** 2) * rng.standard_normal((len(a), 3))
        d = (um + V[a]) * dt
        # 축별 반사: x, y, z 순서로 한 축씩 움직이며 건물 칸에 들어가면 그 축을 뒤집는다
        q = p.copy()
        for ax in range(3):
            trial = q.copy(); trial[:, ax] += d[:, ax]
            hit = field.is_solid(trial)
            if hit.any():
                d[hit, ax] *= -1; V[a[hit], ax] *= -1
                trial[hit, ax] = q[hit, ax] + d[hit, ax]
                still = field.is_solid(trial)
                trial[still, ax] = q[still, ax]          # 반사해도 막히면 그 축은 정지
            q = trial
        neg = q[:, 2] < 0
        q[neg, 2] *= -1; V[a[neg], 2] *= -1
        P[a] = q
        # --- 크기 성장
        AGE[a] += dt
        S2[a] += 2 * (1 - beta) * su ** 2 * TL * (1 - np.exp(-AGE[a] / TL)) * dt
        # --- 침적
        if vd > 0:
            s = np.sqrt(S2[a])
            g0 = 2 * np.exp(-q[:, 2] ** 2 / (2 * S2[a])) / (np.sqrt(2 * np.pi) * s)
            lost = M[a] * (1 - np.exp(-vd * g0 * dt))
            M[a] -= lost
            np.add.at(dep, SRC[a], lost)
        # --- 도메인 이탈
        gone = (q[:, 0] < x0) | (q[:, 0] > x1) | (q[:, 1] < y0) | (q[:, 1] > y1) | (q[:, 2] > ztop)
        alive[a[gone]] = False
    return C, D


def concentration(P, S2, M, SRC, alive, Ns, xc, yc, zlev, solid_masks, cut=4.0):
    a = np.where(alive)[0]
    out = np.zeros((len(zlev), Ns, len(xc), len(yc)), np.float32)
    if len(a) == 0:
        return out
    p, s2, m, sr = P[a], S2[a], M[a], SRC[a]
    s = np.sqrt(s2)
    norm = m / ((2 * np.pi) ** 1.5 * s ** 3)
    gx = np.exp(-(xc[None, :] - p[:, :1]) ** 2 / (2 * s2[:, None]))
    gy = np.exp(-(yc[None, :] - p[:, 1:2]) ** 2 / (2 * s2[:, None]))
    order = np.argsort(sr, kind='stable'); bounds = np.searchsorted(sr[order], np.arange(Ns + 1))
    for iz, z in enumerate(zlev):
        gz = np.exp(-(z - p[:, 2]) ** 2 / (2 * s2)) + np.exp(-(z + p[:, 2]) ** 2 / (2 * s2))
        A = norm * gz
        for k in range(Ns):
            idx = order[bounds[k]:bounds[k + 1]]
            if len(idx):
                out[iz, k] = (gx[idx] * A[idx, None]).T @ gy[idx]
        out[iz][:, solid_masks[iz]] = 0
    return out
