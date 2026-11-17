#!/usr/bin/env python3
"""낱장 그림용 스타일. 원칙 네 가지.

1. 그림 하나에 그래프 하나. 패널 합성 금지.
2. 범례는 축 바깥(위쪽)에 둔다 — 데이터와 겹칠 수 없다.
3. 선은 얇게(축 0.7 / 데이터 1.1), 마커는 흰 테두리로 분리.
4. 지도는 흰 바탕 + 연회색 건물 + 밝→진한 순차색. 검은 바탕 금지.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter

MM = 1/25.4
W_ONE = 88*MM          # 단일단 그림 기본 폭
H_ONE = 66*MM
C = dict(blue='#20618D', orange='#C97B1E', green='#2E7D5B', red='#B4442A',
         purple='#7D5BA6', sky='#5C9BC7', grey='#6B7580', light='#D9DEE3', ink='#1A2027')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 8, 'axes.labelsize': 8.5, 'axes.titlesize': 8.5,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
    'axes.linewidth': 0.7, 'axes.edgecolor': C['ink'], 'axes.labelcolor': C['ink'],
    'axes.labelpad': 3.5, 'text.color': C['ink'],
    'xtick.color': C['ink'], 'ytick.color': C['ink'],
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3.0, 'ytick.major.size': 3.0,
    'xtick.minor.size': 1.7, 'ytick.minor.size': 1.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.minor.width': 0.6, 'ytick.minor.width': 0.6,
    'xtick.major.pad': 3, 'ytick.major.pad': 3,
    'xtick.top': True, 'ytick.right': True,
    'axes.spines.top': True, 'axes.spines.right': True,
    'lines.linewidth': 1.1, 'lines.markersize': 4.6, 'lines.markeredgewidth': 0.5,
    'legend.frameon': True, 'legend.fancybox': False, 'legend.framealpha': 1.0,
    'legend.facecolor': 'white', 'legend.edgecolor': '#9AA6B0',
    'legend.handlelength': 1.5, 'legend.handletextpad': .5, 'legend.borderpad': .45,
    'legend.columnspacing': 1.1, 'legend.borderaxespad': 0, 'legend.labelspacing': .45,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 130,
})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper', 'figs')


def newfig(w=W_ONE, h=H_ONE):
    fig, ax = plt.subplots(figsize=(w, h))
    return fig, ax


def finish(fig, ax, name, minor=True, minor_y=True, tight=True):
    """tight=False 는 지정한 figure 크기 그대로 저장한다 — 여러 장의 축 박스를
    픽셀 단위로 맞춰야 하는 지도들에 쓴다(tight 크롭은 라벨 길이에 따라 폭이 달라진다)."""
    if minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if minor_y and ax.get_yscale() == 'linear':
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    os.makedirs(OUT, exist_ok=True)
    # bbox_inches=None 은 rcParams['savefig.bbox']='tight' 로 되돌아가므로
    # 전체 캔버스 bbox 를 명시해야 한다.
    from matplotlib.transforms import Bbox
    kw = {} if tight else dict(bbox_inches=Bbox([[0, 0], fig.get_size_inches()]), pad_inches=0)
    fig.savefig(os.path.join(OUT, name + '.pdf'), **kw)
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=600, **kw)
    plt.close(fig)


# ---------------------------------------------------------------- 축 안쪽 빈자리 찾기
# 범례·통계문구는 축 **안쪽**에 둔다. 자리는 실제로 그려진 화면을 읽어서 고른다:
# 캔버스를 래스터로 떠 잉크(흰 바탕과 다른 픽셀)의 적분영상을 만들고, 들어갈 상자를
# 축 안에서 격자 탐색해 잉크가 가장 적은 곳에 놓는다. 이미 놓은 상자는 피한다.
def _taken(ax):
    """이미 놓은 요소들(라벨·문구·범례). 상자가 아니라 **아티스트**로 들고 있는다 —
    축 범위를 넓히면 데이터에 붙은 라벨이 따라 움직이므로 좌표를 그때그때 다시 읽어야 한다."""
    if not hasattr(ax, '_placed_artists'):
        ax._placed_artists = []
    return ax._placed_artists


def _obstacles(ax):
    out = []
    for a in _taken(ax):
        try:
            e = a.get_window_extent()
            if e.width > 0 and e.height > 0: out.append((e.x0, e.y0, e.x1, e.y1))
        except Exception:
            pass
    return out


def _inkmap(ax):
    fig = ax.figure
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).astype(np.float32)
    rgb, a = buf[..., :3]/255.0, buf[..., 3:4]/255.0
    rgb = 1.0 - (1.0 - rgb)*a                      # 흰 바탕 위로 합성
    dark = (1.0 - rgb.min(axis=2))[::-1]           # 행 = 디스플레이 y (아래가 0)
    dark = np.clip((dark - 0.09)/0.91, 0, 1)       # 옅은 음영대(배경)는 자리를 막지 않는다
    S = np.zeros((dark.shape[0]+1, dark.shape[1]+1), np.float64)
    S[1:, 1:] = dark.cumsum(0).cumsum(1)
    return S


def _cover(S, x, y, w, h):
    H, W = S.shape[0]-1, S.shape[1]-1
    i0, j0 = int(max(y, 0)), int(max(x, 0))
    i1, j1 = int(min(y+h, H)), int(min(x+w, W))
    if i1 <= i0 or j1 <= j0: return 1.0
    return (S[i1, j1]-S[i0, j1]-S[i1, j0]+S[i0, j0])/((i1-i0)*(j1-j0))


def _best_spot(ax, w, h, prefer=None, S=None, corner=0.06, margin=0.030, clear=0.030):
    """폭 w·높이 h 상자를 축 안쪽 어디에 둘지. 반환: (x, y, 잉크비율) — (x,y)는 좌하단.

    corner  구석 선호 가중치(범례는 크게, 통계 문구는 0 — 자연스러운 빈 공간에 놓는다)
    margin  축 테두리와 띄울 여백(축 짧은 변 대비)
    clear   상자 둘레에 이만큼 더 비어 있기를 바란다 — 빠듯한 틈에 끼는 것을 막는다
    """
    if S is None: S = _inkmap(ax)
    bb = ax.get_window_extent()
    sm = min(bb.width, bb.height)
    m = margin*sm + 2
    x0, x1 = bb.x0+m, bb.x1-m-w
    y0, y1 = bb.y0+m, bb.y1-m-h
    if x1 <= x0 or y1 <= y0:                        # 축보다 크면 좌상단에 붙인다
        return bb.x0+m, bb.y1-m-h, 1.0
    taken = _obstacles(ax)
    cl = clear*sm
    pref = {'upper left': (0, 1), 'upper right': (1, 1),
            'lower left': (0, 0), 'lower right': (1, 0)}.get(prefer)
    best, bc, bcov = (x0, y1), 1e30, 1.0
    for fx in np.linspace(0, 1, 33):
        X = x0 + fx*(x1-x0)
        for fy in np.linspace(0, 1, 33):
            Y = y0 + fy*(y1-y0)
            cover = _cover(S, X, Y, w, h)
            # 잉크를 덮지 않는 것이 최우선. 구석 선호는 어디까지나 동점 처리용이라
            # 가중치를 확실히 낮게 둔다(예전엔 구석 보너스가 5 % 겹침을 이겨 버렸다).
            cost = 10.0*cover
            if cl > 0:                              # 둘레까지 비어 있는 자리를 우대
                cost += 3.0*_cover(S, X-cl, Y-cl, w+2*cl, h+2*cl)
            for (tx0, ty0, tx1, ty1) in taken:      # 이미 놓은 상자는 조금이라도 겹치면 안 된다
                ow = max(0.0, min(X+w, tx1+4)-max(X, tx0-4))
                oh = max(0.0, min(Y+h, ty1+4)-max(Y, ty0-4))
                # 이미 놓인 요소와의 겹침은 사실상 금지. 범례 안쪽은 흰 바탕이라
                # 잉크가 가장 적게 잡히므로, 벌점이 작으면 오히려 범례 위를 고른다.
                if ow*oh > 0: cost += 12.0 + 20.0*ow*oh/(w*h)
            cost += corner*(min(fx, 1-fx) + min(fy, 1-fy))
            if pref: cost += 0.06*(abs(fx-pref[0]) + abs(fy-pref[1]))
            if cost < bc: bc, best, bcov = cost, (X, Y), cover
    return best[0], best[1], bcov


def _to_axes(ax, x, y):
    return tuple(ax.transAxes.inverted().transform((x, y)))


def _headroom(ax, frac):
    """위쪽에 빈 띠를 만든다 — 그래프 안에 범례 자리가 없을 때의 정공법."""
    y0, y1 = ax.get_ylim()
    if ax.get_yscale() == 'log':
        ax.set_ylim(y0, y1*(y1/y0)**frac)
    else:
        ax.set_ylim(y0, y1 + frac*(y1-y0))


def legend_inside(ax, ncol=2, clean=0.0015, prefer=None, **kw):
    """범례를 축 **안쪽** 빈 공간에 놓는다. 열 수도 자리가 가장 깨끗한 쪽으로 고른다.

    빈 자리가 없으면 위쪽 여백을 늘려 자리를 만든다. 반투명으로 데이터를 가리는 것보다
    축을 넓히는 쪽이 낫다.
    """
    fig = ax.figure
    best = None
    for attempt in range(4):
        sizes = []
        for n in sorted({max(1, ncol), 2, 1}, reverse=True):
            leg = ax.legend(loc='upper left', ncol=n, **kw)
            fig.canvas.draw()
            e = leg.get_window_extent()
            leg.remove()
            if e.width < ax.get_window_extent().width - 6:
                sizes.append((n, e.width, e.height))
        fig.canvas.draw()
        S = _inkmap(ax)                              # 범례 없는 상태의 잉크 지도
        best = None
        for n, w, h in (sizes or [(1, 0.0, 0.0)]):
            x, y, cover = _best_spot(ax, w, h, S=S, prefer=prefer)
            if best is None or cover < best[0] - 1e-4:
                best = (cover, n, x, y, w, h)
        if best[0] <= clean or attempt == 3:
            break
        # 범례가 딱 들어갈 만큼만 위를 연다. 축을 f 배 늘리면 새로 생기는 빈 띠의 높이는
        # f/(1+f)·H 이므로, 필요한 띠 높이 r·H 에 대해 f = r/(1-r) 로 한 번에 맞춘다.
        # 여백은 **가장 납작한** 후보 기준으로 연다. 이번 판에서 세로로 긴 ncol=1 이
        # 뽑혔더라도 자리가 생기면 다음 판에서는 납작한 쪽이 선택되기 때문이다.
        ah = ax.get_window_extent().height
        hmin = min([h for _, _, h in sizes] or [best[5]])
        r = min(0.55, (hmin + 0.02*ah)/max(ah, 1))
        _headroom(ax, r/(1-r))
    cover, n, x, y, w, h = best
    leg = ax.legend(loc='lower left', bbox_to_anchor=_to_axes(ax, x, y), ncol=n, **kw)
    _taken(ax).append(leg)
    return leg


# 이전 이름 호환
legend_above = legend_inside


def note(ax, text, loc=None, dy=0.0, fontsize=8, corner=0.0):
    """통계 문구. 박스 없이, 축 안쪽에서 **자연스럽게** 비어 있는 자리에.

    범례와 달리 구석에 밀어붙이지 않는다(corner=0) — 데이터 옆 여백 한가운데가 낫다.

    여러 줄은 줄마다 별도 Text 로 그린다. matplotlib 은 한 문자열 안에 수식 위첨자가
    섞인 줄의 높이를 잘못 잡아 글자 아래를 잘라 먹는다($u_*$ = 0.83 m s$^{-1}$).
    """
    from matplotlib.transforms import IdentityTransform
    fig = ax.figure
    fig.canvas.draw()
    lines = str(text).split('\n')
    ts = [ax.text(0, 0, s, transform=IdentityTransform(), ha='left', va='bottom',
                  fontsize=fontsize, color=C['ink'], zorder=8) for s in lines]
    fig.canvas.draw()
    ex = [t.get_window_extent() for t in ts]
    lh = 1.42*fontsize*fig.dpi/72.0                  # 줄 간격
    w, h = max(e.width for e in ex), lh*len(lines)
    for t in ts: t.set_visible(False)                # 자기 자신은 잉크에서 빼고 찾는다
    x, y, _ = _best_spot(ax, w, h, prefer=loc, corner=corner, margin=0.055, clear=0.035)
    for k, t in enumerate(ts):
        t.set_visible(True); t.set_clip_on(False)
        # 디스플레이 좌표를 그대로 쓰면 600 dpi 저장 때 위치가 어긋난다 → 축 비율로 고정
        t.set_transform(ax.transAxes)
        t.set_position(_to_axes(ax, x, y + (len(lines)-1-k)*lh))
        _taken(ax).append(t)
    return ts


def log_minor(ax, axis='y'):
    a = ax.yaxis if axis == 'y' else ax.xaxis
    a.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)*0.1, numticks=20))
    a.set_minor_formatter(NullFormatter())


def label_points(ax, xs, ys, labels, fontsize=7, pad=4.2, color='#55606B'):
    """점 라벨을 겹치지 않게 배치. 축 한계를 먼저 정한 뒤 호출할 것."""
    fig = ax.figure
    fig.canvas.draw()
    s = fig.dpi/72.0                                 # pt → 디스플레이 px (offset points 와 단위를 맞춘다)
    P = ax.transData.transform(np.column_stack([xs, ys]))
    ab = ax.get_window_extent()
    placed, cand = [], []
    for f in (1.0, 1.45, 2.0, 2.8):                 # 가까운 자리가 막히면 조금씩 멀리
        cand += [(pad*f, pad*.8*f), (-pad*f, pad*.8*f), (pad*f, -pad*.8*f), (-pad*f, -pad*.8*f),
                 (0, pad*1.15*f), (0, -pad*1.15*f), (pad*1.2*f, 0), (-pad*1.2*f, 0)]
    M = 2.0*s                                        # 라벨끼리 최소 간격 — 맞닿으면 숫자가 붙어 읽힌다
    mk = 3.0*s                                       # 마커 반경 여유
    for idx in np.argsort(-np.asarray(ys, float)):
        x0, y0 = P[idx]
        w = (0.62*fontsize*len(str(labels[idx])) + 1.5)*s
        h = 1.15*fontsize*s
        best, bc, bb = cand[0], 1e18, None
        for dx, dy in cand:
            dxp, dyp = dx*s, dy*s
            bx = x0 + dxp - (w if dx < 0 else 0); by = y0 + dyp - (h if dy < 0 else 0)
            box = (bx, by, bx+w, by+h); cost = 0.0
            for q in P:
                if box[0]-mk < q[0] < box[2]+mk and box[1]-mk < q[1] < box[3]+mk: cost += 5
            for b in placed:
                if not (box[2]+M < b[0] or box[0]-M > b[2] or box[3]+M < b[1] or box[1]-M > b[3]): cost += 7
            if box[0] < ab.x0+2 or box[2] > ab.x1-2 or box[1] < ab.y0+2 or box[3] > ab.y1-2: cost += 15
            cost += .28*(abs(dx)+abs(dy))/pad       # 점에서 멀어지는 것을 강하게 억제
            if cost < bc: bc, best, bb = cost, (dx, dy), box
        placed.append(bb)
        far = abs(best[0]) + abs(best[1]) > 2.4*pad          # 멀리 놓인 라벨은 지시선을 붙인다
        an = ax.annotate(str(labels[idx]), (xs[idx], ys[idx]), textcoords='offset points', xytext=best,
                         fontsize=fontsize, color=color, zorder=7,
                         ha='left' if best[0] >= 0 else 'right', va='bottom' if best[1] >= 0 else 'top',
                         arrowprops=dict(arrowstyle='-', color=color, lw=.4,
                                         shrinkA=0.5, shrinkB=3.5) if far else None)
        _taken(ax).append(an)                                # 범례·문구가 라벨을 피하도록 등록


def sources_scatter(ax, x, y, ids, court=(3, 6, 10), s=26, label=True):
    court = set(court)
    k = np.array([i not in court for i in ids])
    x = np.asarray(x, float); y = np.asarray(y, float)
    ax.scatter(x[k], y[k], s=s, c=C['blue'], edgecolors='white', linewidths=.5, zorder=5,
               label='street-level source')
    ax.scatter(x[~k], y[~k], s=s+4, facecolors='white', edgecolors=C['red'], linewidths=1.0, zorder=5,
               label='enclosed courtyard')
    if label:
        label_points(ax, x, y, [str(i) for i in ids])
    return k


def expand(ax, fx=0.08, fy=0.12, logy=False, logx=False):
    """데이터가 축 가장자리에 붙지 않도록 여백을 준다(라벨 배치 전에 호출)."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    if logx:
        ax.set_xlim(x0*(1-fx*2), x1*(1+fx*4))
    else:
        d = x1-x0; ax.set_xlim(x0-fx*d, x1+fx*d)
    if logy:
        ax.set_ylim(y0*10**(-fy), y1*10**(fy*1.6))
    else:
        d = y1-y0; ax.set_ylim(y0-fy*d, y1+fy*1.3*d)
