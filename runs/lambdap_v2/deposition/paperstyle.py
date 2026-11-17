#!/usr/bin/env python3
"""논문 그림 공통 스타일 — 무게감 있는 박스 프레임, 안쪽 눈금, 충돌 없는 라벨 배치.

이전 판(dep_paper_figs*.py)의 문제: 7 pt 얇은 선, 점과 겹치는 숫자 라벨, 데이터 위에 놓인 통계 문구.
여기서는 (1) 본문 8 pt·선 0.9~1.4 pt·사방 박스 프레임, (2) 라벨을 8 방향 후보 중 충돌 없는 곳에
자동 배치, (3) 통계 문구는 반투명 박스로 데이터와 분리한다.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MM = 1/25.4
W1, W2 = 89*MM, 183*MM              # 단일단 / 양단
# Okabe–Ito 계열, 인쇄용으로 약간 어둡게
C = dict(blue='#005B99', orange='#D97706', green='#00785E', red='#C2410C',
         purple='#A3468A', sky='#3F9BD1', grey='#4A5560', ink='#111820')

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans', 'Arial'],
    'font.size': 8, 'axes.labelsize': 8.5, 'axes.titlesize': 8.5,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7.5,
    'axes.linewidth': 0.9, 'axes.edgecolor': C['ink'], 'axes.labelcolor': C['ink'],
    'text.color': C['ink'], 'xtick.color': C['ink'], 'ytick.color': C['ink'],
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.size': 3.6, 'ytick.major.size': 3.6,
    'xtick.minor.size': 2.0, 'ytick.minor.size': 2.0,
    'xtick.major.width': 0.9, 'ytick.major.width': 0.9,
    'xtick.minor.width': 0.7, 'ytick.minor.width': 0.7,
    'xtick.top': True, 'ytick.right': True,
    'axes.spines.top': True, 'axes.spines.right': True,
    'lines.linewidth': 1.4, 'lines.markersize': 4.5,
    'legend.frameon': True, 'legend.framealpha': .92, 'legend.edgecolor': '#C9D2D9',
    'legend.borderpad': .45, 'legend.handlelength': 1.6,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.015,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 130,
})


def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight='bold',
            va='top', ha='left', color=C['ink'])


def save(fig, name, outdir=None):
    outdir = outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper')
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name + '.pdf'))
    fig.savefig(os.path.join(outdir, name + '.png'), dpi=600)
    plt.close(fig)


def statbox(ax, text, loc='upper right', fontsize=7.5):
    """통계 문구를 반투명 박스에 넣어 데이터와 겹쳐도 읽히게 한다."""
    pos = {'upper right': (.97, .96, 'right', 'top'), 'upper left': (.03, .96, 'left', 'top'),
           'lower right': (.97, .04, 'right', 'bottom'), 'lower left': (.03, .04, 'left', 'bottom')}[loc]
    ax.text(pos[0], pos[1], text, transform=ax.transAxes, ha=pos[2], va=pos[3], fontsize=fontsize,
            bbox=dict(boxstyle='round,pad=0.32', fc='white', ec='#C9D2D9', lw=.7, alpha=.92), zorder=8)


def place_labels(ax, xs, ys, labels, fontsize=6.5, pad=7.0, color='#3B4650'):
    """점 라벨을 8 방향 후보 중 충돌이 가장 적은 위치에 놓는다(디스플레이 좌표 기준 탐욕 배치)."""
    ax.figure.canvas.draw()
    trans = ax.transData
    P = trans.transform(np.column_stack([xs, ys]))
    ab = ax.get_window_extent()                  # 축 경계 (디스플레이 좌표)
    placed = []                                  # 이미 놓인 라벨 bbox
    order = np.argsort(-np.asarray(ys, float))   # 위쪽부터
    cand = [(pad, pad), (-pad, pad), (pad, -pad), (-pad, -pad),
            (0, pad*1.25), (0, -pad*1.25), (pad*1.3, 0), (-pad*1.3, 0)]
    for idx in order:
        x0, y0 = P[idx]
        w = 4.6*fontsize/6.5*len(str(labels[idx])) + 3
        h = 8.5*fontsize/6.5
        best, best_cost = cand[0], 1e18
        for dx, dy in cand:
            bx = x0 + dx - (w if dx < 0 else 0)
            by = y0 + dy - (h if dy < 0 else 0)
            box = (bx, by, bx+w, by+h)
            cost = 0.0
            for q in P:                          # 다른 점과의 충돌
                if box[0]-2 < q[0] < box[2]+2 and box[1]-2 < q[1] < box[3]+2:
                    cost += 4
            for b in placed:                     # 다른 라벨과의 충돌
                if not (box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3]):
                    cost += 6
            if box[0] < ab.x0+1 or box[2] > ab.x1-1 or box[1] < ab.y0+1 or box[3] > ab.y1-1:
                cost += 12                        # 축 밖으로 나가면 큰 벌점
            cost += 0.02*abs(dx) + 0.02*abs(dy)
            if cost < best_cost:
                best_cost, best, best_box = cost, (dx, dy), box
        placed.append(best_box)
        ax.annotate(str(labels[idx]), (xs[idx], ys[idx]), textcoords='offset points',
                    xytext=best, fontsize=fontsize, color=color, zorder=7,
                    ha='left' if best[0] >= 0 else 'right',
                    va='bottom' if best[1] >= 0 else 'top')


def minor_ticks(ax, x=True, y=True):
    from matplotlib.ticker import AutoMinorLocator
    if x: ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    if y and ax.get_yscale() == 'linear': ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def scatter_sources(ax, x, y, ids, court=(3, 6, 10), label=True, s=30):
    """가로변(채움) / 중정(빈 기호) 산점도 + 라벨 자동 배치."""
    court = set(court)
    k = np.array([i not in court for i in ids])
    ax.scatter(np.asarray(x)[k], np.asarray(y)[k], s=s, c=C['blue'], lw=.6,
               edgecolors='white', zorder=5)
    ax.scatter(np.asarray(x)[~k], np.asarray(y)[~k], s=s+6, facecolors='white',
               edgecolors=C['red'], lw=1.3, zorder=5)
    if label:
        place_labels(ax, np.asarray(x, float), np.asarray(y, float), [str(i) for i in ids])
    return k
