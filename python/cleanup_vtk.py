#!/usr/bin/env python3
"""
VTK 출력 정리 도구 — 시뮬레이션 결과 디렉토리의 용량을 줄인다.

기본 동작은 항상 dry-run. 실제 삭제는 --apply 를 명시해야만 일어난다.

  ./python/cleanup_vtk.py                        전체 디렉토리 용량 리포트
  ./python/cleanup_vtk.py --thin 10 CaseA        CaseA 를 10개 중 1개만 남기도록 (미리보기)
  ./python/cleanup_vtk.py --thin 10 --apply CaseA   실제 삭제
  ./python/cleanup_vtk.py --thin 10 --all        보호 디렉토리를 뺀 전체 (미리보기)

솎아내기는 파일 종류(fluid / LDM / vorticity / …)별로 따로 적용하며,
각 종류의 첫 스텝과 마지막 스텝은 항상 보존한다.
VTK 가 아닌 파일(CSV, PNG, 리포트)은 어떤 경우에도 건드리지 않는다.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# 삭제 금지 — ablation 실험 세트. 결론이 아직 문서화되지 않았다.
PROTECTED = {
    "CaseAS_neutral",
    "CaseAS_neutral_test1_no_flts",
    "CaseAS_neutral_test2_no_artdiff",
    "CaseAS_neutral_test3_no_sgs",
    "CaseAS_neutral_test4_no_conduction",
    "CaseAS_neutral_test5_eulert_off",
    "CaseAS_neutral_baseline_with_flts",
    "CaseAS_neutral_diag_eulert_on",
}

# fluid_1200stp.vtk → ("fluid", 1200)
VTK_RE = re.compile(r"^(?P<prefix>.+?)_(?P<step>\d+)stp\.vtk$")

# 스크립트는 python/ 에 있고 스캔 대상은 리포 루트다
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def human(n):
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024:
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}P"


def collect_vtks(d):
    """디렉토리 아래 모든 VTK 를 prefix 별로 모아 [(step, path, size)] 로 반환."""
    groups = defaultdict(list)
    for dirpath, _, filenames in os.walk(d):
        for fn in filenames:
            m = VTK_RE.match(fn)
            if not m:
                continue
            p = os.path.join(dirpath, fn)
            try:
                sz = os.path.getsize(p)
            except OSError:
                continue
            key = (dirpath, m.group("prefix"))
            groups[key].append((int(m.group("step")), p, sz))
    for v in groups.values():
        v.sort()
    return groups


def plan_thin(groups, keep_every):
    """솎아낸 뒤 삭제할 파일 목록을 만든다. 각 그룹의 첫/마지막 스텝은 보존."""
    doomed = []
    for (_, _), items in sorted(groups.items()):
        if len(items) <= 2:
            continue
        for idx, (_, path, size) in enumerate(items):
            if idx == 0 or idx == len(items) - 1:
                continue          # 첫 스텝 / 마지막 스텝 보존
            if idx % keep_every == 0:
                continue          # N개 중 1개 보존
            doomed.append((path, size))
    return doomed


def report():
    """전체 디렉토리를 용량 순으로 나열한다."""
    rows = []
    for name in sorted(os.listdir(ROOT)):
        d = os.path.join(ROOT, name)
        if not os.path.isdir(d) or name.startswith(".git"):
            continue
        vtk_n = vtk_b = other_n = other_b = 0
        for dirpath, _, filenames in os.walk(d):
            for fn in filenames:
                try:
                    sz = os.path.getsize(os.path.join(dirpath, fn))
                except OSError:
                    continue
                if fn.endswith(".vtk"):
                    vtk_n += 1
                    vtk_b += sz
                else:
                    other_n += 1
                    other_b += sz
        if vtk_n or other_n:
            rows.append((vtk_b + other_b, name, vtk_n, vtk_b, other_n, other_b))

    rows.sort(reverse=True)
    print(f"{'DIR':<40} {'TOTAL':>8} {'VTK':>8} {'VTK크기':>9} {'그외':>6} {'그외크기':>9}")
    print("-" * 86)
    tot = tot_vtk = 0
    for total, name, vn, vb, on, ob in rows:
        flag = "  [보호]" if name in PROTECTED else ""
        print(f"{name:<40} {human(total):>8} {vn:>8} {human(vb):>9} "
              f"{on:>6} {human(ob):>9}{flag}")
        tot += total
        tot_vtk += vb
    print("-" * 86)
    print(f"{'합계':<40} {human(tot):>8} {'':>8} {human(tot_vtk):>9}")
    print(f"\nVTK 가 전체의 {100 * tot_vtk / tot:.1f}% 를 차지합니다.")


def thin(targets, keep_every, apply_it):
    grand = 0
    for name in targets:
        if name in PROTECTED:
            print(f"[건너뜀] {name} — 보호 디렉토리", file=sys.stderr)
            continue
        d = os.path.join(ROOT, name)
        if not os.path.isdir(d):
            print(f"[건너뜀] {name} — 디렉토리 없음", file=sys.stderr)
            continue

        groups = collect_vtks(d)
        if not groups:
            print(f"[건너뜀] {name} — VTK 없음", file=sys.stderr)
            continue

        doomed = plan_thin(groups, keep_every)
        freed = sum(s for _, s in doomed)
        kept = sum(len(v) for v in groups.values()) - len(doomed)
        print(f"{name:<40} 삭제 {len(doomed):>6}개 / 보존 {kept:>5}개  → {human(freed):>8} 확보")
        grand += freed

        if apply_it:
            for path, _ in doomed:
                try:
                    os.remove(path)
                except OSError as e:
                    print(f"  삭제 실패 {path}: {e}", file=sys.stderr)

    print("-" * 86)
    verb = "확보함" if apply_it else "확보 예정"
    print(f"총 {human(grand)} {verb}")
    if not apply_it:
        print("\n실제로 삭제하려면 같은 명령에 --apply 를 붙이세요.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="대상 디렉토리 이름")
    ap.add_argument("--thin", type=int, metavar="N",
                    help="VTK 를 N개 중 1개만 남긴다 (첫/마지막 스텝은 항상 보존)")
    ap.add_argument("--all", action="store_true",
                    help="보호 디렉토리를 제외한 모든 출력 디렉토리를 대상으로")
    ap.add_argument("--apply", action="store_true",
                    help="실제로 삭제한다. 없으면 미리보기만 한다")
    args = ap.parse_args()

    if args.thin is None:
        report()
        return

    if args.thin < 2:
        ap.error("--thin 은 2 이상이어야 합니다")

    targets = args.dirs
    if args.all:
        targets = [n for n in sorted(os.listdir(ROOT))
                   if os.path.isdir(os.path.join(ROOT, n))
                   and not n.startswith(".git")
                   and n not in PROTECTED]
    if not targets:
        ap.error("대상 디렉토리를 지정하거나 --all 을 쓰세요")

    if args.apply:
        print(f"*** 실제 삭제 모드 — {len(targets)}개 디렉토리, {args.thin}개 중 1개 보존 ***\n")
    thin(targets, args.thin, args.apply)


if __name__ == "__main__":
    main()
