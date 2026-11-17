#!/bin/bash
# exp/status.sh — 모든 실험의 현재 상태를 한 화면에.
#   exp/status.sh          전체
#   exp/status.sh <run>    하나만 자세히
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

detail() {
  local d; d=$(resolve_run "$1") || exit 1
  local id; id=$(basename "$d")
  echo "== $id  [$(run_state "$d")]"
  echo
  sed -n '/^## 1. 목적/,/^## 2\./p' "$d/RUN.md" 2>/dev/null | sed '1d;$d' | grep -v '^<!--' | sed '/^$/d' | sed 's/^/  목적: /'
  echo "  플래그: $(paste -sd' ' "$d/env/flags.txt" 2>/dev/null)"
  echo "  진행:   step=$(last_step "$d")  sim=$(last_simtime "$d")s  wall=$(last_walltime "$d")s"
  echo "  최근:   $(grep -E '^[0-9]+\s+tu_max' "$d/run.log" 2>/dev/null | tail -1)"
  echo "  출력:   $(ls "$d/out"/*.vtk 2>/dev/null | wc -l) vtk, $(du -sh "$d/out" 2>/dev/null | cut -f1)"
  echo "  로그:   $(basename "$d")/run.log"
  local bad; bad=$(grep -icE "nan|segmentation|cudaError|illegal|error" "$d/run.log" 2>/dev/null)
  [ "${bad:-0}" -gt 0 ] && echo "  ⚠ 경고:  로그에 nan/error 문자열 $bad 건 — grep 확인 필요"
  echo
}

[ $# -gt 0 ] && { detail "$1"; exit 0; }

shopt -s nullglob
runs=("$RUNS_DIR"/*/)
[ ${#runs[@]} -eq 0 ] && { echo "runs/ 가 비어 있다.  exp/run.sh <슬러그> -m \"목적\" 으로 시작한다."; exit 0; }

printf "%-34s %-9s %8s %9s %9s %7s  %s\n" "RUN" "상태" "step" "sim[s]" "wall[s]" "출력" "목적"
printf "%.0s─" {1..120}; echo
for d in "${runs[@]}"; do
  d=${d%/}; id=$(basename "$d")
  st=$(run_state "$d")
  purpose=$(sed -n '/^## 1. 목적/,/^<!--/p' "$d/RUN.md" 2>/dev/null | sed '1d' | grep -v '^<!--' | sed '/^$/d' | head -1 | sed 's/\(.\{0,46\}\).*/\1/')
  printf "%-34s %-9s %8s %9s %9s %7s  %s\n" \
    "$id" "$st" "$(last_step "$d")" "$(last_simtime "$d")" "$(last_walltime "$d")" \
    "$(du -sh "$d/out" 2>/dev/null | cut -f1)" "$purpose"
done
echo
echo "자세히: exp/status.sh <run>   ·   감시: exp/watch.sh <run>   ·   결과기록: exp/report.sh <run>"
