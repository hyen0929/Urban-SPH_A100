#!/bin/bash
# exp/lib.sh — 실험 스크립트 공통 함수. 단독 실행하지 않는다.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNS_DIR="$REPO_ROOT/runs"
INDEX_MD="$RUNS_DIR/INDEX.md"

die()  { echo "error: $*" >&2; exit 1; }
info() { echo "  $*"; }

# runs/ 아래 이름 또는 경로를 실제 run 디렉토리로 해석.
#   resolve_run 2026-09-08_foo   → /abs/runs/2026-09-08_foo
#   resolve_run foo              → 부분일치 (유일할 때만)
#   resolve_run latest           → 가장 최근 생성된 run
resolve_run() {
  local q="$1" hits
  [ -n "$q" ] || die "run 이름이 필요하다. 목록: exp/status.sh"
  if [ "$q" = latest ]; then
    hits=$(ls -1dt "$RUNS_DIR"/*/ 2>/dev/null | head -1)
    [ -n "$hits" ] || die "runs/ 가 비어 있다"
    echo "${hits%/}"; return
  fi
  [ -d "$q" ] && { (cd "$q" && pwd); return; }
  [ -d "$RUNS_DIR/$q" ] && { echo "$RUNS_DIR/$q"; return; }
  hits=$(ls -1d "$RUNS_DIR"/*"$q"*/ 2>/dev/null)
  case $(echo "$hits" | grep -c .) in
    0) die "'$q' 에 맞는 run 이 없다. 목록: exp/status.sh" ;;
    1) echo "$(echo "$hits" | head -1 | sed 's:/$::')" ;;
    *) echo "'$q' 가 여러 run 에 일치한다:" >&2; echo "$hits" >&2; exit 1 ;;
  esac
}

# 실행 상태: RUNNING / DONE / DIED / NEVER_RAN / ARCHIVED(옛 기록 편입분)
run_state() {
  local d="$1" pid
  [ -d "$d/env" ] || { echo ARCHIVED; return; }
  [ -f "$d/env/pid" ] || { echo NEVER_RAN; return; }
  pid=$(cat "$d/env/pid")
  if kill -0 "$pid" 2>/dev/null; then echo RUNNING; return; fi
  if [ -f "$d/run.log" ] && grep -qiE "nan|segmentation|cudaError|illegal|core dumped" "$d/run.log"; then
    echo DIED; return
  fi
  # 마지막 진행 줄이 목표 시간에 도달했는지는 알 수 없으므로 로그 유무로만 판정
  [ -s "$d/run.log" ] && echo DONE || echo DIED
}

# run.log 마지막 진행 줄에서 값 추출: last_metric <dir> <key>
#   진행 줄 형식: "<step>\t tu_max= 8.76\ttt_max=300.00\t..."
last_metric() {
  local d="$1" key="$2"
  grep -oE "$key=[ ]*[-0-9.e+]+" "$d/run.log" 2>/dev/null | tail -1 | sed "s/$key=[ ]*//"
}
last_step() {
  grep -oE '^[0-9]+\s+tu_max' "$1/run.log" 2>/dev/null | tail -1 | awk '{print $1}'
}
last_simtime() {
  grep -oE 'Simulation time:[ ]*[0-9.]+' "$1/run.log" 2>/dev/null | tail -1 | awk '{print $NF}'
}
last_walltime() {
  grep -oE 'Computing time:[ ]*[0-9.]+' "$1/run.log" 2>/dev/null | tail -1 | awk '{print $NF}'
}
