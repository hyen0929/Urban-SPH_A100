#!/bin/bash
# exp/watch.sh — 실험 건전성 감시. 손이 필요할 때만 종료(=알림)한다.
# watch_implicit.sh 를 모든 run 에 대해 일반화한 것.
#
#   exp/watch.sh <run> [--stall SECS] [--interval SECS]
#
# 종료 조건: 프로세스 종료 / 로그의 NaN·error / 진행 정체(기본 20분).
# 그 외에는 매 주기 heartbeat 한 줄만 찍는다.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

RUN=""; STALL=1200; INTERVAL=60
while [ $# -gt 0 ]; do
  case "$1" in
    --stall)    STALL="$2"; shift 2 ;;
    --interval) INTERVAL="$2"; shift 2 ;;
    -h|--help)  sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 0 ;;
    *)          RUN="$1"; shift ;;
  esac
done
D=$(resolve_run "${RUN:-latest}") || exit 1
ID=$(basename "$D"); LOG="$D/run.log"
PID=$(cat "$D/env/pid" 2>/dev/null) || die "$ID: env/pid 가 없다 (아직 실행 안 됨)"

echo "watching $ID  pid=$PID  stall=${STALL}s"
while true; do
  if ! kill -0 "$PID" 2>/dev/null; then
    echo "EVENT=PROCESS_ENDED  $ID  pid $PID  $(date '+%F %T')"
    tail -20 "$LOG"; exit 0
  fi
  if grep -iqE "nan|segmentation|cudaError|illegal|assert|core dumped|out of memory" "$LOG" 2>/dev/null; then
    echo "EVENT=ERROR_IN_LOG  $ID  $(date '+%F %T')"
    grep -inE "nan|segmentation|cudaError|illegal|assert|core dumped|out of memory" "$LOG" | tail -8
    exit 0
  fi
  AGE=$(( $(date +%s) - $(stat -c %Y "$LOG" 2>/dev/null || date +%s) ))
  if [ "$AGE" -gt "$STALL" ]; then
    echo "EVENT=STALLED  $ID  로그 ${AGE}s 무변화  $(date '+%F %T')  step=$(last_step "$D")"
    exit 0
  fi
  echo "[$(date '+%T')] ok  $ID  step=$(last_step "$D")  sim=$(last_simtime "$D")s  tu_max=$(last_metric "$D" tu_max)  tt_max=$(last_metric "$D" tt_max)"
  sleep "$INTERVAL"
done
