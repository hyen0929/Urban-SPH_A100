#!/usr/bin/env bash
# 회전 체인 감시. 이상이 있거나 전체가 끝났을 때만 종료(=알림)한다.
set -uo pipefail
cd "$(dirname "$0")/.."
INT=300; STALL=2400            # 5분 주기, 40분 무변화면 정체로 본다
say(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
active(){ ls -dt runs/2026-09-2*_leipzig-rot*/ 2>/dev/null | head -1; }

while true; do
  if grep -q "270° 체인 완료" runs/chain_rot270.log 2>/dev/null; then
    say "EVENT=ALL_DONE  회전 체인 전체 완료"; tail -6 runs/chain_rot270.log; exit 0; fi
  if ! tmux has-session -t rot270 2>/dev/null; then
    say "EVENT=CHAIN_GONE  rot270 tmux 세션이 완료 표시 없이 사라짐"; tail -6 runs/chain_rot270.log 2>/dev/null; exit 0; fi

  d=$(active); log="$d/run.log"
  if [ -f "$log" ]; then
    if grep -qiE 'cudaerror|segmentation|out of memory|abort' "$log"; then
      say "EVENT=ERROR  $(basename $d)"; grep -iE 'cudaerror|segmentation|out of memory|abort' "$log" | tail -3; exit 0; fi
    if grep -qE 'tu_max=(nan|-?inf)|uz_max=(nan|-?inf)' "$log"; then
      say "EVENT=NAN  $(basename $d)"; tail -3 "$log"; exit 0; fi
    age=$(( $(date +%s) - $(stat -c %Y "$log") ))
    # 런 사이 전환(빌드·입력읽기)에는 로그가 길게 조용하므로, 계산이 도는 중일 때만 정체로 본다
    if [ $age -gt $STALL ] && pgrep -f "$d/env/SOPHIA_gpu" >/dev/null 2>&1; then
      say "EVENT=STALL  $(basename $d)  로그 ${age}s 무변화"; tail -3 "$log"; exit 0; fi
  fi
  sleep $INT
done
