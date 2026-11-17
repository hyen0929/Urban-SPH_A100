#!/usr/bin/env bash
# 3창 앙상블 런 감시. 이상이 있거나 완주했을 때만 종료(=알림).
set -uo pipefail
cd "$(dirname "$0")/.."
D=runs/2026-09-30_leipzig-ens3-src29
INT=300; STALL=3600          # 5분 주기, 60분 무변화면 정체 (4KB 블록버퍼 감안)
say(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
while true; do
  if ! pgrep -f "$D/env/SOPHIA_gpu" >/dev/null 2>&1; then
    last=$(grep -oP '^\d+(?=\s+tu_max)' $D/run.log 2>/dev/null | tail -1)
    if [ "${last:-0}" -ge 74000 ] 2>/dev/null; then say "EVENT=DONE  완주 (step $last)"; else say "EVENT=DIED  프로세스 종료 (step ${last:-?}/75000)"; fi
    tail -5 $D/run.log; exit 0; fi
  if grep -qiE 'cudaerror|segmentation|out of memory|abort' $D/run.log 2>/dev/null; then
    say "EVENT=ERROR"; grep -iE 'cudaerror|segmentation|out of memory|abort' $D/run.log | tail -3; exit 0; fi
  if grep -qE 'tu_max=(nan|-?inf)|uz_max=(nan|-?inf)' $D/run.log 2>/dev/null; then
    say "EVENT=NAN"; tail -3 $D/run.log; exit 0; fi
  age=$(( $(date +%s) - $(stat -c %Y $D/run.log) ))
  if [ $age -gt $STALL ]; then say "EVENT=STALL  로그 ${age}s 무변화"; tail -3 $D/run.log; exit 0; fi
  df=$(df --output=avail -BG /home | tail -1 | tr -dc 0-9)
  if [ "$df" -lt 100 ]; then say "EVENT=DISK  /home 여유 ${df} GB"; exit 0; fi
  sleep $INT
done
