#!/usr/bin/env bash
# exp/watch_temp.sh <run> — 온도 건전성 감시. 물리적으로 말이 안 되는 값이 나오면 종료(=알림).
#
# 안정 성층(T_con=285, T_ini=300)에서 온도는 [285, 300] 안에 있어야 한다.
# 지표보다 차갑거나 유입공기보다 뜨거우면 수치 문제다. 폭주는 ±1e21 로 간다.
set -uo pipefail
cd "$(dirname "$0")/.."
D="runs/$1"; LOG="$D/run.log"
LO=${2:-282.0}   # tt_min 하한 (지표 285 보다 3K 아래)
HI=${3:-305.0}   # tt_max 상한 (유입 300 보다 5K 위)
echo "watching temp $1  허용범위 tt_min>$LO  tt_max<$HI"
while true; do
  if [ -f "$LOG" ]; then
    L=$(grep -E "^[0-9]+[[:space:]]+tu_max" "$LOG" 2>/dev/null | tail -1)
    if [ -n "$L" ]; then
      R=$(echo "$L" | awk -v lo="$LO" -v hi="$HI" '{
        tmax=""; tmin="";
        for(i=1;i<=NF;i++){
          if($i ~ /^tt_max=/){ sub(/tt_max=/,"",$i); tmax=$i }
          if($i ~ /^tt_min=/){ sub(/tt_min=/,"",$i); tmin=$i }
        }
        if(tmax=="" || tmin=="") { print "SKIP"; exit }
        if(tmax ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || tmin ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/){ print "TRIP non-numeric"; exit }
        if(tmin+0 < lo+0){ printf "TRIP tt_min=%s < %s\n", tmin, lo; exit }
        if(tmax+0 > hi+0){ printf "TRIP tt_max=%s > %s\n", tmax, hi; exit }
        print "OK"
      }')
      case "$R" in
        TRIP*) echo "EVENT=TEMP_OUT_OF_RANGE  $1  $(date '+%F %T')"; echo "  $R"; echo "  $L"; exit 0 ;;
      esac
    fi
  fi
  # 프로세스가 끝나면 이 감시도 역할이 끝난다 (완주는 watch.sh 가 알린다)
  P=$(cat "$D/env/pid" 2>/dev/null) || true
  [ -n "${P:-}" ] && ! kill -0 "$P" 2>/dev/null && { echo "프로세스 종료 — 온도 감시 종료 (마지막: $L)"; exit 0; }
  sleep 300
done
