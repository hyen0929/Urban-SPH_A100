#!/bin/bash
# 감시견 — Run A(+체인으로 이어지는 Run B)를 5분마다 점검하고
# 사람이 봐야 할 것만 analysis/watchdog.log 에 ALERT/INFO 한 줄로 남긴다.
#   setsid nohup runs/2026-09-10_leipzig-full-bcfix_2/env/watchdog.sh >/dev/null 2>&1 &
set -u
cd "$(dirname "$0")/../../.." || exit 3

A="runs/2026-09-10_leipzig-full-bcfix_2"
B="runs/2026-09-10_leipzig-full-meanflow"
OUT="$A/analysis/watchdog.log"
CSV="$A/out/PPE_diag_explicit.csv"
MEAN="$A/out/mean_field.bin"
INTERVAL=300
STALL=900              # CSV(10 step 마다 fflush)가 이만큼 멈추면 정체
DISK_MIN_GB=200
declare -A SEEN

say(){ echo "[$(date '+%F %T')] $*" >> "$OUT"; }
once(){ local k="$1"; shift; [ -n "${SEEN[$k]:-}" ] && return 1; SEEN[$k]=1; say "$@"; return 0; }

APID=$(cat "$A/env/pid" 2>/dev/null)
say "INFO  watchdog 시작 (pid=$$) — Run A pid=$APID, ${INTERVAL}s 주기"

BASE_DIV=""; BASE_AVG=""; LAST_STEP=-1
while true; do
  NOW=$(date +%s)

  # ── Run A 진행 ──
  if [ -f "$CSV" ]; then
    L=$(tail -1 "$CSV")
    STEP=$(echo "$L" | cut -d, -f1); T=$(echo "$L" | cut -d, -f2)
    BI=$(echo "$L" | cut -d, -f3)
    # col5(div_un1_max, KGC 보정 후)는 t=1.5 s 부터 1.052148e+02 에 고정돼 있다 — 이상입자
    # 하나가 max 를 물고 있어 field 를 못 따라간다. 그래서 col7(보정 전 max)·col6(평균)을 본다.
    DIVP=$(echo "$L" | cut -d, -f7); DAVG=$(echo "$L" | cut -d, -f6)
    AGE=$(( NOW - $(stat -c %Y "$CSV") ))
    case "$BI$DIVP$DAVG" in *nan*|*NAN*|*inf*|*INF*) once nan "ALERT Run A NaN/Inf — step=$STEP t=$T bi_max=$BI div_plain_max=$DIVP div_avg=$DAVG";; esac
    if [ -z "$BASE_DIV" ] && [ -n "$T" ] && awk "BEGIN{exit !($T>=5)}" 2>/dev/null; then
      BASE_DIV="$DIVP"; BASE_AVG="$DAVG"
      say "INFO  기준 div_un1_plain_max=$BASE_DIV, div_un1_avg=$BASE_AVG (t=$T)"
    fi
    if [ -n "$BASE_DIV" ] && awk "BEGIN{exit !($DIVP > 3*$BASE_DIV || $DAVG > 3*$BASE_AVG)}" 2>/dev/null; then
      once div "ALERT Run A 발산 의심 — div_plain_max=$DIVP (기준 $BASE_DIV), div_avg=$DAVG (기준 $BASE_AVG), 3배 초과. step=$STEP t=$T"
    fi
    if kill -0 "$APID" 2>/dev/null && [ "$AGE" -gt "$STALL" ]; then
      once stall "ALERT Run A 정체 — CSV ${AGE}s 무갱신 (step=$STEP t=$T)"
    fi
    # 마일스톤
    if [ -n "$T" ] && awk "BEGIN{exit !($T>=150)}" 2>/dev/null; then once t150 "INFO  Run A t=150 s 도달 — 평균 누적/LDM 방출 시작 (step=$STEP)"; fi
    LAST_STEP="$STEP"
  fi

  # ── 로그의 치명 문자열 / 클립 포화 ──
  if [ -f "$A/run.log" ]; then
    if grep -aiqE "segmentation|cudaError|illegal memory|out of memory|core dumped|assert" "$A/run.log"; then
      once logerr "ALERT Run A 로그 에러: $(grep -aiE 'segmentation|cudaError|illegal memory|out of memory|core dumped|assert' "$A/run.log" | tail -1)"
    fi
    TU=$(grep -ao 'tu_max=[0-9.]*' "$A/run.log" | tail -1 | cut -d= -f2)
    if [ -n "${TU:-}" ] && awk "BEGIN{exit !($TU>=19.0)}" 2>/dev/null; then
      once clip "ALERT Run A 속도클립 포화 임박 — tu_max=$TU (limit 20 m/s), step≈$LAST_STEP"
    fi
  fi

  # ── mean_field.bin ──
  if [ -f "$MEAN" ]; then once mean "INFO  mean_field.bin 생성 ($(stat -c %s "$MEAN") B)"; fi

  # ── Run A 종료 ──
  if [ -n "$APID" ] && ! kill -0 "$APID" 2>/dev/null; then
    if once adead "INFO  Run A 종료 — 마지막 step=$LAST_STEP, mean_field.bin=$([ -f "$MEAN" ] && stat -c %s "$MEAN" || echo 없음)"; then
      [ -f "$MEAN" ] || say "ALERT Run A 가 mean_field.bin 없이 종료했다 — Run B 기동 불가"
    fi
  fi

  # ── Run B ──
  BPID=$(cat "$B/env/pid" 2>/dev/null || true)
  if [ -n "${BPID:-}" ]; then
    if kill -0 "$BPID" 2>/dev/null; then
      once bstart "INFO  Run B 기동 확인 (pid=$BPID)"
    else
      once bdead "INFO  Run B 종료 (pid=$BPID) — $(tail -c 200 "$B/run.log" 2>/dev/null | tr '\n' ' ')"
    fi
  fi
  grep -aq FATAL "$B/env/chain.log" 2>/dev/null && once chainfatal "ALERT 체인 실패: $(grep -a FATAL "$B/env/chain.log" | tail -1)"

  # ── 디스크 ──
  FREE=$(df -BG --output=avail /home | tail -1 | tr -dc '0-9')
  [ -n "$FREE" ] && [ "$FREE" -lt "$DISK_MIN_GB" ] && once disk "ALERT 디스크 여유 ${FREE} GB < ${DISK_MIN_GB} GB"

  # 둘 다 끝났으면 감시 종료
  if [ -n "${SEEN[adead]:-}" ] && [ -n "${SEEN[bdead]:-}" ]; then say "INFO  두 런 모두 종료 — watchdog 종료"; exit 0; fi
  sleep "$INTERVAL"
done
