#!/bin/bash
# 감시견 — Run A(+체인으로 이어지는 Run B)를 5분마다 점검하고
# 사람이 봐야 할 것만 analysis/watchdog.log 에 ALERT/WARN/INFO 한 줄로 남긴다.
#   setsid nohup runs/2026-09-15_leipzig-dt01-sponge-uy/env/watchdog.sh >/dev/null 2>&1 &
# 2026-09-10 판에서 바뀐 점: pid 를 매 주기 다시 읽어 resume.sh 재개를 따라가고, 체크포인트 실패와
# Run A(정전 전, dt=0.005) 대비 이탈(analysis/track_vs_runA.py)을 본다. tu_max 클립 경보는 뺐다
# (t≈5 s 부터 20 m/s 포화가 알려진 상태라 신호가 아니다).
set -u
cd "$(dirname "$0")/../../.." || exit 3

A="runs/2026-09-15_leipzig-dt01-sponge-uy"
B="runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow"
OUT="$A/analysis/watchdog.log"
CSV="$A/out/PPE_diag_explicit.csv"
MEAN="$A/out/mean_field.bin"
INTERVAL=300
STALL=900              # CSV(10 step 마다 fflush, ~48 s)가 이만큼 멈추면 정체
DISK_MIN_GB=200
declare -A SEEN

say(){ echo "[$(date '+%F %T')] $*" >> "$OUT"; }
once(){ local k="$1"; shift; [ -n "${SEEN[$k]:-}" ] && return 1; SEEN[$k]=1; say "$@"; return 0; }

say "INFO  watchdog 시작 (pid=$$), ${INTERVAL}s 주기"
LAST_TRACK=""; LAST_CKPT=""; A_DONE=0
while true; do
  NOW=$(date +%s)
  APID=$(cat "$A/env/pid" 2>/dev/null)
  ALIVE=0; ETIME=0
  if [ -n "$APID" ] && kill -0 "$APID" 2>/dev/null; then ALIVE=1; ETIME=$(ps -o etimes= -p "$APID" | tr -d ' '); fi

  # ── Run A 진행 ──
  if [ -f "$CSV" ]; then
    L=$(tail -1 "$CSV")
    STEP=$(echo "$L" | cut -d, -f1); T=$(echo "$L" | cut -d, -f2)
    case "$L" in *nan*|*NAN*|*inf*|*INF*) once nan "ALERT Run A NaN/Inf — $L";; esac
    AGE=$(( NOW - $(stat -c %Y "$CSV") ))
    # 재개 직후엔 입력 읽기·로드로 수 분간 CSV 가 안 움직이므로 기동 20 분 뒤부터만 본다
    if [ "$ALIVE" = 1 ] && [ "${ETIME:-0}" -gt 1200 ] && [ "$AGE" -gt "$STALL" ]; then
      once "stall$APID" "ALERT Run A 정체 — CSV ${AGE}s 무갱신 (step=$STEP t=$T, pid=$APID)"
    fi
    for M in 50 100 150 200 237 250 299; do
      if [ -n "$T" ] && awk "BEGIN{exit !($T>=$M)}" 2>/dev/null; then once "t$M" "INFO  Run A t=$M s 도달 (step=$STEP)"; fi
    done
  fi

  # ── Run A 대 정전 전 Run A 이탈 ──
  TR=$(python3 "$A/analysis/track_vs_runA.py" --check 2>&1 | tail -1)
  LV=${TR%% *}
  if [ "$LV" != "$LAST_TRACK" ] && { [ "$LV" = WARN ] || [ "$LV" = FAIL ] || [ -n "$LAST_TRACK" ]; }; then
    case "$LV" in
      FAIL) say "ALERT 이탈 판정 FAIL — $TR";;
      WARN) say "WARN  이탈 경고 — $TR";;
      OK)   say "INFO  이탈 판정 복귀 — $TR";;
    esac
  fi
  LAST_TRACK="$LV"

  # ── 로그: 치명 문자열 / 체크포인트 ──
  if [ -f "$A/run.log" ]; then
    E=$(grep -aiE "segmentation|cudaError|illegal memory|out of memory|core dumped|assert|\[ckpt\] (ERROR|FATAL)" "$A/run.log" | tail -1)
    [ -n "$E" ] && once "err:$E" "ALERT Run A 로그 에러: $E"
    C=$(grep -a "\[ckpt\] saved" "$A/run.log" | tail -1)
    if [ -n "$C" ] && [ "$C" != "$LAST_CKPT" ]; then LAST_CKPT="$C"; fi
  fi

  # ── mean_field.bin / Run A 종료 ──
  [ -f "$MEAN" ] && once mean "INFO  mean_field.bin 생성 ($(stat -c %s "$MEAN") B)"
  if [ "$ALIVE" = 0 ] && [ -n "$APID" ]; then
    if [ -f "$MEAN" ] && grep -aqE "^Simulation time: (299\.[5-9]|300\.)" "$A/run.log"; then
      once adone "INFO  Run A 정상 종료 — mean_field.bin $(stat -c %s "$MEAN") B" && A_DONE=1
      A_DONE=1
    else
      once "adead$APID" "ALERT Run A 비정상 종료 (pid=$APID, 마지막 t=${T:-?}) — 마지막 체크포인트: ${LAST_CKPT:-없음}. 재개: $A/env/resume.sh"
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

  if [ "$A_DONE" = 1 ] && [ -n "${SEEN[bdead]:-}" ]; then say "INFO  두 런 모두 종료 — watchdog 종료"; exit 0; fi
  sleep "$INTERVAL"
done
