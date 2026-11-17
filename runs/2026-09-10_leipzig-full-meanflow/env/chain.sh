#!/bin/bash
# Run B 자동 기동 체인 (2026-09-10)
#
# Run A(2026-09-10_leipzig-full-bcfix_2)가 t=299 s 에 mean_field.bin 을 떨구고
# t=300 s 에 종료하면, 그 파일이 온전한지 확인한 뒤 Run B 를 같은 GPU 에 띄운다.
#
#   setsid nohup runs/2026-09-10_leipzig-full-meanflow/env/chain.sh >/dev/null 2>&1 &
#
# 종료 코드: 0 = Run B 기동, 1 = Run A 가 파일 없이 죽음(기동 안 함), 2 = 타임아웃
set -u
cd "$(dirname "$0")/../../.." || exit 3

A_DIR="runs/2026-09-10_leipzig-full-bcfix_2"
B_DIR="runs/2026-09-10_leipzig-full-meanflow"
MEAN="$A_DIR/out/mean_field.bin"
LOG="$B_DIR/env/chain.log"
GPU=4
MIN_SIZE=800000000          # 861 MB 예상. 이보다 작으면 미완성으로 본다
INTERVAL=60
DEADLINE=$(( $(date +%s) + 120*3600 ))   # 120 h

log(){ echo "[$(date '+%F %T')] $*" >> "$LOG"; }

log "chain 시작 (pid=$$). A=$A_DIR  B=$B_DIR  GPU=$GPU"

APID=$(cat "$A_DIR/env/pid" 2>/dev/null)
if [ -z "$APID" ]; then log "FATAL: Run A 의 env/pid 가 없다"; exit 1; fi
log "Run A pid=$APID 대기 시작"

# ── 1. Run A 종료 대기 ────────────────────────────────────────────────
while kill -0 "$APID" 2>/dev/null; do
  if [ "$(date +%s)" -gt "$DEADLINE" ]; then
    log "FATAL: 120 h 타임아웃 — Run A 가 아직 살아있다. Run B 기동 안 함"; exit 2
  fi
  sleep "$INTERVAL"
done
log "Run A 종료 감지"

# ── 2. mean_field.bin 검증 (존재·크기·안정) ───────────────────────────
if [ ! -f "$MEAN" ]; then
  log "FATAL: Run A 가 끝났는데 $MEAN 이 없다. 덤프가 안 된 것 — Run B 기동 안 함"
  log "  확인: grep -a '\[mean\]' $A_DIR/run.log   / AVG_END_TIME < end-time 인지"
  exit 1
fi
S1=$(stat -c %s "$MEAN"); sleep 30; S2=$(stat -c %s "$MEAN")
if [ "$S1" != "$S2" ]; then
  log "쓰기 진행 중($S1 -> $S2). 안정될 때까지 대기"
  while S1=$(stat -c %s "$MEAN"); sleep 30; S2=$(stat -c %s "$MEAN"); [ "$S1" != "$S2" ]; do :; done
fi
if [ "$S2" -lt "$MIN_SIZE" ]; then
  log "FATAL: $MEAN 크기 $S2 B < $MIN_SIZE B — 잘린 파일. Run B 기동 안 함"; exit 1
fi
MAGIC=$(head -c 14 "$MEAN")
if [ "$MAGIC" != "SOPHIA_MEAN_v1" ]; then
  log "FATAL: magic 불일치 ('$MAGIC'). Run B 기동 안 함"; exit 1
fi
log "mean_field.bin OK ($S2 B, magic=$MAGIC)"

# ── 3. GPU 비었는지 확인 ──────────────────────────────────────────────
for i in $(seq 1 20); do
  USED=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$GPU" 2>/dev/null | tr -d ' ')
  [ -n "$USED" ] && [ "$USED" -lt 2000 ] && break
  log "GPU $GPU 아직 ${USED:-?} MiB 사용 중 — 30 s 후 재확인 ($i/20)"; sleep 30
done

# ── 4. Run B 기동 ─────────────────────────────────────────────────────
CUDA_VISIBLE_DEVICES=$GPU setsid nohup "$B_DIR/env/SOPHIA_gpu" 1 > "$B_DIR/run.log" 2>&1 &
BPID=$!
echo "$BPID" > "$B_DIR/env/pid"
sleep 10
if kill -0 "$BPID" 2>/dev/null; then
  log "Run B 기동 완료 (pid=$BPID, GPU=$GPU) → $B_DIR/run.log"
  exit 0
else
  log "FATAL: Run B 가 기동 직후 죽었다. 로그 마지막:"; tail -20 "$B_DIR/run.log" >> "$LOG"
  exit 1
fi
