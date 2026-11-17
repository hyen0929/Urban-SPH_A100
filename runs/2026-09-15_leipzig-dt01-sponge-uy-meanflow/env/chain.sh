#!/bin/bash
# Run B 자동 기동 체인 (2026-09-14)
#
# Run A(2026-09-15_leipzig-dt01-sponge-uy)가 t=299 s 에 mean_field.bin 을 떨구고 t=300 s 에 끝나면,
# 그 파일이 온전한지 확인한 뒤 Run B 를 같은 GPU 에 띄운다.
#
#   setsid nohup runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow/env/chain.sh >/dev/null 2>&1 &
#
# 2026-09-10 판과 달리 Run A 가 mean_field.bin 없이 죽어도 끝내지 않는다 — resume.sh 로
# 재개되면 env/pid 가 바뀌므로 매 주기 pid 를 다시 읽고 계속 기다린다.
# 종료 코드: 0 = Run B 기동, 1 = 파일 검증/기동 실패, 2 = 타임아웃
set -u
cd "$(dirname "$0")/../../.." || exit 3

A_DIR="runs/2026-09-15_leipzig-dt01-sponge-uy"
B_DIR="runs/2026-09-15_leipzig-dt01-sponge-uy-meanflow"
MEAN="$A_DIR/out/mean_field.bin"
LOG="$B_DIR/env/chain.log"
GPU=4
MIN_SIZE=800000000          # 861 MB 예상. 이보다 작으면 미완성으로 본다
INTERVAL=60
DEADLINE=$(( $(date +%s) + 14*24*3600 ))   # 14 일

log(){ echo "[$(date '+%F %T')] $*" >> "$LOG"; }
log "chain 시작 (pid=$$). A=$A_DIR  B=$B_DIR  GPU=$GPU"

# ── 1. Run A 정상 종료 대기 (재개 추적) ──────────────────────────────
WARNED=""
while true; do
  if [ "$(date +%s)" -gt "$DEADLINE" ]; then log "FATAL: 14 일 타임아웃 — Run B 기동 안 함"; exit 2; fi
  APID=$(cat "$A_DIR/env/pid" 2>/dev/null)
  if [ -n "$APID" ] && kill -0 "$APID" 2>/dev/null; then WARNED=""; sleep "$INTERVAL"; continue; fi
  if [ -f "$MEAN" ] && grep -aqE "^Simulation time: (299\.[5-9]|300\.)" "$A_DIR/run.log"; then break; fi
  if [ "$WARNED" != "$APID" ]; then
    log "Run A(pid=${APID:-?}) 가 mean_field.bin 없이 멈춰 있음 — resume.sh 재개를 기다린다"
    WARNED="$APID"
  fi
  sleep "$INTERVAL"
done
log "Run A 정상 종료 감지 (pid=$APID)"

# ── 2. mean_field.bin 검증 (크기 안정·크기·magic) ─────────────────────
S1=$(stat -c %s "$MEAN"); sleep 30; S2=$(stat -c %s "$MEAN")
while [ "$S1" != "$S2" ]; do log "쓰기 진행 중($S1 -> $S2)"; S1=$S2; sleep 30; S2=$(stat -c %s "$MEAN"); done
if [ "$S2" -lt "$MIN_SIZE" ]; then log "FATAL: $MEAN 크기 $S2 B < $MIN_SIZE B — 잘린 파일. Run B 기동 안 함"; exit 1; fi
MAGIC=$(head -c 14 "$MEAN")
if [ "$MAGIC" != "SOPHIA_MEAN_v1" ]; then log "FATAL: magic 불일치 ('$MAGIC'). Run B 기동 안 함"; exit 1; fi
log "mean_field.bin OK ($S2 B, magic=$MAGIC)"

# ── 3. GPU 비었는지 확인 ──────────────────────────────────────────────
for i in $(seq 1 20); do
  USED=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$GPU" 2>/dev/null | tr -d ' ')
  [ -n "$USED" ] && [ "$USED" -lt 2000 ] && break
  log "GPU $GPU 아직 ${USED:-?} MiB 사용 중 — 30 s 후 재확인 ($i/20)"; sleep 30
done
if [ -z "$USED" ] || [ "$USED" -ge 2000 ]; then log "FATAL: GPU $GPU 가 10 분째 사용 중 — Run B 기동 안 함"; exit 1; fi

# ── 4. Run B 기동 ─────────────────────────────────────────────────────
CUDA_VISIBLE_DEVICES=$GPU setsid nohup stdbuf -oL "$B_DIR/env/SOPHIA_gpu" 1 > "$B_DIR/run.log" 2>&1 &
BPID=$!
echo "$BPID" > "$B_DIR/env/pid"
sleep 10
if kill -0 "$BPID" 2>/dev/null; then
  log "Run B 기동 완료 (pid=$BPID, GPU=$GPU) → $B_DIR/run.log"; exit 0
else
  log "FATAL: Run B 가 기동 직후 죽었다. 로그 마지막:"; tail -20 "$B_DIR/run.log" >> "$LOG"; exit 1
fi
