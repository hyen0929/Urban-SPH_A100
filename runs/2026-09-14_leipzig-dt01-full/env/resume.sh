#!/bin/bash
# Run A 재개 — 정전·재부팅·비정상 종료 뒤 이것 하나만 실행하면 된다 (리포 어디서든).
#
#   runs/2026-09-14_leipzig-dt01-full/env/resume.sh [체크포인트 경로]
#
# 1. Run A 가 이미 살아 있으면 아무것도 안 함
# 2. out/restart.bin 이 온전한지(끝 magic) 확인, 깨졌으면 out/restart.bin.prev 사용
# 3. GPU 4 가 비었는지 확인 (남의 작업이 있으면 중단하고 알림)
# 4. SOPHIA_RESTART=<체크포인트> 로 같은 바이너리 재기동, run.log 에 이어 씀
# 5. 감시견(watchdog.sh, groundcheck.py)·Run B 체인(chain.sh)이 죽어 있으면 다시 띄움
set -u
cd "$(dirname "$0")/../../.." || exit 3

A="runs/2026-09-14_leipzig-dt01-full"
B="runs/2026-09-14_leipzig-dt01-meanflow"
GPU=4

say(){ echo "[$(date '+%F %T')] $*" | tee -a "$A/analysis/watchdog.log"; }

APID=$(cat "$A/env/pid" 2>/dev/null)
if [ -n "$APID" ] && kill -0 "$APID" 2>/dev/null; then
  echo "Run A 가 이미 실행 중 (pid=$APID) — 재개하지 않음"; exit 1
fi
if grep -aq "\[mean\] saved" "$A/run.log" 2>/dev/null && [ -f "$A/out/mean_field.bin" ] && \
   grep -aqE "^Simulation time: (299\.[5-9]|300\.)" "$A/run.log"; then
  echo "Run A 는 이미 끝났다 (mean_field.bin 있음, t≈300 s) — 재개하지 않음"; exit 1
fi

intact(){ [ -f "$1" ] && [ "$(tail -c 16 "$1" | head -c 15)" = "SOPHIA_CKPT_END" ]; }
CK="${1:-}"
if [ -z "$CK" ]; then
  if intact "$A/out/restart.bin"; then CK="$A/out/restart.bin"
  elif intact "$A/out/restart.bin.prev"; then CK="$A/out/restart.bin.prev"; say "WARN  restart.bin 이 깨져 restart.bin.prev 로 재개"
  else echo "온전한 체크포인트가 없다 ($A/out/restart.bin{,.prev}) — 처음부터 다시 돌려야 함"; exit 1
  fi
elif ! intact "$CK"; then
  echo "$CK 가 온전하지 않다 (끝 magic 없음)"; exit 1
fi

USED=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$GPU" 2>/dev/null | tr -d ' ')
if [ -z "$USED" ] || [ "$USED" -ge 2000 ]; then
  echo "GPU $GPU 사용 중 (${USED:-?} MiB) — 누구 작업인지 확인 필요:"
  nvidia-smi --query-compute-apps=pid,gpu_uuid,used_memory --format=csv,noheader
  exit 1
fi

echo "" >> "$A/run.log"
echo "=== RESTART $(date '+%F %T') from $CK ===" >> "$A/run.log"
SOPHIA_RESTART="$CK" CUDA_VISIBLE_DEVICES=$GPU setsid nohup stdbuf -oL "$A/env/SOPHIA_gpu" 1 >> "$A/run.log" 2>&1 &
NPID=$!
echo "$NPID" > "$A/env/pid"
say "INFO  Run A 재개 (pid=$NPID, GPU=$GPU, from $CK)"

pgrep -f "$A/env/watchdog.sh" >/dev/null || { setsid nohup "$A/env/watchdog.sh" >/dev/null 2>&1 & say "INFO  watchdog 재기동"; }
pgrep -f "$A/env/groundcheck.py" >/dev/null || { setsid nohup "$A/env/groundcheck.py" >/dev/null 2>&1 & say "INFO  groundcheck 재기동"; }
pgrep -f "$B/env/chain.sh"    >/dev/null || { setsid nohup "$B/env/chain.sh"    >/dev/null 2>&1 & say "INFO  Run B 체인 재기동"; }

echo "입력 읽기+체크포인트 로드에 수 분 걸린다. 확인: grep -a '\[ckpt\]' $A/run.log"
