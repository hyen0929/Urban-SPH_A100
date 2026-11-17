#!/usr/bin/env bash
# exp/swap_localust.sh — 안정 성층 런을 국소 u* 바이너리로 교체한다.
#
# step >= WANT 인 체크포인트가 기록되기를 기다렸다가, 본 런을 세우고
# 새 바이너리로 그 체크포인트에서 재개한다. 실패하면 옛 바이너리를 되돌린다.
#
#   nohup bash exp/swap_localust.sh > runs/swap_localust.log 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/.."
D=runs/2026-09-21_leipzig-most-stable
SCR=/tmp/claude-1014/-home-kdh-Urban-SPH/90fe7215-f210-4eca-bb99-e4c2226c7eef/scratchpad
NEWBIN=$SCR/SOPHIA_gpu_localust
TARGET=$D/env/SOPHIA_gpu
GPU=4
WANT=${1:-6000}
LOG(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

ckpt_step(){ python3 - <<PY 2>/dev/null || echo -1
import struct
try:
    b=open("$D/out/restart.bin","rb").read(52)
    print(struct.unpack("<i", b[48:52])[0] if len(b)==52 else -1)
except Exception:
    print(-1)
PY
}

MARK=runs/.localust_regression_ok
[ -x "$NEWBIN" ] || { LOG "!! 새 바이너리 없음: $NEWBIN"; exit 1; }
PID=$(cat "$D/env/pid" 2>/dev/null) || { LOG "!! env/pid 없음"; exit 1; }
LOG "대기 시작 — 체크포인트 step >= $WANT (본 런 PID $PID)"

while :; do
  kill -0 "$PID" 2>/dev/null || { LOG "!! 본 런이 먼저 종료됨 — 교체 중단"; exit 1; }
  S=$(ckpt_step)
  case "$S" in ''|*[!0-9-]*) S=-1 ;; esac
  if [ "$S" -ge "$WANT" ] && [ -f "$MARK" ]; then LOG "체크포인트 step=$S + 회귀 통과 확인 — 교체 진행"; break; fi
  sleep 60
done

# 1) Run B 체인을 먼저 세운다 (본 런 종료를 완주로 오인하지 않도록)
touch runs/chain_stable.stop
pkill -f "chain_stable.sh" 2>/dev/null && LOG "체인 정지" || LOG "체인 이미 없음"

# 2) 본 런 종료
LOG "본 런 종료 (PID $PID)"
kill "$PID" 2>/dev/null
for _ in $(seq 1 60); do kill -0 "$PID" 2>/dev/null || break; sleep 2; done
kill -0 "$PID" 2>/dev/null && { LOG "!! 종료 실패 — 중단"; exit 1; }
LOG "종료 확인"

# 3) 바이너리 교체 + 재현정보 보존
cp -p "$TARGET" "$D/env/SOPHIA_gpu.pre-localust" && LOG "옛 바이너리 보존: env/SOPHIA_gpu.pre-localust"
cp -p "$NEWBIN" "$TARGET" || { LOG "!! 바이너리 복사 실패"; exit 1; }
git diff -- Parameters.cuh functions_MOST.cuh ISPH_Calc.cuh functions_DD.cuh > "$D/env/localust.diff" 2>/dev/null
LOG "새 바이너리 설치 + env/localust.diff 저장"

# 4) 재개
CKPT=$D/out/restart.bin
LOG "재개: SOPHIA_RESTART=$CKPT  GPU=$GPU"
cd "$(pwd)"
nohup env CUDA_VISIBLE_DEVICES=$GPU SOPHIA_RESTART="$CKPT" \
      stdbuf -oL -eL "$TARGET" 1 >> "$D/run.log" 2>&1 &
NEWPID=$!
echo "$NEWPID" > "$D/env/pid"
LOG "새 PID $NEWPID (env/pid 갱신)"

# 5) 재개 확인
for _ in $(seq 1 120); do
  grep -q "ckpt\] RESTART" "$D/run.log" && { LOG "재개 확인: $(grep 'ckpt\] RESTART' "$D/run.log" | tail -1)"; break; }
  kill -0 "$NEWPID" 2>/dev/null || { LOG "!! 새 프로세스가 죽었다"; tail -20 "$D/run.log"; exit 1; }
  sleep 10
done

# 6) Run B 체인 재기동
rm -f runs/chain_stable.stop
nohup bash exp/chain_stable.sh >> runs/chain_stable.log 2>&1 &
LOG "Run B 체인 재기동"
LOG "교체 완료"
