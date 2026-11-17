#!/usr/bin/env bash
# Run B 2건이 끝나기를 기다렸다가 270° A → 270° B 를 잇는다. 전부 GPU 4 (A100, sm_80).
# 중단: touch runs/chain_rot270.stop  (현재 런이 끝난 뒤 멈춘다)
set -u
cd "$(dirname "$0")/.."
GPU=4
LOG(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
STOP(){ [ -f runs/chain_rot270.stop ] && { LOG "중단 요청 — 체인 종료"; exit 0; }; }
COMMON="-D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2
        -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 -D LDM_DUMP_EVERY=100"

LOG "Run B 2건 완료 대기"
while ! grep -q "Run B 2건 완료" runs/chain_runB.log 2>/dev/null; do
  tmux has-session -t runB 2>/dev/null || { grep -q "Run B 2건 완료" runs/chain_runB.log 2>/dev/null || { LOG "!! runB 세션이 완료 표시 없이 사라짐 — 중단"; exit 1; }; }
  sleep 60
done
LOG "Run B 확인 — GPU 4 해제됨"
STOP

LOG "=== 270° 전이 LES (Run A) 시작 ==="
exp/run.sh leipzig-rot270 --fg -g $GPU -t ldm -i Leipzig_rot270 \
  -m "풍향 민감도 270°: 도시를 270° 회전시킨 중립 전이 LES. 건물(IBM 마커)만 도메인 중심 (657.5,0) 기준 회전하고 유체·지면·천장은 미변경. 평균창 150~299s 로 기준런 2026-09-15_leipzig-dt01-sponge-uy 및 90°·180° 와 동일. LDM 출력 100 step. GPU 4 (A100, sm_80). 목적: lambda_p-발자국 관계의 풍향 일반성 검증 (0/90/180/270 네 방향 완성)." \
  $COMMON -D AVG_START_TIME=150.0 -D AVG_END_TIME=299.0 -D CKPT_EVERY=2000 > /dev/null 2>&1
A270=$(ls -dt runs/*_leipzig-rot270 2>/dev/null | head -1)
LOG "270° Run A 종료 → $A270  마지막 step $(grep -oP '^\d+(?=\s+tu_max)' "$A270/run.log" 2>/dev/null | tail -1)"
STOP

if [ -s "$A270/out/mean_field.bin" ]; then
  LOG "=== 270° 평균장 LDM (Run B) 시작 ==="
  exp/run.sh leipzig-rot270-meanflow --fg -g $GPU -t ldm -i Leipzig_rot270 \
    -m "풍향 민감도 270° 평균장 LDM(Run B). $A270 의 mean_field.bin 을 얼려 LDM 만 재계산. GPU 4 (A100, sm_80) — 90°·180° Run B 와 같은 아키텍처." \
    $COMMON -D LDM_MEAN_MODE=1 -D MEAN_FIELD_PATH="\"./$A270/out/mean_field.bin\"" > /dev/null 2>&1
  LOG "270° Run B 종료 → $(ls -dt runs/*_leipzig-rot270-meanflow 2>/dev/null | head -1)"
else
  LOG "!! mean_field.bin 없음 — 270° Run B 취소"
fi
LOG "=== 270° 체인 완료 ==="
