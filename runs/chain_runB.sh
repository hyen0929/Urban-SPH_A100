#!/usr/bin/env bash
# 90°·180° 평균장 LDM (Run B) 순차 실행. GPU 4 (A100, sm_80).
# 두 각도의 Run B 를 같은 아키텍처에서 돌려 각도 간 비교를 깨끗하게 유지한다.
set -u
cd "$(dirname "$0")/.."
GPU=4
LOG(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
COMMON="-D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2
        -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 -D LDM_DUMP_EVERY=100"

run_B(){ local a=$1 ad=$2; local slug="leipzig-rot$a-meanflow"
  if [ ! -s "$ad/out/mean_field.bin" ]; then LOG "!! ${a}° mean_field.bin 없음 — 취소"; return 1; fi
  LOG "=== ${a}° 평균장 LDM 시작 ==="
  exp/run.sh "$slug" --fg -g $GPU -t ldm -i "Leipzig_rot$a" \
    -m "풍향 민감도 ${a}° 평균장 LDM(Run B). ${ad} 의 mean_field.bin 을 얼려 LDM 만 재계산하여 전이 대비 비정상성 기여를 분리한다. LDM 출력 100 step. GPU 4 (A100, sm_80) — 180° Run A 는 H100 이었으나 Run B 는 두 각도 모두 A100 으로 맞춰 각도 간 비교를 일치시켰다." \
    $COMMON -D LDM_MEAN_MODE=1 -D MEAN_FIELD_PATH="\"./$ad/out/mean_field.bin\"" > /dev/null 2>&1
  local rc=$?
  LOG "${a}° Run B 종료 rc=$rc  →  $(ls -dt runs/*_${slug} 2>/dev/null | head -1)"
  return $rc; }

run_B 90  runs/2026-09-23_leipzig-rot90
run_B 180 runs/2026-09-23_leipzig-rot180
LOG "=== Run B 2건 완료 ==="
