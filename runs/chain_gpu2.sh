#!/usr/bin/env bash
# GPU 2 (H100 NVL, sm_90) 체인 — 180° A(실행중) → 180° B → 270° A → 270° B
# 중단: touch runs/chain_gpu2.stop  (현재 런이 끝난 뒤 멈춘다)
set -u
cd "$(dirname "$0")/.."
GPU=2; export ARCH=sm_90
LOG(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
STOP(){ [ -f runs/chain_gpu2.stop ] && { LOG "중단 요청 — 체인 종료"; exit 0; }; }
COMMON="-D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2
        -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 -D LDM_DUMP_EVERY=100"

wait_run(){ local d=$1 p
  for _ in $(seq 1 200); do [ -f "$d/env/pid" ] && break; sleep 5; done
  p=$(cat "$d/env/pid" 2>/dev/null) || return 1
  LOG "대기: $(basename "$d")  PID $p"
  while kill -0 "$p" 2>/dev/null; do sleep 60; done
  if grep -qEi 'cudaError|segmentation|abort|nan' "$d/run.log" 2>/dev/null; then
    LOG "!! $(basename "$d") 이상 종료: $(grep -Ei 'cudaError|segmentation|abort|nan' "$d/run.log" | tail -1)"; return 1; fi
  LOG "완료: $(basename "$d")"; return 0; }

run_A(){ local a=$1; local slug="leipzig-rot$a"
  LOG "=== ${a}° 전이 LES (Run A) 시작 ==="
  exp/run.sh "$slug" -g $GPU -t ldm -i "Leipzig_rot$a" \
    -m "풍향 민감도 ${a}°: 도시를 ${a}° 회전시킨 중립 전이 LES. 건물(IBM 마커)만 도메인 중심 (657.5,0) 기준 회전하고 유체·지면·천장은 미변경. 평균창 150~299s 로 기준런 2026-09-15_leipzig-dt01-sponge-uy 와 동일. LDM 출력 100 step. 이 런은 GPU 2 (H100 NVL, sm_90) 에서 돌아 A100 각도들과 반올림이 달라 난류 실현이 다르다. 목적: lambda_p-발자국 관계의 풍향 일반성 검증." \
    $COMMON -D AVG_START_TIME=150.0 -D AVG_END_TIME=299.0 -D CKPT_EVERY=2000 >/dev/null 2>&1
  sleep 10; wait_run "$(ls -dt runs/*_${slug} | head -1)"; }

run_B(){ local a=$1 ad=$2; local slug="leipzig-rot$a-meanflow"
  if [ ! -s "$ad/out/mean_field.bin" ]; then LOG "!! ${a}° mean_field.bin 없음 — Run B 취소"; return 1; fi
  LOG "=== ${a}° 평균장 LDM (Run B) 시작 ==="
  exp/run.sh "$slug" -g $GPU -t ldm -i "Leipzig_rot$a" \
    -m "풍향 민감도 ${a}° 평균장 LDM. ${ad} 의 mean_field.bin 을 얼려 LDM 만 재계산. 전이 대비 비정상성 기여 분리용." \
    $COMMON -D LDM_MEAN_MODE=1 -D MEAN_FIELD_PATH="\"./$ad/out/mean_field.bin\"" >/dev/null 2>&1
  sleep 10; wait_run "$(ls -dt runs/*_${slug} | head -1)"; }

A180=runs/2026-09-23_leipzig-rot180
LOG "=== GPU 2 체인 시작 (180° A 는 이미 실행 중) ==="
wait_run "$A180" && { STOP; run_B 180 "$A180"; }
STOP
if run_A 270; then A270=$(ls -dt runs/*_leipzig-rot270 | head -1); STOP; run_B 270 "$A270"; fi
LOG "=== GPU 2 체인 완료 ==="
