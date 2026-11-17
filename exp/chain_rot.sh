#!/usr/bin/env bash
# 풍향 회전 실험 체인 — 각 각도마다 전이 LES(Run A) 완주 후 평균장 LDM(Run B)을 잇고,
# 다음 각도로 넘어간다. 90° 는 이미 실행 중이므로 그것이 끝나기를 기다리는 것부터 시작한다.
#
#   nohup bash exp/chain_rot.sh > runs/chain_rot.log 2>&1 &
#
# 중단: runs/chain_rot.stop 파일을 만들면 현재 런이 끝난 뒤 멈춘다.
set -u
cd "$(dirname "$0")/.."
GPU=4
LOG() { echo "[$(date '+%m-%d %H:%M')] $*"; }
STOP() { [ -f runs/chain_rot.stop ] && { LOG "중단 요청 확인 — 체인 종료"; exit 0; }; }

# 런이 끝날 때까지 대기. 성공이면 0, 실패면 1.
wait_run() {          # $1 = 런 디렉터리
  local d=$1 p
  for _ in $(seq 1 200); do [ -f "$d/env/pid" ] && break; sleep 5; done
  p=$(cat "$d/env/pid" 2>/dev/null) || return 1
  LOG "대기: $(basename "$d")  PID $p"
  while kill -0 "$p" 2>/dev/null; do sleep 60; done
  if grep -qEi 'error|abort|segmentation|failed' "$d/run.log" 2>/dev/null; then
    LOG "!! $(basename "$d") 에러 종료: $(grep -Ei 'error|abort|segmentation|failed' "$d/run.log" | tail -1)"
    return 1
  fi
  LOG "완료: $(basename "$d")"
  return 0
}

run_A() {             # $1 = 각도
  local a=$1 slug="leipzig-rot$a"
  LOG "=== ${a}° 전이 LES 시작 ==="
  exp/run.sh "$slug" \
    -m "풍향 민감도 ${a}°: 도시를 ${a}° 회전시킨 전이 LES. 건물(IBM 마커)만 도메인 중심 (657.5,0) 기준 회전하고 유체·지면·천장·경계는 미변경(도메인 정사각·격자 자기사상·지형 기복 없음). 평균창 150~299s 로 기준런 2026-09-15_leipzig-dt01-sponge-uy 와 동일. 목적: lambda_p-발자국 관계의 풍향 일반성 검증." \
    -i "Leipzig_rot$a" -g $GPU \
    -D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2 \
    -D SPONGE_DAMP_UY=1 -D AVG_START_TIME=150.0 -D AVG_END_TIME=299.0 \
    -D INIT_CROSSFLOW_UY=0.0 -D CKPT_EVERY=2000 \
    -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 >/dev/null 2>&1
  sleep 10
  wait_run "$(ls -dt runs/*_${slug} | head -1)"
}

run_B() {             # $1 = 각도, $2 = Run A 디렉터리
  local a=$1 ad=$2 slug="leipzig-rot$a-meanflow"
  if [ ! -s "$ad/out/mean_field.bin" ]; then
    LOG "!! $ad 에 mean_field.bin 이 없다 — ${a}° Run B 건너뜀"; return 1
  fi
  LOG "=== ${a}° 평균장 LDM 시작 ==="
  exp/run.sh "$slug" \
    -m "풍향 ${a}° 의 Run B(frozen mean-flow LDM). 같은 각도의 Run A 가 t=299s 에 덤프한 평균장으로 LDM 만 재계산. 방출·난수 동일, 차이는 분해된 비정상성의 유무뿐. 목적: 전이 대 평균장 차이가 풍향에 따라 달라지는가." \
    -i "Leipzig_rot$a" -g $GPU \
    -D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2 \
    -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0 -D LDM_MEAN_MODE=1 \
    -D "MEAN_FIELD_PATH=\"./$ad/out/mean_field.bin\"" \
    -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 >/dev/null 2>&1
  sleep 10
  wait_run "$(ls -dt runs/*_${slug} | head -1)"
}

LOG "체인 시작 — 90°(실행 중) → 180° → 270°, 각 각도마다 Run A 후 Run B"

# 90° 는 이미 돌고 있다
A90=$(ls -dt runs/*_leipzig-rot90 2>/dev/null | head -1)
if [ -n "$A90" ]; then
  wait_run "$A90" && run_B 90 "$A90"
else
  LOG "!! 실행 중인 90° 런을 찾지 못했다"
fi

for a in 180 270; do
  STOP
  run_A "$a" || { LOG "!! ${a}° Run A 실패 — 체인 중단"; exit 1; }
  ad=$(ls -dt runs/*_leipzig-rot$a | head -1)
  run_B "$a" "$ad"
done

LOG "체인 완료"
