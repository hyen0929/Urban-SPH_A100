#!/usr/bin/env bash
# 50 소스 × 4 풍향, 평균장 LDM + 벽면침적 ON (DD_WALL_DEP=1).
# 기존 2026-10-06_leipzig-src50-rot* (벽면 off) 와 별도 경로(slug src50w)에 저장한다.
# solv 에 buffer_type 출력을 넣어 0=공중 / 1=지면침적 / 2=벽면·지붕침적 을 직접 구분한다.
# 중단: touch runs/chain_src50w.stop
set -u
cd "$(dirname "$0")/.."
GPU=${GPU:-4}
LOG(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
declare -A MF=(
  [0]="./runs/2026-09-15_leipzig-dt01-sponge-uy/out/mean_field.bin"
  [90]="./runs/2026-09-23_leipzig-rot90/out/mean_field.bin"
  [180]="./runs/2026-09-23_leipzig-rot180/out/mean_field.bin"
  [270]="./runs/2026-09-26_leipzig-rot270/out/mean_field.bin" )
COMMON="-D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2
        -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0 -D LDM_MEAN_MODE=1
        -D DD_WALL_DEP=1
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500
        -D LDM_DUMP_EVERY=100"

for a in 0 90 180 270; do
  [ -f runs/chain_src50w.stop ] && { LOG "중단 요청 — 종료"; exit 0; }
  LOG "=== 50소스 벽면침적 ${a}° 시작 ==="
  t0=$(date +%s)
  exp/run.sh "leipzig-src50w-rot$a" --fg -g $GPU -t ldm -i "Leipzig_src50w_rot$a" \
    -m "50 소스 × 풍향 ${a}° 평균장 LDM, 벽면침적 ON. 2026-10-06_leipzig-src50-rot${a}(벽면 off)와 소스·유동장이 동일하고 DD_WALL_DEP=1 만 다르다 — 벽면침적 효과를 같은 소스집합에서 직접 비교하기 위함(09-17 벽면 런은 13소스였다). solv 에 buffer_type 출력을 추가해 0=공중/1=지면/2=벽면·지붕 을 직접 구분한다. DD_wall_dist 기본 2.0 m, 반사 임계 기본(0.5*dx_sp=1.25 m). 소스당 50,000 입자, LDM 출력 1 s." \
    $COMMON -D MEAN_FIELD_PATH="\"${MF[$a]}\"" > /dev/null 2>&1
  rc=$?; d=$(ls -dt runs/*_leipzig-src50w-rot$a 2>/dev/null | head -1)
  LOG "${a}° 종료 rc=$rc  $(( ($(date +%s)-t0)/60 ))분  LDM $(ls $d/out 2>/dev/null | grep -c '^LDM_')개  → $d"
done
LOG "=== 50소스 벽면침적 4방향 완료 ==="
