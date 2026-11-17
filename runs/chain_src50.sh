#!/usr/bin/env bash
# 50 소스 × 4 풍향, 평균장 LDM(Run B). 29소스 재실행이 끝난 뒤 이어서 돈다.
# 소스 선정: 회전중심 r<=410 m(네 방향 풍하여유 246 m), 최소간격 80 m,
#   방출운(반경 15 m) 열린비율 >= 0.90(건물 안으로 새지 않게), λp 8구간 균등, 기존 13개 유지.
# 중단: touch runs/chain_src50.stop
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
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500
        -D LDM_DUMP_EVERY=100"

LOG "29소스 재실행 완료 대기"
while ! grep -q "29소스 4방향 완료" runs/chain_src29b.log 2>/dev/null; do
  tmux has-session -t src29b 2>/dev/null || { grep -q "29소스 4방향 완료" runs/chain_src29b.log 2>/dev/null || { LOG "!! src29b 가 완료 표시 없이 끝남 — 중단"; exit 1; }; }
  sleep 60
done
LOG "확인 — GPU $GPU 해제됨"

for a in 0 90 180 270; do
  [ -f runs/chain_src50.stop ] && { LOG "중단 요청 — 종료"; exit 0; }
  LOG "=== 50소스 ${a}° 시작 ==="
  t0=$(date +%s)
  exp/run.sh "leipzig-src50-rot$a" --fg -g $GPU -t ldm -i "Leipzig_src50_rot$a" \
    -m "50 소스 × 풍향 ${a}° 평균장 LDM. 29소스(n=29)에서 소스를 50개로 더 늘려 검정력을 올린다. 선정: 회전중심 r<=410 m, 최소간격 80 m, 방출운(반경 15 m) 열린비율 >= 0.90, λp 8구간 균등 층화, 기존 13개 유지. 소스당 50,000 입자(총 250만), 질량 31.25 로 총 방출량 고정. 유동장은 각 방향 Run A 의 mean_field.bin 동결. LDM 출력 1 s(100 step). 주의: 간격 80 m 라 λp(75) 창이 35% 겹친다 — 명목 n=50 이지만 유효표본은 그보다 작다." \
    $COMMON -D MEAN_FIELD_PATH="\"${MF[$a]}\"" > /dev/null 2>&1
  rc=$?; d=$(ls -dt runs/*_leipzig-src50-rot$a 2>/dev/null | head -1)
  LOG "${a}° 종료 rc=$rc  $(( ($(date +%s)-t0)/60 ))분  LDM $(ls $d/out 2>/dev/null | grep -c '^LDM_')개  → $d"
done
LOG "=== 50소스 4방향 완료 ==="
