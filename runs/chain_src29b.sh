#!/usr/bin/env bash
# 29 소스 × 4 풍향, 평균장 LDM(Run B). 각 방향의 Run A 에서 얻은 mean_field.bin 을 얼려 쓴다.
# 소스는 회전 중심 (657.5,0) 기준 r<=411 m 안에 있어 네 방향 모두 풍하여유 246 m 이상을 확보한다.
# 소스당 50,000 개 — 입자수 수렴성 검증(runs/lambdap/pconv)에서 안전선으로 확인된 값.
# 중단: touch runs/chain_src29b.stop
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

for a in 0 90 180 270; do
  [ -f runs/chain_src29b.stop ] && { LOG "중단 요청 — 종료"; exit 0; }
  [ -s "${MF[$a]#./}" ] || { LOG "!! ${a}° mean_field 없음 — 건너뜀"; continue; }
  LOG "=== ${a}° 시작 ==="
  t0=$(date +%s)
  exp/run.sh "leipzig-src29b-rot$a" --fg -g $GPU -t ldm -i "Leipzig_src29_rot$a" \
    -m "29 소스 × 풍향 ${a}° 평균장 LDM — 재실행본. 09-30 판은 LDM 출력 주기를 5000 step(50 s)으로 잘못 잡아 덤프가 7개뿐이었다(입자수 수렴성 실험 플래그를 그대로 복사한 실수). LDM 은 1 s(100 step) 간격이어야 플룸 나이에 따른 상관 변화를 볼 수 있다. 기존 13 소스의 n 부족(ρs 신뢰구간 폭 0.89)을 보완하려고 소스를 29개로 늘렸다. 선정 기준: 회전 중심에서 r<=410 m(네 방향 모두 풍하여유 246 m 확보), 최소간격 120 m(λp(75) 창 겹침 10%), λp 8구간 균등 층화, 기존 13개 유지. 소스당 50,000 입자 — 수렴성 검증에서 지표가 이미 수렴한 값(runs/lambdap/pconv). 총 방출량은 질량 31.25 로 고정. 유동장은 각 방향 Run A 의 mean_field.bin 을 동결해 쓰므로 방향 간 유동 차이만 남는다. 목적: 형태-확산 상관이 풍향을 넘어 살아남는지 n=29 로 검정." \
    $COMMON -D MEAN_FIELD_PATH="\"${MF[$a]}\"" > /dev/null 2>&1
  rc=$?; d=$(ls -dt runs/*_leipzig-src29b-rot$a 2>/dev/null | head -1)
  LOG "${a}° 종료 rc=$rc  $(( ($(date +%s)-t0)/60 ))분  LDM $(ls $d/out 2>/dev/null | grep -c '^LDM_')개  → $d"
done
LOG "=== 29소스 4방향 완료 ==="
