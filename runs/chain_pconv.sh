#!/usr/bin/env bash
# 소스 입자 개수 수렴성 — 중립 기준 평균장(Run B) 위에서 소스당 입자수만 바꾼다.
# 유동장은 2026-09-15 기준런의 mean_field.bin 으로 고정되므로 다섯 런의 유동은 완전히 동일하고,
# 달라지는 것은 LDM 입자 개수뿐이다. 소스별 총 방출량은 질량을 1/N 로 줄여 고정했다.
# GPU 2 (H100 NVL, sm_90) — 다른 작업과 공유. 솔버가 결정론적이라 공유해도 결과는 불변.
# 중단: touch runs/chain_pconv.stop
set -u
cd "$(dirname "$0")/.."
GPU=2; export ARCH=sm_90
MF="./runs/2026-09-15_leipzig-dt01-sponge-uy/out/mean_field.bin"
LOG(){ echo "[$(date '+%m-%d %H:%M')] $*"; }
COMMON="-D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2
        -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0 -D LDM_MEAN_MODE=1
        -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500
        -D LDM_DUMP_EVERY=5000"

for N in 10000 25000 50000 100000 200000; do
  [ -f runs/chain_pconv.stop ] && { LOG "중단 요청 — 종료"; exit 0; }
  LOG "=== N=$N (소스당) 시작 ==="
  t0=$(date +%s)
  exp/run.sh "leipzig-pconv-N$N" --fg -g $GPU -t ldm -i "Leipzig_pconv_N$N" \
    -m "소스 입자수 수렴성 검증: 소스당 $N 개. 중립 기준런(2026-09-15_leipzig-dt01-sponge-uy)의 평균장을 얼려 LDM 만 재계산하므로 다섯 런의 유동장이 완전히 동일하고 입자 개수만 다르다. 소스별 총 방출량은 입자질량을 1/N 로 줄여 고정(N=100000 에서 15.625). 커널 길이 h=3.75 고정. 출력은 실제시간 50초 간격(LDM_DUMP_EVERY=5000 → t=150/200/250/300s). GPU 2 (H100, sm_90) 공유 실행 — 결정론적이라 결과 불변, 벽시계 시간만 늘어남." \
    $COMMON -D MEAN_FIELD_PATH="\"$MF\"" > /dev/null 2>&1
  rc=$?; d=$(ls -dt runs/*_leipzig-pconv-N$N 2>/dev/null | head -1)
  LOG "N=$N 종료 rc=$rc  $(( ($(date +%s)-t0)/60 ))분  LDM $(ls $d/out 2>/dev/null | grep -c '^LDM_')개  → $d"
done
LOG "=== 수렴성 5개 런 완료 ==="
