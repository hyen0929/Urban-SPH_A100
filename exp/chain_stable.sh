#!/usr/bin/env bash
# exp/chain_stable.sh — 안정 성층 Run A 완주를 기다렸다가 평균장 LDM(Run B)을 잇는다.
#
#   nohup bash exp/chain_stable.sh > runs/chain_stable.log 2>&1 &
#
# 중단: runs/chain_stable.stop 파일을 만들면 Run B 를 띄우지 않고 끝난다.
set -u
cd "$(dirname "$0")/.."
GPU=4
AD=runs/2026-09-21_leipzig-most-stable
LOG() { echo "[$(date '+%m-%d %H:%M')] $*"; }

PID=$(cat "$AD/env/pid" 2>/dev/null) || { LOG "!! $AD/env/pid 없음"; exit 1; }
LOG "대기: $(basename "$AD")  PID $PID"
while kill -0 "$PID" 2>/dev/null; do sleep 60; done

if grep -qiE "nan|segmentation|cudaError|illegal|core dumped|out of memory" "$AD/run.log" 2>/dev/null; then
  LOG "!! Run A 가 에러로 끝났다 — Run B 취소"
  grep -inE "nan|segmentation|cudaError|illegal|core dumped|out of memory" "$AD/run.log" | tail -5
  exit 1
fi
LOG "완료: $(basename "$AD")"

[ -f runs/chain_stable.stop ] && { LOG "중단 요청 확인 — Run B 생략"; exit 0; }

if [ ! -s "$AD/out/mean_field.bin" ]; then
  LOG "!! mean_field.bin 이 없다 — Run B 건너뜀 (AVG_END_TIME 이 end-time 보다 앞서야 덤프된다)"
  exit 1
fi

LOG "=== 안정 성층 평균장 LDM(Run B) 시작 ==="
exp/run.sh leipzig-most-stable-meanflow \
  -m "안정 성층(T_con=285K) Run A 의 Run B — frozen mean-flow LDM. 같은 런이 t=299s 에 덤프한 평균장으로 LDM 만 재계산한다. 방출·난수 동일, 차이는 분해된 비정상성의 유무뿐. 건성침적이 쓰는 us_update/OB_length_global 은 mean_field.bin 헤더에 함께 저장돼 ISPH_Calc.cuh:235 에서 복원되므로 MOST_boundary=0 으로도 안정 성층의 u*·Obukhov 길이를 물려받는다(단 t=299s 순간값이지 시간평균은 아님). 목적: 중립에서 얻은 전이 대 평균장 차이(총 노출 x1.03, 셀 순위상관 0.63)가 안정 성층에서도 같은지." \
  -i Leipzig_stable -t ldm -g $GPU \
  -D CASE_ID=1002 -D DEACTIVATE_IBM=0 -D MOST_boundary=0 -D SPONGE_GAMMA_MAX=0.2 \
  -D SPONGE_DAMP_UY=1 -D INIT_CROSSFLOW_UY=0.0 -D LDM_MEAN_MODE=1 \
  -D "MEAN_FIELD_PATH=\"./$AD/out/mean_field.bin\"" \
  -D SLICE_DUMP_EVERY=50 -D FULL_DUMP_EVERY=5000 -D VORT_DUMP_EVERY=500 >/dev/null 2>&1

sleep 10
BD=$(ls -dt runs/*_leipzig-most-stable-meanflow 2>/dev/null | head -1)
LOG "Run B 기동: ${BD:-실패}"
