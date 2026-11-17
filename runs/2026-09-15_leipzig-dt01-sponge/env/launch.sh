#!/bin/bash
# 2026-09-15_leipzig-dt01-sponge 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
# 이 런은 분기 재개: 2026-09-14_leipzig-dt01-full 의 count=11000 체크포인트에서 시작한다.
# 중간에 죽으면 이 launch.sh 가 아니라 env/resume.sh 로 (자기 out/restart.bin 에서) 재개할 것.
CUDA_VISIBLE_DEVICES=4 SOPHIA_RESTART=runs/2026-09-14_leipzig-dt01-full/out/restart_t110.bin exec stdbuf -oL runs/2026-09-15_leipzig-dt01-sponge/env/SOPHIA_gpu 1
