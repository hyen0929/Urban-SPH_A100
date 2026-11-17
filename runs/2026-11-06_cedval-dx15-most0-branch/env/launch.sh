#!/bin/bash
# 2026-11-06_cedval-dx15-most0-branch 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=4 SOPHIA_RESTART=runs/2026-11-03_cedval-c2-puff-dx15/out/restart_t1000.bin exec runs/2026-11-06_cedval-dx15-most0-branch/env/SOPHIA_gpu 1
