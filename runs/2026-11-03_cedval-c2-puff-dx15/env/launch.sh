#!/bin/bash
# 2026-11-03_cedval-c2-puff-dx15 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=4 exec runs/2026-11-03_cedval-c2-puff-dx15/env/SOPHIA_gpu 1
