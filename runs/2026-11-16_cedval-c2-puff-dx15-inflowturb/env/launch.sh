#!/bin/bash
# 2026-11-16_cedval-c2-puff-dx15-inflowturb 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=3 exec runs/2026-11-16_cedval-c2-puff-dx15-inflowturb/env/SOPHIA_gpu 1
