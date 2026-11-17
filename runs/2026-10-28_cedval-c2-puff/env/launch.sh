#!/bin/bash
# 2026-10-28_cedval-c2-puff 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=2 exec runs/2026-10-28_cedval-c2-puff/env/SOPHIA_gpu 1
