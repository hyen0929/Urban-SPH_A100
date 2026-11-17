#!/bin/bash
# 2026-09-10_interior-leak-diag 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=0 exec runs/2026-09-10_interior-leak-diag/env/SOPHIA_gpu 1
