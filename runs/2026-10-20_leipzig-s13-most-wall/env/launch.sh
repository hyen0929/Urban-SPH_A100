#!/bin/bash
# 2026-10-20_leipzig-s13-most-wall 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=1 exec runs/2026-10-20_leipzig-s13-most-wall/env/SOPHIA_gpu 1
