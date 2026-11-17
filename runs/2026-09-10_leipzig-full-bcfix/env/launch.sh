#!/bin/bash
# 2026-09-10_leipzig-full-bcfix 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=4 exec runs/2026-09-10_leipzig-full-bcfix/env/SOPHIA_gpu 1
