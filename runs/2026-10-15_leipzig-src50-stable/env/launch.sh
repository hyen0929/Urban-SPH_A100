#!/bin/bash
# 2026-10-15_leipzig-src50-stable 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=3 exec runs/2026-10-15_leipzig-src50-stable/env/SOPHIA_gpu 1
