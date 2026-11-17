#!/bin/bash
# 2026-10-06_leipzig-src29b-rot270 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=4 exec runs/2026-10-06_leipzig-src29b-rot270/env/SOPHIA_gpu 1
