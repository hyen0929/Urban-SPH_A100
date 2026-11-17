#!/bin/bash
# 2026-09-28_leipzig-rot270-meanflow 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=4 exec runs/2026-09-28_leipzig-rot270-meanflow/env/SOPHIA_gpu 1
