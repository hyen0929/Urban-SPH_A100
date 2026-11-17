#!/bin/bash
# 2026-09-23_leipzig-rot180 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=2 exec runs/2026-09-23_leipzig-rot180/env/SOPHIA_gpu 1
