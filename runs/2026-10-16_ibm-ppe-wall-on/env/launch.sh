#!/bin/bash
# 2026-10-16_ibm-ppe-wall-on 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
CUDA_VISIBLE_DEVICES=1 SOPHIA_RESTART=runs/2026-09-15_leipzig-dt01-sponge-uy/out/restart.bin exec runs/2026-10-16_ibm-ppe-wall-on/env/SOPHIA_gpu 1
