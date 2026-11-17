#!/bin/bash
# 2026-10-16_ibm-ppe-wall-off 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=0 -DSPONGE_GAMMA_MAX=0.2 -DSPONGE_DAMP_UY=1 -DINIT_CROSSFLOW_UY=0.0 -DSLICE_DUMP_EVERY=50 -DFULL_DUMP_EVERY=500 -DVORT_DUMP_EVERY=500 -DINPUT_CASE_NAME=\"Leipzig_wallbc\" -DOUTPUT_CASE_NAME=\"runs/2026-10-16_ibm-ppe-wall-off/out\"  SOPHIA_gpu.cu -o runs/2026-10-16_ibm-ppe-wall-off/env/SOPHIA_gpu -lpthread
