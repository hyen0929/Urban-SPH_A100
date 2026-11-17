#!/bin/bash
# 2026-09-23_stable-meanflow-localust 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=0 -DSPONGE_GAMMA_MAX=0.2 -DSPONGE_DAMP_UY=1 -DINIT_CROSSFLOW_UY=0.0 -DLDM_MEAN_MODE=1 -DSLICE_DUMP_EVERY=50 -DFULL_DUMP_EVERY=5000 -DVORT_DUMP_EVERY=500 -DL_MOST_z_lo=0.0 -DL_MOST_z_hi=2.5 -DMEAN_FIELD_PATH=\"./runs/2026-09-21_leipzig-most-stable/out/mean_field.bin\" -DINPUT_CASE_NAME=\"Leipzig_stable\" -DOUTPUT_CASE_NAME=\"runs/2026-09-23_stable-meanflow-localust/out\"  SOPHIA_gpu.cu -o runs/2026-09-23_stable-meanflow-localust/env/SOPHIA_gpu -lpthread
