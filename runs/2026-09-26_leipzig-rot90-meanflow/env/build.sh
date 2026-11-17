#!/bin/bash
# 2026-09-26_leipzig-rot90-meanflow 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=0 -DSPONGE_GAMMA_MAX=0.2 -DSPONGE_DAMP_UY=1 -DINIT_CROSSFLOW_UY=0.0 -DSLICE_DUMP_EVERY=50 -DFULL_DUMP_EVERY=5000 -DVORT_DUMP_EVERY=500 -DLDM_DUMP_EVERY=100 -DLDM_MEAN_MODE=1 -DMEAN_FIELD_PATH=\"./runs/2026-09-23_leipzig-rot90/out/mean_field.bin\" -DINPUT_CASE_NAME=\"Leipzig_rot90\" -DOUTPUT_CASE_NAME=\"runs/2026-09-26_leipzig-rot90-meanflow/out\"  SOPHIA_gpu.cu -o runs/2026-09-26_leipzig-rot90-meanflow/env/SOPHIA_gpu -lpthread
