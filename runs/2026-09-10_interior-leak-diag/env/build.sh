#!/bin/bash
# 2026-09-10_interior-leak-diag 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=0 -DSPONGE_GAMMA_MAX=0 -DINIT_CROSSFLOW_UY=0.0 -DAVG_START_TIME=6.0 -DAVG_END_TIME=10.0 -DLDM_DURATION_AFTER_AVG=4.0 -DFULL_DUMP_EVERY=200 -DVORT_DUMP_EVERY=200 -DINPUT_CASE_NAME=\"Leipzig_surface_s10\" -DOUTPUT_CASE_NAME=\"runs/2026-09-10_interior-leak-diag/out\"  SOPHIA_gpu.cu -o runs/2026-09-10_interior-leak-diag/env/SOPHIA_gpu -lpthread
