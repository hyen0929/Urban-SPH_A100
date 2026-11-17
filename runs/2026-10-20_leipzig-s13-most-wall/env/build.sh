#!/bin/bash
# 2026-10-20_leipzig-s13-most-wall 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=1 -DL_MOST_z_lo=0.0 -DL_MOST_z_hi=2.5 -DABL_h=1.0e9 -DSPONGE_GAMMA_MAX=0.2 -DSPONGE_DAMP_UY=1 -DINIT_CROSSFLOW_UY=0.0 -DDD_WALL_DEP=1 -DAVG_START_TIME=150.0 -DAVG_END_TIME=299.0 -DCKPT_EVERY=1000 -DLDM_DUMP_EVERY=100 -DSLICE_DUMP_EVERY=50 -DFULL_DUMP_EVERY=5000 -DVORT_DUMP_EVERY=500 -DINPUT_CASE_NAME=\"Leipzig_s10_wall\" -DOUTPUT_CASE_NAME=\"runs/2026-10-20_leipzig-s13-most-wall/out\"  SOPHIA_gpu.cu -o runs/2026-10-20_leipzig-s13-most-wall/env/SOPHIA_gpu -lpthread
