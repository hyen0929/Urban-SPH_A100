#!/bin/bash
# 2026-09-14_ckpt-restart-verify 재현 빌드
set -e
cd "$(dirname "$0")/../../.."
nvcc -w -use_fast_math -arch=sm_80 -O3 -expt-relaxed-constexpr -I./cub-1.8.0/ -DCASE_ID=1002 -DDEACTIVATE_IBM=0 -DMOST_boundary=0 -DSPONGE_GAMMA_MAX=0 -DAVG_START_TIME=0.2 -DAVG_END_TIME=0.6 -DINIT_CROSSFLOW_UY=0.0 -DCKPT_EVERY=40 -DSLICE_DUMP_EVERY=25 -DFULL_DUMP_EVERY=75 -DINPUT_CASE_NAME=\"CkptTest\" -DOUTPUT_CASE_NAME=\"runs/2026-09-14_ckpt-restart-verify/out\"  SOPHIA_gpu.cu -o runs/2026-09-14_ckpt-restart-verify/env/SOPHIA_gpu -lpthread
