#!/bin/bash
# 2026-09-14_leipzig-dt01-full 재현 실행 (리포 루트에서)
cd "$(dirname "$0")/../../.."
# stdout 줄 버퍼링(stdbuf -oL): 블록 버퍼링이면 run.log 가 ~40 분 늦게 찍혀 감시가 늦는다
CUDA_VISIBLE_DEVICES=4 exec stdbuf -oL runs/2026-09-14_leipzig-dt01-full/env/SOPHIA_gpu 1
