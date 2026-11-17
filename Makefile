# SOPHIA_gpu build
#
#   make                     기존과 동일 (Parameters.cuh 기본값으로 ./SOPHIA_gpu 생성)
#   make BIN=path EXP="-DA=1 -DB=2"
#                            실험용 빌드. EXP 의 -D 플래그가 Parameters.cuh 의
#                            #ifndef 가드를 덮어쓴다 (exp/run.sh 가 사용).
#   make ARCH=sm_86          GPU 아키텍처 변경
#   make clean

NVCC      ?= nvcc
ARCH      ?= sm_80
BIN       ?= SOPHIA_gpu
EXP       ?=
NVCCFLAGS ?= -w -use_fast_math -arch=$(ARCH) -O3 -expt-relaxed-constexpr -I./cub-1.8.0/
LDFLAGS   ?= -lpthread

.PHONY: make all clean
make: all

all:
	$(NVCC) $(NVCCFLAGS) $(EXP) SOPHIA_gpu.cu -o $(BIN) $(LDFLAGS)

clean:
	rm -f SOPHIA_gpu
