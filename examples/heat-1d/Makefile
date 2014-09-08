SRC=heat-1d

include ../common.mk

NPROCS=4

ifeq ($(CC), icc)
	OPT_FLAGS=-O3 -fp-model precise
	PAR_FLAGS := -parallel
	OMP_FLAGS := -openmp
else
	# for gcc
	OPT_FLAGS := -O3 -ftree-vectorize -msse3 
	PAR_FLAGS := -ftree-parallelize-loops=4
	OMP_FLAGS := -fopenmp
endif

CFLAGS = -DTIME -DVERIFY
LDFLAGS = -lm

PLC=../../polycc

PLCFLAGS+=--pet
