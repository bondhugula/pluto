SRC=fpc_d2q9

PLCFLAGS_OLD_PLUTO = ${PLCFLAGS} --pet --islsolve
PLCFLAGS+=--pet --iss --islsolve

CFLAGS += -std=c99

include ../../../common.mk

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

$(SRC).par-pluto.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS_OLD_PLUTO)  -o $@

par-pluto: $(SRC).par-pluto.c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par-pluto.c -o $@ $(LDFLAGS)
