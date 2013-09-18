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

all: orig tiled lbpar

lbpar:tile.sizes $(SRC).c
	$(PLC) --tile --parallel --lbtile $(SRC).c -o temp1.c
	bash replace.sh
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) temp2.c -o temp2 $(LDFLAGS)
	bash filename.sh


tiled:tile.sizes $(SRC).c
	$(PLC) --tile --parallel $(SRC).c -o temp1.c
	bash replace1.sh
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) temp2.c -o temp2 $(LDFLAGS)
	bash filename1.sh

orig:tile.sizes $(SRC).c
	cp $(SRC).c original.c
	sed -in '/scop/,/endscop/ s/A\[\([a-z0-9()/+ \t]\{,20\}\)\]/A\[(\1)%2\]/g' original.c 
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) original.c -o original $(LDFLAGS)

lbtest: par lbpar
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	OMP_NUM_THREADS=$(NTHREADS) ./lbpar 2> out_lbpar4
	diff -q out_par4 out_lbpar4
	rm -f .test
	@echo Success!


clean:
	rm -f temp* lb-* tp-* original*
