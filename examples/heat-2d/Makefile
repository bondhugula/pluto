
CC=icc

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

lbpar:tile.sizes NToriginal.c
	$(PLC) --tile --parallel --lbtile NToriginal.c -o temp1.c
	bash replace.sh
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) temp2.c -o temp2 $(LDFLAGS)
	bash filename.sh


tiled:tile.sizes NToriginal.c
	$(PLC) --tile --parallel NToriginal.c -o temp1.c
	bash replace1.sh
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) temp2.c -o temp2 $(LDFLAGS)
	bash filename1.sh

orig:tile.sizes NToriginal.c
	cp NToriginal.c original.c
	sed -in '/scop/,/endscop/ s/A\[\([a-z0-9()/+ \t]\{,20\}\)\]/A\[(\1)%2\]/g' original.c 
	$(CC) $(OPT_FLAGS) $(PAR_FLAGS) $(OMP_FLAGS) $(CFLAGS) original.c -o original $(LDFLAGS)

clean:
	rm -f temp* lb-* tp-* original*
