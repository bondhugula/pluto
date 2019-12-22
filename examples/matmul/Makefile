SRC=matmul

PLCFLAGS = --second-level-tile
TILEFLAGS =
EXTRA_CLEAN = atlas mkl netlib atlaspt goto

include ../common.mk

atlaspt: dgemm.f
	gfortran -m64 -O dgemm.f -L /usr/local/ATLAS/lib -o atlaspt -llapack -lptcblas -lptf77blas -latlas -lpthread

atlas: dgemm.f
	gfortran -m64 -O dgemm.f -L /usr/local/ATLAS/lib -o atlas -llapack -lcblas -lf77blas -latlas -lpthread

openblas: matmul.blas.c
	$(CC) -O3 -DOPENBLAS $(CFLAGS) $(OPENBLAS_CFLAGS) matmul.blas.c $(OPENBLAS_LDFLAGS) -o $@

blis: matmul.blas.c
	$(CC) -O3 -DBLIS $(CFLAGS) $(BLIS_CFLAGS) matmul.blas.c -lblis  -o $@

mkl: matmul.blas.c
	$(CC) -O3 -DMKL $(CFLAGS) $(MKL_CFLAGS) matmul.blas.c $(MKL_LDFLAGS) -o $@

netlib: dgemm.f
	ifort -O3 dgemm.f -lblas -o netlib

goto: dgemm.f
	 ifort -o goto dgemm.f -L ~/tmp/GotoBLAS -lgoto_core2p-r1.23 -lpthread

blas_test: tiled mkl openblas blis
	touch .test
	./tiled 2> out_tiled
	OMP_NUM_THREADS=1 ./mkl 2> out_mkl
	diff -q out_tiled out_mkl
	BLIS_NUM_THREADS=1 ./blis 2> out_blis
	diff -q out_tiled out_blis
	OPENBLAS_NUM_THREADS=1 ./openblas 2> out_openblas
	diff -q out_tiled out_mkl out_blis out_openblas
	rm -f .test
