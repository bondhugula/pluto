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
	$(CC) -O3 -DMKL $(CFLAGS) $(MKL_CFLAGS) matmul.blas.c -L $(MKL_LDFLAGS) -o $@

netlib: dgemm.f
	ifort -O3 dgemm.f -lblas -o netlib

goto: dgemm.f
	 ifort -o goto dgemm.f -L ~/tmp/GotoBLAS -lgoto_core2p-r1.23 -lpthread

mkl_test: matmul.blas.c
	icc -DTEST -I/usr/local/mkl/include  -L/usr/local/mkl/lib/em64t $(SRC).blas.c -o mkl_test -lmkl_em64t -lmkl -lguide -lpthread

ptest: tiled_test par_test mkl_test  orig_test
	./orig_test 2> out_orig
	./tiled_test 2> out_tiled
	diff -q out_tiled out_orig
	export OMP_NUM_THREADS=4
	./par_test 2> out_par4
	diff -q out_orig out_par4
	./mkl_test 2> out_mkl
	diff -q out_tiled out_mkl
