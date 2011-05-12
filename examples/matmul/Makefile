SRC=matmul

PLCFLAGS = --unroll --ufactor=8
TILEFLAGS = --l2tile
EXTRA_CLEAN = atlas mkl netlib atlaspt goto

include ../common.mk

atlaspt: dgemm.f decls.h
	gfortran -m64 -O dgemm.f -L /usr/local/ATLAS/lib -o atlaspt -llapack -lptcblas -lptf77blas -latlas -lpthread

atlas: dgemm.f decls.h
	gfortran -m64 -O dgemm.f -L /usr/local/ATLAS/lib -o atlas -llapack -lcblas -lf77blas -latlas -lpthread

mkl: matmul.blas.c decls.h
	$(CC) -O3 $(CFLAGS) -I$(MKL)/include matmul.blas.c -L $(MKL)/lib/ia32 -lmkl_intel -lpthread -o mkl

acml: matmul.blas.c decls.h
	$(CC) -O3 -c -I $(ACML)/ifort64_mp/include -openmp matmul.acml.c -DTIME
	ifort -openmp -nofor-main  -L $(ACML)/ifort64_mp/lib matmul.acml.o -lacml_mp -o acml

netlib: dgemm.f decls.h
	ifort -O3 dgemm.f -lblas -o netlib

goto: dgemm.f decls.h
	 ifort -o goto dgemm.f -L ~/tmp/GotoBLAS -lgoto_core2p-r1.23 -lpthread

mkl_test: matmul.blas.c decls.h
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
