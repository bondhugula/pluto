SRC= gemver

PLCFLAGS = --context=1 
TILEFLAGS = --tile --l2tile

include ../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS)

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

par2d_test: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d_test -DTEST

mkl: $(SRC).blas.c decls.h
	$(CC) -lm $(OPT_FLAGS) $(CFLAGS) -I/usr/local/mkl/include $(SRC).blas.c -L /usr/local/mkl/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread -o mkl

mkl_test: $(SRC).blas.c decls.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) -DTEST -I/usr/local/mkl/include  -L/usr/local/mkl/lib/em64t $(SRC).blas.c -o mkl_test -lmkl_em64t -lmkl -lguide -lpthread

ptest: orig_test mkl_test tiled_test par_test
	./orig_test 2> out_orig
	./tiled_test 2> out_tiled
	diff -q out_tiled out_orig
	export OMP_NUM_THREADS=4 
	./par_test 2> out_par4
	diff -q out_par4 out_orig
	./mkl_test 2> out_mkl
	diff -q out_mkl out_orig
