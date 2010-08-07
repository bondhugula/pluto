SRC= doitgen

PLCFLAGS = --context=10 --prevector #--unroll --ufactor=8 #--maxfuse
TILEFLAGS = #--l2tile

include ../common.mk

mkl: doitgen.mkl.c decls.h
	$(CC) $(CFLAGS) -I$(MKL)/include doitgen.mkl.c -L $(MKL)/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread -o mkl

acml: doitgen.acml.c decls.h
	$(CC) -c $(CFLAGS) -I $(ACML)/ifort64_mp_int64/include doitgen.acml.c
	ifort -openmp -nofor-main  -L $(ACML)/ifort64_mp_int64/lib doitgen.acml.o -lacml_mp -o acml

mkl2: doitgen.dgemm.c decls.h
	$(CC) $(CFLAGS) -I $(MKL)/include doitgen.dgemv.c -L $(MKL)/lib/em64t -lmkl_em64t -lmkl -lguide -lpthread -o mkl2

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS)

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d $(CFLAGS)

par2d_test: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d_test -DTEST

ptest: tiled_test par2d_test par_test
	./tiled_test 2> out_tiled
	export OMP_NUM_THREADS=4 
	./par_test 2> out_par4
	diff -q out_tiled out_par4
	export OMP_NUM_THREADS=4; ./par2d_test 2> out_par2d
	diff -q out_tiled out_par2d
