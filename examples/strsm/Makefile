SRC= strsm

PLCFLAGS = #--unroll
TILEFLAGS = #--l2tile

include ../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS)

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

par2d_test: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d_test -DTEST

ptest: tiled_test par2d_test par_test
	./tiled_test 2> out_tiled
	export OMP_NUM_THREADS=4 
	./par_test 2> out_par4
	diff -q out_tiled out_par4
	export OMP_NUM_THREADS=4; ./par2d_test 2> out_par2d
	diff -q out_tiled out_par2d
