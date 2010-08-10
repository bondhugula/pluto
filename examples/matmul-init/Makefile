SRC=matmul-init

PLCFLAGS = --context=1000
TILEFLAGS = --l2tile #--unroll

include ../common.mk

pouchet: $(SRC).pouchet.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).pouchet.c -o pouchet  $(CFLAGS)

pouchet_test: $(SRC).pouchet.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).pouchet.c -o pouchet_test  $(CFLAGS) -DTEST


ptest: opt_test tiled_test orig_test par_test pouchet_test
	./orig_test 2> out_orig
	./tiled_test 2> out_tiled
	diff -q out_orig out_tiled
	export OMP_NUM_THREADS=4; ./par_test 2> out_par4
	diff -q out_orig out_par4
	export OMP_NUM_THREADS=4; ./pouchet_test 2> out_pouchet4
	diff -q out_orig out_pouchet4
	@echo Success!

clean:
