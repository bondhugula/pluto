SRC=fdtd-2d

PLCFLAGS += #--context=32  #--unroll
TILEFLAGS = #--l2tile

include ../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS) 

$(SRC).sched.c: 
	$(PLC) $(SRC).c --sched --tile --parallel $(TILEFLAGS) $(PLCFLAGS) 

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

sched: $(SRC).sched.c decls.h util.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) -openmp -lm $(SRC).sched.c -o sched

par2d_test: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d_test -DTEST

sched_test: $(SRC).sched.c decls.h util.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) -openmp -lm $(SRC).sched.c -o sched_test -DTEST

ptest: tiled_test sched_test par2d_test par_test
	./tiled_test 2> out_tiled
	export OMP_NUM_THREADS=4 
	./par_test 2> out_par4
	diff -q out_tiled out_par4
	./par2d_test 2> out_par2d
	diff -q out_tiled out_par2d
	./sched_test 2> out_sched
	diff -q out_tiled out_sched
