SRC=jacobi-1d-imper

PLCFLAGS = #--unroll
TILEFLAGS = #--l2tile

include ../common.mk


$(SRC).sched.c: 
	$(PLC) $(SRC).c --sched --tile --parallel $(TILEFLAGS) $(PLCFLAGS) -o $(SRC).sched.c

sched: $(SRC).sched.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).sched.c -o sched

sched_test: $(SRC).sched.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).sched.c -o sched_test -DTEST

limlam: $(SRC).limlam.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).limlam.c -o limlam

