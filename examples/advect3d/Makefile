SRC=advect3d

PLCFLAGS = --context=100
TILEFLAGS = #--l2tile

include ../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS) 

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d


hopt: $(SRC).hopt.c data.decls util.h
	$(CC) $(OPT_FLAGS) -lm $(SRC).hopt.c -o hopt

hopt_test: $(SRC).hopt.c data.decls util.h
	$(CC) $(OPT_FLAGS) -lm $(SRC).hopt.c -o hopt_test -DTEST
