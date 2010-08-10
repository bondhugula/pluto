SRC=adi

PLCFLAGS = --context=50
TILEFLAGS =  #--l2tile

include ../common.mk

hopt: $(SRC).hopt.c data.decls util.h
	$(CC) $(OPT_FLAGS) -lm $(SRC).hopt.c -o hopt $(CFLAGS)

hopt_test: $(SRC).hopt.c data.decls util.h
	$(CC) $(OPT_FLAGS) -lm $(SRC).hopt.c -o hopt_test -DTEST $(CFLAGS)

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS) --output $(SRC).par2d.c

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

