SRC=dct

PLCFLAGS = --context=32
TILEFLAGS = --l2tile

include ../common.mk

pouchet: $(SRC).pouchet.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).pouchet.c -o pouchet  $(CFLAGS)

