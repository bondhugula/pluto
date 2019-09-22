SRC=dct

PLCFLAGS = --codegen-context=32
TILEFLAGS = --second-level-tile

include ../common.mk

pouchet: $(SRC).pouchet.c decls.h papi_defs.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).pouchet.c -o pouchet  $(CFLAGS)

