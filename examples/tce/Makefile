SRC=tce-4index-transform

PLCFLAGS=--context=32
TILEFLAGS = 

include ../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS) --output $(SRC).kernel.par2d.c

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

