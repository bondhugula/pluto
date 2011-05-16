#
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#

#CC=icc
CC=gcc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

ifeq ($(CC), icc)
	OPT_FLAGS=-O3 -I/usr/include -fp-model precise
	PAR_FLAGS := -parallel
	OMP_FLAGS := -openmp
else
	# for gcc
	OPT_FLAGS := -O3 -ftree-vectorize -msse3 
	PAR_FLAGS := -ftree-parallelize-loops=4
	OMP_FLAGS := -fopenmp
endif

CFLAGS = -DTIME
LDFLAGS = -lm
PLCFLAGS +=
TILEFLAGS += 

#PERFCTR=perfctr

ifdef PERFCTR
	CFLAGS += -DPERFCTR -L/usr/local/lib64 -lpapi
endif

PLC=../../polycc

all: orig tiled par

orig: $(SRC).c decls.h  util.h
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c -o orig $(LDFLAGS)

orig_par: decls.h  util.h $(SRC).c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(PAR_FLAGS) $(SRC).c -o orig_par $(LDFLAGS)

$(SRC).opt.c: 
	$(PLC) $(SRC).c $(PLCFLAGS) 

$(SRC).tiled.c: 
	$(PLC) $(SRC).c --tile $(TILEFLAGS) $(PLCFLAGS) 

$(SRC).par.c: 
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS) 

opt: $(SRC).opt.c decls.h  util.h
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).opt.c -o opt $(LDFLAGS)

tiled: $(SRC).tiled.c decls.h  util.h
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).tiled.c -o tiled $(LDFLAGS)

par: $(SRC).par.c decls.h  util.h
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o par  $(LDFLAGS)

perf: orig tiled par orig_par
	rm -f .test
	./orig
	OMP_NUM_THREADS=4 ./orig_par
	./tiled
	OMP_NUM_THREADS=4 ./par 


test: orig tiled par
	touch .test
	./orig > out_orig
	./tiled > out_tiled
	diff -q out_orig out_tiled
	export OMP_NUM_THREADS=4; ./par > out_par4
	diff -q out_orig out_par4
	@echo Success!
	rm -f .test

opt-test: orig opt
	touch .test
	./orig > out_orig
	./opt > out_opt
	diff -q out_orig out_opt

clean:
	rm -f out_* *.tiled.c *.opt.c *.par.c orig opt tiled par sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog

exec-clean:
	rm -f out_* opt orig tiled  sched sched hopt hopt par orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
