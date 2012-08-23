#
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#

CC=icc
#CC=gcc

NPROCS=4

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

ifeq ($(CC), icc)
	OPT_FLAGS=-O3 -fp-model precise
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

$(SRC).opt.c:  $(SRC).c
	$(PLC) $(SRC).c $(PLCFLAGS)  -o $@

$(SRC).tiled.c:  $(SRC).c
	$(PLC) $(SRC).c --tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

orig: $(SRC).c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c -o $@ $(LDFLAGS)

orig_par: $(SRC).c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(PAR_FLAGS) $(SRC).c -o $@ $(LDFLAGS)

opt: $(SRC).opt.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).opt.c -o $@ $(LDFLAGS)

tiled: $(SRC).tiled.c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).tiled.c -o $@ $(LDFLAGS)

par: $(SRC).par.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o $@  $(LDFLAGS)

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
	rm -f .test
	diff -q out_orig out_par4
	@echo Success!

opt-test: orig opt
	touch .test
	./orig > out_orig
	./opt > out_opt
	rm -f .test
	diff -q out_orig out_opt
	@echo Success!
	rm -f .test

clean:
	rm -f out_* *.tiled.c *.opt.c *.par.c orig opt tiled par sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog

exec-clean:
	rm -f out_* opt orig tiled  sched sched hopt hopt par orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
