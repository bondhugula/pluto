#
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#
BASEDIR=$(dir $(lastword $(MAKEFILE_LIST)))

#CC=icc
CC=gcc

NPROCS=4
NTHREADS=4
PLC=$(BASEDIR)../polycc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

ifeq ($(CC), icc)
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -openmp
else
	# for gcc
	OPT_FLAGS     := -O3 -march=native -mtune=native -ftree-vectorize
	PAR_FLAGS     := -ftree-parallelize-loops=$(NTHREADS)
	OMP_FLAGS     := -fopenmp
endif

CFLAGS += -DTIME
LDFLAGS += -lm
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

$(SRC).lbpar.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel --partlbtile $(TILEFLAGS) $(PLCFLAGS) -o $@


orig: $(SRC).c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c -o $@ $(LDFLAGS)

orig_par: $(SRC).c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(PAR_FLAGS) $(SRC).c -o $@ $(LDFLAGS)

opt: $(SRC).opt.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).opt.c -o $@ $(LDFLAGS)

tiled: $(SRC).tiled.c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).tiled.c -o $@ $(LDFLAGS)

lbpar: $(SRC).lbpar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).lbpar.c -o $@  $(LDFLAGS)

par: $(SRC).par.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o $@  $(LDFLAGS)

perf: orig tiled par orig_par
	rm -f .test
	./orig
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par
	./tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 

lbperf: par lbpar
	rm -f .test
	OMP_NUM_THREADS=$(NTHREADS) ./par
	OMP_NUM_THREADS=$(NTHREADS) ./lbpar 


test: orig tiled par
	touch .test
	./orig 2> out_orig
	./tiled 2> out_tiled
	diff -q out_orig out_tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	rm -f .test
	diff -q out_orig out_par4
	@echo Success!

lbtest: par lbpar
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	OMP_NUM_THREADS=$(NTHREADS) ./lbpar 2> out_lbpar4
	rm -f .test
	diff -q out_par4 out_lbpar4
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
	rm -f out_* *.lbpar.c *.tiled.c *.opt.c *.par.c orig opt tiled par sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog *.pluto.cloog

exec-clean:
	rm -f out_* opt orig tiled lbtile lbpar  sched sched hopt hopt par orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
