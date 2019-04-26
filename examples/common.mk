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
POLYBENCHINCDIR=$(BASEDIR)polybench/utilities
POLYBENCHSRC=$(BASEDIR)polybench/utilities/polybench.c
PLC=$(BASEDIR)../polycc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

ifeq ($(CC), icc)
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -qopenmp
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

ifdef POLYBENCH
	CFLAGS += -DPOLYBENCH_USE_SCALAR_LB -DPOLYBENCH_TIME -I $(POLYBENCHINCDIR) $(POLYBENCHSRC)
	DISTOPT_FLAGS += --variables_not_global
endif

all: orig tiled par

$(SRC).opt.c:  $(SRC).c
	$(PLC) $(SRC).c --notile --noparallel $(PLCFLAGS)  -o $@

$(SRC).tiled.c:  $(SRC).c
	$(PLC) $(SRC).c --noparallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --full-diamond-tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).mlbpar.c:  $(SRC).c
	$(PLC) $(SRC).c $(TILEFLAGS) $(PLCFLAGS)  -o $@

# Version that doesn't use diamond tiling
$(SRC).pipepar.c:  $(SRC).c
	$(PLC) $(SRC).c --nodiamond-tile $(TILEFLAGS) $(PLCFLAGS) -o $@

orig: $(SRC).c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c -o $@ $(LDFLAGS)

orig_par: $(SRC).c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(PAR_FLAGS) $(SRC).c -o $@ $(LDFLAGS)

opt: $(SRC).opt.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).opt.c -o $@ $(LDFLAGS)

tiled: $(SRC).tiled.c 
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).tiled.c -o $@ $(LDFLAGS)

mlbpar: $(SRC).mlbpar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).mlbpar.c -o $@  $(LDFLAGS)

# Version that doesn't use diamond tiling
pipepar: $(SRC).pipepar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).pipepar.c -o $@  $(LDFLAGS)

perf: orig tiled par orig_par
	rm -f .test
	./orig
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par
	./tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 

# Compare performance with and without diamond tiling.
pipeperf: par pipepar
	rm -f .test
	OMP_NUM_THREADS=$(NTHREADS) ./par
	OMP_NUM_THREADS=$(NTHREADS) ./pipepar 


test: orig tiled par
	touch .test
	./orig 2> out_orig
	./tiled 2> out_tiled
	diff -q out_orig out_tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	rm -f .test
	diff -q out_orig out_par4
	@echo Success!

lbtest: par pipepar
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	OMP_NUM_THREADS=$(NTHREADS) ./pipepar 2> out_pipepar4
	rm -f .test
	diff -q out_par4 out_pipepar4
	diff -q out_par4 out_fulldiamondtile4
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
	rm -f out_* *.pipepar.c *.tiled.c *.opt.c *.par.c orig opt tiled par sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog *.pluto.cloog

exec-clean:
	rm -f out_* opt orig tiled sched sched hopt hopt par pipepar orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
