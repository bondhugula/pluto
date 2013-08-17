#
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#
BASEDIR=$(dir $(lastword $(MAKEFILE_LIST)))

CC=icc
#CC=gcc

NPROCS=4
NTHREADS=4
NPROCS_LIST = 1 2 4 8 16 32
TILE_SIZES = 16 32 64 128 256
POLYRTINCDIR=$(BASEDIR)../polyrt
POLYRTLIBDIR=$(BASEDIR)../polyrt/.libs
POLYBENCHINCDIR=$(BASEDIR)polybench/utilities
POLYBENCHSRC=$(BASEDIR)polybench/utilities/polybench.c
HOSTS_FILE=$(BASEDIR)hosts

PLC=$(BASEDIR)../polycc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

ifeq ($(CC), icc)
	OPT_FLAGS     := -O3 -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -openmp
	MPICC         := OMPI_CC=icc mpicc # OPENMPI
	#uncomment the following line and comment the previous line when running on fist (MVAPICH/MPICH)
	#MPICC         := mpicc -cc=icc
else
	# for gcc
	OPT_FLAGS     := -O3 -ftree-vectorize -msse3 
	PAR_FLAGS     := -ftree-parallelize-loops=4
	OMP_FLAGS     := -fopenmp
	MPICC         := OMPI_CC=gcc mpicc # OPENMPI
	#uncomment the following line and comment the previous line when running on fist (MVAPICH/MPICH)
	#MPICC         := mpicc -cc=gcc
endif

CFLAGS += -DTIME
LDFLAGS += -lm
PLCFLAGS += --isldep --lastwriter
DISTOPT_FLAGS += --cloogsh
TILEFLAGS += 

#PERFCTR=perfctr

ifdef PERFCTR
	CFLAGS += -DPERFCTR -L/usr/local/lib64 -lpapi
endif

PLC=../../polycc

all: orig tiled par distopt

$(SRC).opt.c:  $(SRC).c
	$(PLC) $(SRC).c $(PLCFLAGS)  -o $@

$(SRC).tiled.c:  $(SRC).c
	$(PLC) $(SRC).c --tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).distopt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --commreport --mpiomp --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@


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

distopt: $(SRC).distopt.c sigma.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI $(CFLAGS) $(SRC).distopt.c sigma.c pi.c \
		-o $@ -L $(POLYRTLIBDIR) -lpolyrt $(LDFLAGS)
		
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

dist_test: orig_par distopt
	touch .test
	OMP_NUM_THREADS=8 ./orig_par > out_orig_par
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt > out_distopt
	rm -f .test
	diff -q out_orig_par out_distopt
	@echo Success!

dist_test2: orig_par distopt
	touch .test
	OMP_NUM_THREADS=$(NPROCS) ./orig_par > out_orig_par
	OMP_NUM_THREADS=1 mpirun -np $(NPROCS) ./distopt > out_distopt
	rm -f .test
	diff -q out_orig_par out_distopt
	@echo Success!



clean:
	rm -f out_* *.tiled.c *.opt.c *.par.c *.distopt.c orig opt tiled par distopt sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog

exec-clean:
	rm -f out_* opt orig tiled  sched sched hopt hopt par orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
