
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#
BASEDIR=$(dir $(lastword $(MAKEFILE_LIST)))

#CC=icc
CC=gcc

MPI=openmpi
#MPI=mvapich
#MPI=intelmpi

#NON_CLUSTER=true

NPROCS=1
NTHREADS=4
NTHREADS_WITH_MPI=2
NTHREADS_LIST = 1 2 4 8 16 32
NPROCS_LIST = 1 2 4 8 16 32
TILE_SIZES = 16 32 64 128 256
POLYRTINCDIR=$(BASEDIR)../polyrt
POLYBENCHINCDIR=$(BASEDIR)polybench/utilities
POLYBENCHSRC=$(BASEDIR)polybench/utilities/polybench.c
HOSTS_FILE=$(BASEDIR)hosts

POLYBENCHINCDIR=$(BASEDIR)polybench/utilities
POLYBENCHSRC=$(BASEDIR)polybench/utilities/polybench.c
PLC=$(BASEDIR)../polycc

# Intel MKL and AMD ACML library paths
MKL=/opt/intel/mkl
ACML=/usr/local/acml

#Metis libray path
METIS=/usr/local/lib/libmetis.a

#Intel scalapack library paths and compiler flags
S_CFLAGS +=  -DMKL_ILP64 -I$(MKLROOT)/include 
S_LFLAGS +=  $(MKLROOT)/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_cdft_core.a $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm  


ifeq ($(CC), icc)
	CXX           := icpc
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -qopenmp
	ifeq ($(MPI), mvapich)
		MPICC        := mpicc -cc=icc -D__MPI
		MPICXX       := mpicc -cc=icpc -D__MPI
	else ifeq ($(MPI), intelmpi)
		MPICC        := mpiicc -D__MPI
		MPICXX       := mpiicpc -D__MPI -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
	else #openmpi
		MPICC        := OMPI_CC=icc mpicc -D__MPI
		MPICXX       := OMPI_CC=icpc mpicc -D__MPI
	endif
else
	# for gcc
	CXX           := g++
	OPT_FLAGS     := -O3 -march=native -mtune=native -ftree-vectorize
	PAR_FLAGS     := -ftree-parallelize-loops=$(NTHREADS)
	OMP_FLAGS     := -fopenmp
	ifeq ($(MPI), mvapich)
		MPICC        := mpicc -cc=gcc -D__MPI
		MPICXX       := mpicc -cc=g++ -D__MPI
	else ifeq ($(MPI), intelmpi)
		MPICC        := mpicc -D__MPI
		MPICXX       := mpicxx -D__MPI
	else #openmpi
		MPICC        := OMPI_CC=gcc mpicc -D__MPI
		MPICXX       := OMPI_CC=g++ mpic++ -D__MPI
	endif
endif

ifeq ($(NON_CLUSTER), true)
	MPI = openmpi
endif

ifeq ($(MPI), mvapich)
	MPIARGS       = -np $(NPROCS) -hostfile $(HOSTS_FILE)
	MPIENV        = MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS_WITH_MPI)
	MPIRUN        = mpirun_rsh $(MPIARGS) $(MPIENV)
else ifeq ($(MPI), intelmpi)
	MPIARGS       = -ppn=1 -np $(NPROCS) -hostfile $(HOSTS_FILE)
	MPIENV        = -env I_MPI_FABRICS ofa -env OMP_NUM_THREADS $(NTHREADS_WITH_MPI)
	MPIRUN        = mpirun $(MPIARGS) $(MPIENV)
else #openmpi
	MPIARGS       = -np $(NPROCS)
	MPIENV        = OMP_NUM_THREADS=$(NTHREADS_WITH_MPI)
	MPIRUN        = $(MPIENV) mpirun $(MPIARGS)
endif

CFLAGS += -DTIME
LDFLAGS += -lm
PLCFLAGS += --isldep --lastwriter --indent
DISTOPT_FLAGS += --cloogsh --timereport
TILEFLAGS += 

#PERFCTR=perfctr

ifdef PERFCTR
	CFLAGS += -DPERFCTR -L/usr/local/lib64 -lpapi
endif

ifdef POLYBENCH
	CFLAGS += -DPOLYBENCH_USE_SCALAR_LB -DPOLYBENCH_TIME -I $(POLYBENCHINCDIR) $(POLYBENCHSRC)
	DISTOPT_FLAGS += --variables_not_global
endif

all: orig tiled par distopt

$(SRC).opt.c:  $(SRC).c
	$(PLC) $(SRC).c --notile --noparallel $(PLCFLAGS)  -o $@

$(SRC).tiled.c:  $(SRC).c
	$(PLC) $(SRC).c --noparallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).idt.c:  $(SRC).c
	$(PLC) $(SRC).c --identity $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).mlbpar.c:  $(SRC).c
	$(PLC) $(SRC).c --full-diamond-tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

# Version that doesn't use diamond tiling
$(SRC).pipepar.c:  $(SRC).c
	$(PLC) $(SRC).c --nodiamond-tile $(TILEFLAGS) $(PLCFLAGS) -o $@

$(SRC).dyn_graph_idt.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --dynschedule_graph --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dyn_graph.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --dynschedule_graph $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dyn_idt.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --dynschedule --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dyn.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --dynschedule $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dyn_data_dist.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --dynschedule --data_dist $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).data_dist.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel --data_dist --identity_data_dist $(TILEFLAGS) $(PLCFLAGS)  $(DISTOPT_FLAGS) -o $@

$(SRC).dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --nocommopt --tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).distopt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_data_dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --mpiomp --tile --data_dist $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_idt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_idnt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_foifi.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_foifi --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dhpf.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_foifi --tile $(TILEFLAGS) --isldep $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_fop.c: $(SRC).c

$(SRC).distopt_fop.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_fop --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_fop_idt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile --commopt_fop --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@
	
$(SRC).distopt_dsfo_data_dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --mpiomp --commopt_foifi --tile --data_dist $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_dsfo_data_dist_verify.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --mpiomp --commopt_foifi --tile --data_dist --verify_output $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_foifi_data_dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --mpiomp --commopt_foifi --tile --data_dist $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_foifi_data_dist_verify.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --mpiomp --commopt_foifi --tile --data_dist --verify_output $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_fop_idnt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_fop --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_foifi_idt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile --commopt_foifi --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distopt_foifi_idnt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --commopt_foifi --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).distomp.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --nocommopt --mpiomp --tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).distrecv.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --timereport --recvpar --tile --mpiomp $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).dist_dynsched.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp  --commopt_fop --dynschedule --tile $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dist_dynsched_data_dist.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp  --commopt_fop --dynschedule --tile --data_dist  --verify_output $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dist_dynsched_idt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile --commopt_fop --dynschedule --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

$(SRC).dist_dynsched_foifi_idt.c: $(SRC).c
	$(PLC) $(SRC).c --distmem --mpiomp --tile --commopt_foifi --dynschedule --identity $(TILEFLAGS) $(PLCFLAGS) $(DISTOPT_FLAGS)  -o $@

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

mlbpar: $(SRC).mlbpar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).mlbpar.c -o $@  $(LDFLAGS)

# Version that doesn't use diamond tiling
pipepar: $(SRC).pipepar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).pipepar.c -o $@  $(LDFLAGS)

parcxx: $(SRC).par.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o $@  $(LDFLAGS) -ltbb

idt: $(SRC).idt.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).idt.c -o $@  $(LDFLAGS)

idtcxx: $(SRC).idt.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).idt.c -o $@  $(LDFLAGS) -ltbb

data_dist: $(SRC).data_dist.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).data_dist.c -D__DATA_DIST_DECLS $(POLYRTINCDIR)/buffer_manager.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb

dyn_graph_idt: $(SRC).dyn_graph_idt.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).dyn_graph_idt.c sigma_$(SRC).dyn_graph_idt.c -o $@ $(LDFLAGS) -ltbb

dyn_graph: $(SRC).dyn_graph.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).dyn_graph.c sigma_$(SRC).dyn_graph.c -o $@ $(LDFLAGS) -ltbb

dyn_nopr_idt: $(SRC).dyn_idt.c sigma_$(SRC).dyn_idt.c pi_$(SRC).dyn_idt.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) -D__DYNSCHEDULER_NO_PRIORITY $(OMP_FLAGS) $(SRC).dyn_idt.c sigma_$(SRC).dyn_idt.c pi_$(SRC).dyn_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb

dyn_idt: $(SRC).dyn_idt.c sigma_$(SRC).dyn_idt.c pi_$(SRC).dyn_idt.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).dyn_idt.c sigma_$(SRC).dyn_idt.c pi_$(SRC).dyn_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb

dyn_nopr: $(SRC).dyn.c sigma_$(SRC).dyn.c pi_$(SRC).dyn.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) -D__DYNSCHEDULER_NO_PRIORITY $(OMP_FLAGS) $(SRC).dyn.c sigma_$(SRC).dyn.c pi_$(SRC).dyn.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb

dyn: $(SRC).dyn.c sigma_$(SRC).dyn.c pi_$(SRC).dyn.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).dyn.c sigma_$(SRC).dyn.c pi_$(SRC).dyn.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb 

dyn_data_dist: $(SRC).dyn_data_dist.c sigma_$(SRC).dyn_data_dist.c pi_$(SRC).dyn_data_dist.c
	$(CXX) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) -D__DATA_DIST_DECLS $(SRC).dyn_data_dist.c sigma_$(SRC).dyn_data_dist.c pi_$(SRC).dyn_data_dist.c \
		$(POLYRTINCDIR)/polyrt.c $(POLYRTINCDIR)/buffer_manager.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS) -ltbb

dist: $(SRC).dist.c pi_$(SRC).dist.c sigma_$(SRC).dist.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).dist.c pi_$(SRC).dist.c sigma_$(SRC).dist.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt: $(SRC).distopt.c sigma_$(SRC).distopt.c pi_$(SRC).distopt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt.c sigma_$(SRC).distopt.c pi_$(SRC).distopt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_data_dist: $(SRC).distopt_data_dist.c  sigma_$(SRC).distopt_data_dist.c pi_$(SRC).distopt_data_dist.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS)  -D__DATA_DIST_DECLS $(CFLAGS) $(SRC).distopt_data_dist.c sigma_$(SRC).distopt_data_dist.c pi_$(SRC).distopt_data_dist.c \
		$(POLYRTINCDIR)/polyrt.c $(POLYRTINCDIR)/buffer_manager.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)
		
distopt_idt: $(SRC).distopt_idt.c sigma.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI $(CFLAGS) $(SRC).distopt_idt.c sigma.c pi.c \
		-o $@ -L $(POLYRTLIBDIR) -lpolyrt $(LDFLAGS)
		
distopt_idnt: $(SRC).distopt_idnt.c sigma_$(SRC).distopt_idnt.c pi_$(SRC).distopt_idnt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_idnt.c sigma_$(SRC).distopt_idnt.c pi_$(SRC).distopt_idnt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_foifi: $(SRC).distopt_foifi.c sigma_$(SRC).distopt_foifi.c pi_$(SRC).distopt_foifi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_foifi.c sigma_$(SRC).distopt_foifi.c pi_$(SRC).distopt_foifi.c\
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)
		
dhpf: $(SRC).dhpf.c sigma_$(SRC).dhpf.c pi_$(SRC).dhpf.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).dhpf.c sigma_$(SRC).dhpf.c pi_$(SRC).dhpf.c\
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)
		
distopt_fop_idt: $(SRC).distopt_fop_idt.c sigma_$(SRC).distopt_fop_idt.c pi_$(SRC).distopt_fop_idt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop_idt.c sigma_$(SRC).distopt_fop_idt.c pi_$(SRC).distopt_fop_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_fop_idnt: $(SRC).distopt_fop_idnt.c sigma_$(SRC).distopt_fop_idnt.c pi_$(SRC).distopt_fop_idnt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop_idnt.c sigma_$(SRC).distopt_fop_idnt.c pi_$(SRC).distopt_fop_idnt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)
		
distopt_foifi_idt: $(SRC).distopt_foifi_idt.c sigma_$(SRC).distopt_foifi_idt.c pi_$(SRC).distopt_foifi_idt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_foifi_idt.c sigma_$(SRC).distopt_foifi_idt.c pi_$(SRC).distopt_foifi_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)
		
distopt_foifi_idnt: $(SRC).distopt_foifi_idnt.c sigma_$(SRC).distopt_foifi_idnt.c pi_$(SRC).distopt_foifi_idnt.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_foifi_idnt.c sigma_$(SRC).distopt_foifi_idnt.c pi_$(SRC).distopt_foifi_idnt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distopt_fop: $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c\
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)
	
distopt_dsfo_data_dist: $(SRC).distopt_dsfo_data_dist.c sigma_dsfo.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI -DUSE_LOCAL_ARRAYS  $(CFLAGS) $(SRC).distopt_dsfo_data_dist.c sigma_dsfo.c pi.c\
        -o $@ -L $(POLYRTLIBDIR) -I $(POLYRTINCDIR)  -lpolyrt $(LDFLAGS)
        	

distrecv: $(SRC).distrecv.c sigma_$(SRC).distrecv.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distrecv.c sigma_$(SRC).distrecv.c pi.c\
		-o $@ -I $(POLYRTINCDIR) -L $(POLYRTLIBDIR) -lpolyrt $(LDFLAGS)
		
distopt_dsfo_data_dist_verify: $(SRC).distopt_dsfo_data_dist_verify.c sigma_dsfo.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI   $(CFLAGS) $(SRC).distopt_dsfo_data_dist_verify.c sigma_dsfo.c pi.c\
		-o $@ -L $(POLYRTLIBDIR) -I $(POLYRTINCDIR)  -lpolyrt $(LDFLAGS)
		
distopt_foifi_data_dist: $(SRC).distopt_foifi_data_dist.c sigma_dsfo.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI -DUSE_LOCAL_ARRAYS  $(CFLAGS) $(SRC).distopt_foifi_data_dist.c sigma_dsfo.c pi.c\
		-o $@ -L $(POLYRTLIBDIR) -I $(POLYRTINCDIR)  -lpolyrt $(LDFLAGS)
		
distopt_foifi_data_dist_verify: $(SRC).distopt_foifi_data_dist_verify.c sigma_dsfo.c pi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) -DMPI   $(CFLAGS) $(SRC).distopt_foifi_data_dist_verify.c sigma_dsfo.c pi.c\
		-o $@ -L $(POLYRTLIBDIR) -I $(POLYRTINCDIR)  -lpolyrt $(LDFLAGS)
	

mpi: $(SRC).mpi.c sigma_$(SRC).mpi.c pi_$(SRC).mpi.c
	$(MPICC) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).mpi.c sigma_$(SRC).mpi.c pi_$(SRC).mpi.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) $(LDFLAGS)

distcxx: $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop.c sigma_$(SRC).distopt_fop.c pi_$(SRC).distopt_fop.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)
		
distcxx_idt: $(SRC).distopt_fop_idt.c sigma_$(SRC).distopt_fop_idt.c pi_$(SRC).distopt_fop_idt.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).distopt_fop_idt.c sigma_$(SRC).distopt_fop_idt.c pi_$(SRC).distopt_fop_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)

dist_dynsched_nopr: $(SRC).dist_dynsched.c sigma_$(SRC).dist_dynsched.c pi_$(SRC).dist_dynsched.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) -D__DYNSCHEDULER_NO_PRIORITY $(SRC).dist_dynsched.c sigma_$(SRC).dist_dynsched.c pi_$(SRC).dist_dynsched.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)
		
dist_dynsched_nopr_idt: $(SRC).dist_dynsched_idt.c sigma_$(SRC).dist_dynsched_idt.c pi_$(SRC).dist_dynsched_idt.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) -D__DYNSCHEDULER_NO_PRIORITY $(SRC).dist_dynsched_idt.c sigma_$(SRC).dist_dynsched_idt.c pi_$(SRC).dist_dynsched_idt.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)

dist_dynsched: $(SRC).dist_dynsched.c sigma_$(SRC).dist_dynsched.c pi_$(SRC).dist_dynsched.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).dist_dynsched.c sigma_$(SRC).dist_dynsched.c pi_$(SRC).dist_dynsched.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS) 

dist_dynsched_data_dist: $(SRC).dist_dynsched_data_dist.c sigma_$(SRC).dist_dynsched_data_dist.c pi_$(SRC).dist_dynsched_data_dist.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) -D__DATA_DIST_DECLS $(CFLAGS) $(SRC).dist_dynsched_data_dist.c sigma_$(SRC).dist_dynsched_data_dist.c pi_$(SRC).dist_dynsched_data_dist.c \
		$(POLYRTINCDIR)/polyrt.c $(POLYRTINCDIR)/buffer_manager.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)

dist_dynsched_data_dist_local: $(SRC).dist_dynsched_data_dist.c sigma_$(SRC).dist_dynsched_data_dist.c pi_$(SRC).dist_dynsched_data_dist.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS)  $(SRC).dist_dynsched_data_dist.c sigma_$(SRC).dist_dynsched_data_dist.c pi_$(SRC).dist_dynsched_data_dist.c \
		$(POLYRTINCDIR)/polyrt.c $(POLYRTINCDIR)/buffer_manager.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS) 

dist_dynsched_idt: $(SRC).dist_dynsched_foifi_idt.c sigma_$(SRC).dist_dynsched_foifi_idt.c pi.c
	$(MPICXX) $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(SRC).dist_dynsched_foifi_idt.c sigma_$(SRC).dist_dynsched_foifi_idt.c pi.c \
		$(POLYRTINCDIR)/polyrt.c -o $@ -I $(POLYRTINCDIR) -ltbb $(LDFLAGS)

scalapack: $(SRC).scalapack.c 
	$(MPICC) -DMPI  $(OPT_FLAGS) $(OMP_FLAGS) $(CFLAGS) $(S_CFLAGS) $(SRC).scalapack.c -o scalapack  $(S_LFLAGS) 

scalapack_run: scalapack 
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./scalapack  2>out_scalapack
	touch .test
	rm -f .test

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

distperf: par dist distomp distopt
	rm -f .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 
	mpirun -np 4 ./dist 
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun -np 4 ./distomp
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun -np 4 ./distopt 


test: orig tiled par
	touch .test
	./orig 2> out_orig
	./tiled 2> out_tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	rm -f .test
	diff -q out_orig out_tiled
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

data_dist_test: par data_dist
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	OMP_NUM_THREADS=$(NTHREADS) ./data_dist 2> out_data_dist4
	diff -q out_par4 out_data_dist4
	rm -f .test
	@echo Success!

dyn_run: dyn
	OMP_NUM_THREADS=$(NTHREADS) ./dyn 2> out_dyn

	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./parcxx 2> out_par
	OMP_NUM_THREADS=$(NTHREADS) ./dyn 2> out_dyn
	rm -f .test
	diff -q out_par out_dyn
	@echo Success!

dyn_data_dist_test: par dyn_data_dist
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_data_dist 2> out_dyn_data_dist
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par
	diff -q out_par out_dyn_data_dist
	rm -f .test
	@echo Success!


	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./idtcxx 2> out_idt
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_graph_idt 2> out_dyn_graph_idt
	rm -f .test
	diff -q out_idt out_dyn_graph_idt
	@echo Success!

dyn_graph_test: parcxx dyn_graph
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./parcxx 2> out_par
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_graph 2> out_dyn_graph
	diff -q out_par out_dyn_graph
	rm -f .test
	@echo Success!

dyn_compare_idt_test: idtcxx dyn_graph_idt dyn_nopr_idt dyn_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./idtcxx 2> out_idt
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_graph_idt 2> out_dyn_graph_idt
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_nopr_idt 2> out_dyn_nopr_idt
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_idt 2> out_dyn_idt
	diff -q out_idt out_dyn_graph_idt
	diff -q out_idt out_dyn_nopr_idt
	diff -q out_idt out_dyn_idt
	rm -f .test
	@echo Success!

dyn_compare_test: parcxx dyn_graph dyn_nopr dyn
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./parcxx 2> out_par
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_graph 2> out_dyn_graph
	OMP_NUM_THREADS=$(NTHREADS) ./dyn_nopr 2> out_dyn_nopr
	OMP_NUM_THREADS=$(NTHREADS) ./dyn 2> out_dyn
	diff -q out_par out_dyn_graph
	diff -q out_par out_dyn_nopr
	diff -q out_par out_dyn
	rm -f .test
	@echo Success!

# does not verify results
dyn_compare_idt_threads: idtcxx dyn_graph_idt dyn_nopr_idt dyn_idt
	$(foreach t, $(NTHREADS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running idtcxx for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./idtcxx 2> out_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn_graph_idt for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn_graph_idt 2> out_dyn_graph_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn_nopr_idt for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn_nopr_idt 2> out_dyn_nopr_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn_idt 2> out_dyn_idt;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "" >&1;\
		echo "" >&1;\
		echo "" >&1;)

# does not verify results
dyn_compare_threads: parcxx dyn_graph dyn_nopr dyn
	$(foreach t, $(NTHREADS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running parcxx for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./parcxx 2> out_par;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn_graph for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn_graph 2> out_dyn_graph;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn_nopr for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn_nopr 2> out_dyn_nopr;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running dyn for $(SRC) with $(t) threads" >&1;\
		echo "" >&1; \
		OMP_NUM_THREADS=$(t) ./dyn 2> out_dyn;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "" >&1;\
		echo "" >&1;\
		echo "" >&1;)

opt-test: orig opt
	touch .test
	./orig 2> out_orig
	./opt 2> out_opt
	rm -f .test
	diff -q out_orig out_opt
	@echo Success!

dist_test: orig_par distopt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	$(MPIRUN) ./distopt 2> out_distopt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	diff -q out_orig_par out_distopt
	@echo Success!


dist_idt_test: orig_par distopt_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	$(MPIRUN) ./distopt_idt 2> out_distopt_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	diff -q out_orig_par out_distopt_idt
	@echo Success!

dist_idt_test2: orig_par distopt_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_idt 2> out_distopt_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	diff -q out_orig_par out_distopt_idt
	@echo Success!

dist_idnt_test: orig_par distopt_idnt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	$(MPIRUN) ./distopt_idnt 2> out_distopt_idnt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idnt
	diff -q out_orig_par out_distopt_idnt
	@echo Success!

dist_dhpf_test: distopt dhpf
	touch .test
	$(MPIRUN) ./distopt 2> out_distopt
	$(MPIRUN) ./dhpf 2> out_dhpf
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_dhpf
	diff -q out_distopt out_dhpf
	@echo Success!

dist_foifi_test: distopt distopt_foifi
	touch .test
	$(MPIRUN) ./distopt 2> out_distopt
	$(MPIRUN) ./distopt_foifi 2> out_distopt_foifi
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_foifi
	diff -q out_distopt out_distopt_foifi
	@echo Success!

dist_fop_idt_test: distopt_idt distopt_fop_idt
	touch .test
	$(MPIRUN) ./distopt_idt 2> out_distopt_idt
	$(MPIRUN) ./distopt_fop_idt 2> out_distopt_fop_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_fop_idt
	diff -q out_distopt_idt out_distopt_fop_idt
	@echo Success!

dist_fop_idt_test2: distopt_idt distopt_fop_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_idt 2> out_distopt_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_fop_idt 2> out_distopt_fop_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_fop_idt
	diff -q out_distopt_idt out_distopt_fop_idt
	@echo Success!

dist_foifi_idt_test: distopt_idt distopt_foifi_idt
	touch .test
	$(MPIRUN) ./distopt_idt 2> out_distopt_idt
	$(MPIRUN) ./distopt_foifi_idt 2> out_distopt_foifi_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_foifi_idt
	diff -q out_distopt_idt out_distopt_foifi_idt
	@echo Success!

dist_foifi_idt_test2: distopt_idt distopt_foifi_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_idt 2> out_distopt_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_foifi_idt 2> out_distopt_foifi_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_foifi_idt
	diff -q out_distopt_idt out_distopt_foifi_idt
	@echo Success!

dist_compare_idt_test: distopt_idt distopt_fop_idt distopt_foifi_idt
	touch .test
	$(MPIRUN) ./distopt_idt 2> out_distopt_idt
	$(MPIRUN) ./distopt_fop_idt 2> out_distopt_fop_idt
	$(MPIRUN) ./distopt_foifi_idt 2> out_distopt_foifi_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_fop_idt
	sed -i '/librdmacm/d' out_distopt_foifi_idt
	diff -q out_distopt_idt out_distopt_fop_idt
	diff -q out_distopt_idt out_distopt_foifi_idt
	@echo Success!

dist_compare_idt_test2: distopt_idt distopt_fop_idt distopt_foifi_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_idt 2> out_distopt_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_fop_idt 2> out_distopt_fop_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_foifi_idt 2> out_distopt_foifi_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_idt
	sed -i '/librdmacm/d' out_distopt_fop_idt
	sed -i '/librdmacm/d' out_distopt_foifi_idt
	diff -q out_distopt_idt out_distopt_fop_idt
	diff -q out_distopt_idt out_distopt_foifi_idt
	@echo Success!

dist_fop_idnt_test: orig_par distopt_fop_idnt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	$(MPIRUN) ./distopt_fop_idnt 2> out_distopt_fop_idnt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_fop_idnt
	diff -q out_orig_par out_distopt_fop_idnt
	@echo Success!

dist_fop_test: distopt distopt_fop
	touch .test
	$(MPIRUN) ./distopt 2> out_distopt
	$(MPIRUN) ./distopt_fop 2> out_distopt_fop
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_fop
	diff -q out_distopt out_distopt_fop
	@echo Success!

dist_dsfo_data_test: orig_par distopt_dsfo_data_dist_verify
	touch .test
	OMP_NUM_THREADS=1 mpirun  -np $(NPROCS) ./distopt_dsfo_data_dist_verify 2> out_distopt_dsfo_data_dist
	OMP_NUM_THREADS=8 ./orig_par 2> out_orig_par
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_dsfo_data_dist
	diff -q out_distopt_dsfo_data_dist out_orig_par 
	@echo Success!

dist_foifi_data_test: orig_par distopt_foifi_data_dist_verify
	touch .test
	OMP_NUM_THREADS=1 mpirun  -np $(NPROCS) ./distopt_foifi_data_dist_verify 2> out_distopt_foifi_data_dist
	OMP_NUM_THREADS=8 ./orig_par 2> out_orig_par
	rm -f .test
	sed -i '/librdmacm/d' out_distopt_foifi_data_dist
	diff -q out_distopt_foifi_data_dist out_orig_par 
	@echo Success!


dist_compare_test: distopt distopt_foifi distopt_fop
	touch .test
	$(MPIRUN) ./distopt 2> out_distopt
	$(MPIRUN) ./distopt_fop 2> out_distopt_fop
	$(MPIRUN) ./distopt_foifi 2> out_distopt_foifi
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_fop
	sed -i '/librdmacm/d' out_distopt_foifi
	diff -q out_distopt out_distopt_fop
	diff -q out_distopt out_distopt_foifi
	@echo Success!

dist_compare_test2: distopt distopt_foifi distopt_fop
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt 2> out_distopt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_fop 2> out_distopt_fop
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distopt_foifi 2> out_distopt_foifi
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_fop
	sed -i '/librdmacm/d' out_distopt_foifi
	diff -q out_distopt out_distopt_fop
	diff -q out_distopt out_distopt_foifi
	@echo Success!

dist_dynsched_data_dist_test: dist_dynsched dist_dynsched_data_dist
	touch .test
	$(MPIRUN) ./dist_dynsched 2> out_dist_dynsched
	$(MPIRUN) ./dist_dynsched_data_dist 2> out_dist_dynsched_data_dist
	rm -f .test
	sed -i '/librdmacm/d' out_dist_dynsched
	sed -i '/librdmacm/d' out_dist_dynsched_data_dist
	diff -q out_dist_dynsched out_dist_dynsched_data_dist
	@echo Success!

dist_dynsched_data_dist_run: dist_dynsched_data_dist_local orig
	touch .test
	$(MPIRUN) ./dist_dynsched_data_dist_local 2> out_dist_dynsched_data_dist_run
	$(MPIRUN) ./dist_dynsched 2> out_dist_dynsched
	./orig 2> out_orig
	rm -f .test
	sed -i '/librdmacm/d' out_dist_dynsched
	sed -i '/librdmacm/d' out_dist_dynsched_data_dist_run
	diff -q out_dist_dynsched out_orig
	diff -q out_dist_dynsched out_dist_dynsched_data_dist_run
	@echo Success!

dist_dynsched_test: distcxx dist_dynsched_nopr dist_dynsched
	touch .test
	$(MPIRUN) ./dist_dynsched 2> out_dist_dynsched
	$(MPIRUN) ./distcxx 2> out_distcxx
	$(MPIRUN) ./dist_dynsched_nopr 2> out_dist_dynsched_nopr
	rm -f .test
	sed -i '/librdmacm/d' out_distcxx
	sed -i '/librdmacm/d' out_dist_dynsched_nopr
	sed -i '/librdmacm/d' out_dist_dynsched
	diff -q out_distcxx out_dist_dynsched_nopr
	diff -q out_distcxx out_dist_dynsched
	@echo Success!

dist_dynsched_test2: distcxx dist_dynsched_nopr dist_dynsched
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distcxx 2> out_distcxx
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./dist_dynsched_nopr 2> out_dist_dynsched_nopr
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./dist_dynsched 2> out_dist_dynsched
	rm -f .test
	sed -i '/librdmacm/d' out_distcxx
	sed -i '/librdmacm/d' out_dist_dynsched_nopr
	sed -i '/librdmacm/d' out_dist_dynsched
	diff -q out_distcxx out_dist_dynsched_nopr
	diff -q out_distcxx out_dist_dynsched
	@echo Success!

dist_dynsched_idt_test: distcxx_idt dist_dynsched_nopr_idt dist_dynsched_idt
	touch .test
	$(MPIRUN) ./distcxx_idt 2> out_distcxx_idt
	$(MPIRUN) ./dist_dynsched_nopr_idt 2> out_dist_dynsched_nopr_idt
	$(MPIRUN) ./dist_dynsched_idt 2> out_dist_dynsched_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distcxx_idt
	sed -i '/librdmacm/d' out_dist_dynsched_nopr_idt
	sed -i '/librdmacm/d' out_dist_dynsched_idt
	diff -q out_distcxx_idt out_dist_dynsched_nopr_idt
	diff -q out_distcxx_idt out_dist_dynsched_idt
	@echo Success!

dist_dynsched_idt_test2: distcxx_idt dist_dynsched_nopr_idt dist_dynsched_idt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./distcxx_idt 2> out_distcxx_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./dist_dynsched_nopr_idt 2> out_dist_dynsched_nopr_idt
	OMP_NUM_THREADS=$(NTHREADS_WITH_MPI) mpirun  -np $(NPROCS) ./dist_dynsched_idt 2> out_dist_dynsched_idt
	rm -f .test
	sed -i '/librdmacm/d' out_distcxx_idt
	sed -i '/librdmacm/d' out_dist_dynsched_nopr_idt
	sed -i '/librdmacm/d' out_dist_dynsched_idt
	diff -q out_distcxx_idt out_dist_dynsched_nopr_idt
	diff -q out_distcxx_idt out_dist_dynsched_idt
	@echo Success!

# does not verify results
dist_dynsched_compare_procs: distcxx dist_dynsched_nopr dist_dynsched
	$(foreach NPROCS, $(NPROCS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running distcxx for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distcxx 2> out_distcxx;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running dist_dynsched_nopr for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./dist_dynsched_nopr 2> out_dist_dynsched_nopr;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running dist_dynsched for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./dist_dynsched 2> out_dist_dynsched;\
		echo "" >&1;)

# does not verify results
dist_dynsched_compare_idt_procs: distcxx_idt dist_dynsched_nopr_idt dist_dynsched_idt
	$(foreach NPROCS, $(NPROCS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running distcxx_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distcxx_idt 2> out_distcxx_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running dist_dynsched_nopr_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./dist_dynsched_nopr_idt 2> out_dist_dynsched_nopr_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running dist_dynsched_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./dist_dynsched_idt 2> out_dist_dynsched_idt;\
		echo "" >&1;)

# does not verify results
dist_compare_procs: distopt distopt_fop distopt_foifi
	$(foreach NPROCS, $(NPROCS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running distopt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt 2> out_distopt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running distopt_fop for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt_fop 2> out_distopt_fop;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running distopt_foifi for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt_foifi 2> out_distopt_foifi;\
		echo "" >&1;)

    
data_dist_compare_scalapack: distopt_data_dist 
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt_data_dist 2> out_distopt_data_dist
	echo "Running distopt data dist for $(SRC) with $(p) procs" >&1
	echo "" >&1 
	make scalapack_run

# does not verify results
data_dist_compare_procs: distopt_data_dist 
	$(foreach p, $(NPROCS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running distrpt data dist for $(SRC) with $(p) procs" >&1;\
		echo "" >&1; \
		mpirun_rsh  -np $(p) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt_data_dist 2> out_distopt_data_dist;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "" >&1;)


# does not verify results
# assumes utmost 6-dimensional tiling
data_dist_compare_tilesizes: 
	$(foreach t, $(TILE_SIZES), \
		echo -e "$(t)\n$(t)\n$(t)\n$(t)\n$(t)\n$(t)" > tile.sizes;\
		make clean;\
		echo "Tile size = $(t)" >&1;\
		make data_dist_compare_procs;\
		)


# does not verify results
# assumes utmost 6-dimensional tiling
dist_compare_tilesizes: 
	$(foreach t, $(TILE_SIZES), \
		echo "$(t)" > tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		make clean;\
		make distopt;\
		make distopt_fop;\
		make distopt_foifi;\
		echo "Tile size = $(t)" >&1;\
		make dist_compare_procs;\
		)

# does not verify results
dist_compare_idt_procs: distopt_idt distopt_fop_idt distopt_foifi_idt
	$(foreach NPROCS, $(NPROCS_LIST),  \
		echo "-----------------------------------------------------------------------------------------------------" >&1;  \
		echo "Running distopt_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt_idt 2> out_distopt_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running distopt_fop_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt_fop_idt 2> out_distopt_fop_idt;\
		echo "" >&1;\
		echo "-----------------------------------------------------------------------------------------------------" >&1;\
		echo "Running distopt_foifi_idt for $(SRC) with $(NPROCS) procs" >&1;\
		echo "" >&1; \
		$(MPIRUN) ./distopt_foifi_idt 2> out_distopt_foifi_idt;\
		echo "" >&1;)

# does not verify results
# assumes utmost 6-dimensional tiling
dist_compare_idt_tilesizes: 
	$(foreach t, $(TILE_SIZES), \
		echo "$(t)" > tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		echo "$(t)" >> tile.sizes;\
		make clean;\
		make distopt_idt;\
		make distopt_fop_idt;\
		make distopt_foifi_idt;\
		echo "Tile size = $(t)" >&1;\
		make dist_compare_idt_procs;\
		)

dist_test_seq: orig_par distopt
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par 2> out_orig_par
	./distopt 2> out_distopt
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	diff -q out_orig_par out_distopt
	@echo Success!

dist_data_test: distopt distopt_data_dist
	touch .test
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt 2> out_distopt
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt_data_dist 2> out_distopt_data_dist
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_data_dist
	diff -q out_distopt out_distopt_data_dist
	@echo Success!

dist_data_test2: distopt distopt_data_dist
	touch .test
	OMP_NUM_THREADS=1 mpirun -np $(NPROCS) ./distopt 2> out_distopt
	OMP_NUM_THREADS=1 mpirun -np $(NPROCS) ./distopt_data_dist 2> out_distopt_data_dist
	rm -f .test
	sed -i '/librdmacm/d' out_distopt
	sed -i '/librdmacm/d' out_distopt_data_dist
	diff -q out_distopt out_distopt_data_dist
	@echo Success!

dist_data_run: distopt_data_dist
	mpirun_rsh  -np $(NPROCS) -hostfile $(HOSTS_FILE) MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=$(NTHREADS) ./distopt_data_dist
	@echo Success!



clean:
	rm -f out_* *.pipepar.* *.tiled.* *.opt.* *.idt.* *.par.* *.mlbpar.c *.dyn_graph_idt.* *.dyn_graph.* *.dyn_idt.* *.dyn.* *.*data_dist.* \
		*.dist.c *.distomp.c *.distopt.c *.distrecv.c *.distopt_foifi.c *.distopt_fop.c *.dhpf.c *_idt.c *_idnt.c \
		*.dist_dynsched*.c \
	   	orig opt tiled idt idtcxx par parcxx dyn_graph_idt dyn_graph dyn_graph_lib dyn_nopr_idt dyn_idt dyn_nopr dyn dist sched orig_par distomp dhpf lbpar \
		distopt distopt_fop distopt_foifi distopt_idt distopt_fop_idt distopt_foifi_idt distopt_idnt distopt_fop_idnt distopt_foifi_idnt \
		distcxx distcxx_idt dist_dynsched_nopr dist_dynsched_nopr_idt dist_dynsched dist_dynsched_idt \
		hopt hopt *.par2d.c *.out.* *.dist*.h \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
		.distmem .srcfilename .outfilename .vectorize par2d parsetab.py *.body.c *.pluto.c \
		*.par.cloog *.tiled.cloog *.pluto.cloog sigma.cloog is_receiver.cloog sigma_fop.cloog sigma_check_fop.cloog is_receiver_fop.cloog packunpack.cloog \
		count_remote_dep_tasks.cloog remote_update_dep_tasks.cloog remote_count_dep_tasks.cloog \
		count_local_dep_tasks.cloog local_update_dep_tasks.cloog add_outgoing_edges.cloog local_init_remote_dep_tasks.cloog \
		count_remote_src_tasks.cloog count_local_src_tasks.cloog count_sending_tasks.cloog \
		init_tasks.cloog write_out.cloog compute_task.cloog \
		pi_*.c sigma.c sigma_*.c packunpack.c pi1.c tau.c pi_defs.h *.append.c .appendfilename *.cloog debug_print_node*

exec-clean:
	rm -f out_* opt orig tiled sched sched hopt hopt idt idtcxx par parcxx dyn_graph_idt dyn_graph dyn_graph_lib dyn_nopr_idt dyn_idt dyn_nopr dyn pipepar orig_par dist distomp \
		distopt distopt_fop distopt_foifi distopt_idt distopt_fop_idt distopt_foifi_idt distopt_idnt distopt_fop_idnt distopt_foifi_idnt \
		distcxx distcxx_idt dist_dynsched_nopr dist_dynsched_nopr_idt dist_dynsched dist_dynsched_idt \
		dhpf *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
