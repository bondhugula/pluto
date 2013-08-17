SUBDIRS = .

ACLOCAL_AMFLAGS = -I m4

#MPI_INCLUDE=/usr/lib/openmpi/include/
#MPI_LIB=/usr/lib/openmpi/lib
#MPI_INCLUDE=/usr/mpi/intel/mvapich2-1.4/include
#MPI_LIB=/usr/mpi/intel/mvapich2-1.4/lib
#MPI_INCLUDE=/usr/include/mvapich2-x86_64/include
#MPI_LIB=/usr/lib64/mvapich2/lib
MPI_INCLUDE=/usr/include/openmpi-x86_64
MPI_LIB=/usr/lib64/openmpi/lib

MAINTAINERCLEANFILES = Makefile.in

lib_LTLIBRARIES = libpolyrt.la libdynscheduler.la libstatscheduler.la

libpolyrt_la_SOURCES = polyrt.c
libdynscheduler_la_SOURCES = scheduler.c
libdynscheduler_la_CFLAGS = -D__USE_DAG_SCHEDULER_AS_LIB__ #-D__DEBUG1__ -D__DEBUG2__
libdynscheduler_la_CFLAGS += -fopenmp
libstatscheduler_la_SOURCES = scheduler.c
libstatscheduler_la_CFLAGS = -D__USE_DAG_SCHEDULER_AS_LIB__ -D__STATIC_SCHEDULE__ #-D__DEBUG1__ -D__DEBUG2__
libstatscheduler_la_CFLAGS += -fopenmp

libpolyrt_la_LIBADD = -L$(MPI_LIB)
INCLUDES = -I$(MPI_INCLUDE)
