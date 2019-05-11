# Pluto
#
# Copyright (C) Uday Bondhugula 2007
#
SUBDIRS = 

MAINTAINERCLEANFILES = Makefile.in

if EXTERNAL_ISL
# ISL_INCLUDE, ISL_LIBADD were already defined
else
ISL_INCLUDE += -I$(top_srcdir)/isl/include
ISL_LIBADD += $(top_srcdir)/isl/libisl.la
endif

bin_PROGRAMS = pluto

pluto_SOURCES = \
		ast_transform.c \
		constraints.c \
		constraints.h \
	       	constraints_polylib.c \
		constraints_isl.c \
		ddg.c \
		distmem.c \
		framework.c \
		framework-dfp.c \
		gurobi-plugin.c \
		iss.c \
		main.c \
		math_support.c \
		math_support.h \
		pluto.c \
		pluto_codegen_if.c \
		pluto.h \
		polyflow.c \
		polyloop.c \
		post_transform.c \
		program.cpp \
		program.h \
		tile.c \
		transforms.c \
		unit_tests.c 

if PLUTO_DEBUG
OPT_FLAGS = 
DEBUG_FLAGS = -g -O0
else
OPT_FLAGS = -O3 -funroll-loops -fopenmp
DEBUG_FLAGS =
endif

pluto_CXXFLAGS = $(OPT_FLAGS) $(DEBUG_FLAGS) \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP -DPLUTO_OPENCL \
   -I../include \
   -I../piplib/include \
   -I../clan/include \
   -I../pet/include \
   $(ISL_INCLUDE) \
   -I../cloog-isl/include \
   -I../openscop/include \
   -I../candl/include \
   -I../polylib/include
# To avoid same name object from being created with libtool and without it.
pluto_CFLAGS = $(AM_CFLAGS) $(pluto_CXXFLAGS)

pluto_LDFLAGS = -static -fopenmp
pluto_LDADD =  \
   ../piplib/libpiplib_dp.la \
   ../openscop/libosl.la \
   ../clan/libclan.la \
   ../cloog-isl/libcloog-isl.la \
   ../pet/libpet.la \
   $(ISL_LIBADD) \
   ../polylib/libpolylib64.la \
   ../candl/libcandl.la \
   -lm

lib_LTLIBRARIES = libpluto.la

# libpluto doesn't use polylib-based functions
libpluto_la_SOURCES = \
		      constraints.c \
		      constraints.h \
		      constraints_isl.c \
		      ddg.c \
		      framework.c \
		      framework-dfp.c \
		      iss.c \
		      libpluto.c \
		      libpluto_dummy.c \
		      math_support.c \
		      math_support.h \
		      pluto.c \
		      pluto.h \
		      polyloop.c \
		      post_transform.c \
		      program.cpp \
		      program.h \
		      tile.c \
		      transforms.c

libpluto_la_CXXFLAGS = $(OPT_FLAGS) $(DEBUG_FLAGS) \
   -DLIB_PLUTO -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../include \
   -I../piplib/include \
   -I../clan/include \
   -I../pet/include \
   -I../openscop/include \
   -I../candl/include \
   -I../cloog-isl/include \
   $(ISL_INCLUDE)
# To avoid same name object from being created with libtool and without it.
libpluto_la_CFLAGS = $(AM_CFLAGS) $(libpluto_la_CXXFLAGS)

# libpluto doesn't use polylib
libpluto_la_LIBADD = \
   ../piplib/libpiplib_dp.la \
   ../pet/libpet.la \
   ../openscop/libosl.la \
   ../clan/libclan.la \
   ../candl/libcandl.la \
   $(ISL_LIBADD) \
   -lm -lgomp
libpluto_la_LDFLAGS = -version-info 1:0:1
