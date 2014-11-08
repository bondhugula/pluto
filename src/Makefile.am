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

pluto_SOURCES = constraints.c constraints_polylib.c constraints_isl.c math_support.c ddg.c transforms.c pluto.c post_transform.c framework.c pluto_codegen_if.c polyloop.c ast_transform.c iss.c main.c constraints.h math_support.h pluto.h program.c program.h tile.c unit_tests.c

if PLUTO_DEBUG
OPT_FLAGS = 
DEBUG_FLAGS = -g
else
OPT_FLAGS = -O3 -funroll-loops -fopenmp
endif

pluto_CPPFLAGS = $(OPT_FLAGS) $(DEBUG_FLAGS) -funroll-loops -fopenmp \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../include \
   -I../piplib/include \
   -I../clan/include \
   $(ISL_INCLUDE) \
   -I../cloog-isl/include \
   -I../openscop/include \
   -I../candl/include \
   -I../polylib/include

pluto_LDFLAGS = -static -fopenmp
pluto_LDADD =  \
   ../piplib/libpiplib64.la \
   ../openscop/libosl.la \
   ../clan/libclan.la \
   ../cloog-isl/libcloog-isl.la \
   $(ISL_LIBADD) \
   ../polylib/libpolylib64.la \
   ../candl/libcandl.la \
   -lm

lib_LTLIBRARIES = libpluto.la

# libpluto doesn't use polylib-based functions
libpluto_la_SOURCES = constraints.c constraints_isl.c math_support.c ddg.c transforms.c pluto.c post_transform.c framework.c polyloop.c iss.c libpluto.c constraints.h math_support.h pluto.h program.c program.h tile.c
libpluto_la_CPPFLAGS = $(OPT_FLAGS) $(DEBUG_FLAGS) \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../include \
   -I../piplib/include \
   -I../clan/include \
   -I../openscop/include \
   -I../candl/include \
   -I../cloog-isl/include \
   $(ISL_INCLUDE)

# libpluto doesn't use polylib
libpluto_la_LIBADD = \
   ../piplib/libpiplib64.la \
   ../openscop/libosl.la \
   ../clan/libclan.la \
   ../candl/libcandl.la \
   $(ISL_LIBADD) \
   -lm -lgomp
