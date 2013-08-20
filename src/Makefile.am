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

pluto_SOURCES = constraints.c constraints_polylib.c constraints_isl.c math_support.c ddg.c transforms.c pluto.c distmem.c post_transform.c polyloop.c halcyon.c framework.c pluto_codegen_if.c main.c constraints.h math_support.h pluto.h program.c program.h tile.c unit_tests.c ast_transform.h ast_transform.c opencl.c
pluto_CPPFLAGS = -O3 -funroll-loops -fopenmp \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../include \
   -I../piplib/include \
   -I../clan-0.6.0/include \
   -I../clan-0.6.0/scoplib/scoplib-0.2.0/include \
   -I../candl-0.4.0/include \
   $(ISL_INCLUDE) \
   -I../cloog-isl/include \
   -I../polylib/include

pluto_LDFLAGS = -static -fopenmp
pluto_LDADD = ../candl-0.4.0/source/libcandl.la \
   ../clan-0.6.0/source/libclan.la \
   ../piplib/libpiplib64.la \
   ../cloog-isl/libcloog-isl.la \
   $(ISL_LIBADD) \
   ../polylib/libpolylib64.la \
   -lm

lib_LTLIBRARIES = libpluto.la

# libpluto doesn't use polylib-based functions
libpluto_la_SOURCES = constraints.c constraints_isl.c math_support.c ddg.c transforms.c pluto.c post_transform.c framework.c polyloop.c libpluto.c constraints.h math_support.h pluto.h program.c program.h tile.c
libpluto_la_CPPFLAGS = -O3 -funroll-loops -fopenmp \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../include \
   -I../piplib/include \
   -I../clan-0.6.0/include \
   -I../clan-0.6.0/scoplib/scoplib-0.2.0/include \
   -I../candl-0.4.0/include \
   $(ISL_INCLUDE)

# libpluto doesn't use polylib
libpluto_la_LIBADD = ../candl-0.4.0/source/libcandl.la \
   ../clan-0.6.0/source/libclan.la \
   ../piplib/libpiplib64.la \
   $(ISL_LIBADD) \
   -lm -lgomp
