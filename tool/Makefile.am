# Pluto
#
# Copyright (C) Uday Bondhugula 2007
#
SUBDIRS =

MAINTAINERCLEANFILES = Makefile.in

if EXTERNAL_ISL
# ISL_INCLUDE, ISL_LIBADD were already defined
else
ISL_INCLUDE += -I$(top_srcdir)/isl/include -I$(top_builddir)/isl/include
ISL_LIBADD += $(top_builddir)/isl/libisl.la
endif

bin_PROGRAMS = pluto

pluto_SOURCES = \
		ast_transform.h \
		ast_transform.c \
		constraints_polylib.c \
		constraints_polylib.h \
		main.cpp \
		osl_pluto.c \
		pet_to_pluto.cpp \
		pluto_codegen_if.c \
		version.h

pluto_CXXFLAGS = $(OPT_FLAGS) $(DEBUG_FLAGS) \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP -DPLUTO_OPENCL \
   -I../include \
   -I$(top_srcdir)/include \
   -I$(top_srcdir)/piplib/include \
   -I../clan/include \
   -I$(top_srcdir)/clan/include \
   -I$(top_srcdir)/pet/include \
   $(ISL_INCLUDE) \
   -I../cloog-isl/include \
   -I$(top_srcdir)/cloog-isl/include \
   -I../openscop/include \
   -I$(top_srcdir)/openscop/include \
   -I../candl/include \
   -I$(top_srcdir)/candl/include \
   -I../polylib/include \
   -I$(top_srcdir)/polylib/include \
   -I$(top_srcdir)/lib

# To avoid same name object from being created with libtool and without it.
pluto_CFLAGS = $(AM_CFLAGS) $(pluto_CXXFLAGS)

pluto_LDFLAGS = -static
# libpet needs to come later after libisl so that the right isl is linked.
pluto_LDADD =  \
   ../cloog-isl/libcloog-isl.la \
   ../polylib/libpolylib64.la \
   ../clan/libclan.la \
   ../candl/libcandl.la \
   ../openscop/libosl.la \
   ../lib/libpluto.la \
   $(ISL_LIBADD) \
   ../piplib/libpiplib_dp.la \
   ../pet/libpet.la \
   -lm
