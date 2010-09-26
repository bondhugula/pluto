# PLuTo
#
# Copyright (C) Uday Kumar Bondhugula 2007
#

SUBDIRS = 

MAINTAINERCLEANFILES = Makefile.in

bin_PROGRAMS = pluto
pluto_SOURCES = constraints.c math_support.c ddg.c pluto.c post_transform.c framework.c main.c constraints.h math_support.h pluto.h program.c program.h tile.c
pluto_CPPFLAGS = -I../piplib-1.4.0/include \
   -I../clan-0.6.0/include \
   -I../clan-0.6.0/scoplib/scoplib-0.2.0/include \
   -I../candl-0.4.0/include \
   -I../isl/include -I../cloog-isl/include
pluto_LDADD = ../candl-0.4.0/source/libcandl.la \
   ../clan-0.6.0/source/libclan.la \
   ../piplib-1.4.0/libpiplib64.la \
   ../cloog-isl/libcloog-isl.la ../isl/libisl.la -lm
INCLUDES = 

