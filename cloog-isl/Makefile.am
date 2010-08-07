#
#   /**-------------------------------------------------------------------**
#    **                              CLooG                                **
#    **-------------------------------------------------------------------**
#    **                           makefile.in                             **
#    **-------------------------------------------------------------------**
#    **                   First version: october 25th 2001                **
#    **-------------------------------------------------------------------**/
#
# makefile.in (or makefile if generated) of CLooG, the Chunky LOOp Generator.
# makefile.in is not a makefile, you must run the 'configure' shellscript to
# generate the makefile thanks to this file.

#/*****************************************************************************
# *               CLooG : the Chunky Loop Generator (experimental)            *
# *****************************************************************************
# *                                                                           *
# * Copyright (C) 2001 Cedric Bastoul                                         *
# *                                                                           *
# * This library is free software; you can redistribute it and/or             *
# * modify it under the terms of the GNU Lesser General Public                *
# * License as published by the Free Software Foundation; either              *
# * version 2.1 of the License, or (at your option) any later version.        *
# *                                                                           *
# * This library is distributed in the hope that it will be useful,           *
# * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
# * Lesser General Public License for more details.                           *
# *                                                                           *
# * You should have received a copy of the GNU Lesser General Public          *
# * License along with this library; if not, write to the Free Software       *
# * Foundation, Inc., 51 Franklin Street, Fifth Floor,                        *
# * Boston, MA  02110-1301  USA                                               *
# *                                                                           *
# * CLooG, the Chunky Loop Generator                                          *
# * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                        *
# *                                                                           *
# *****************************************************************************/
if BUNDLED_ISL
    MAYBE_ISL = isl
    ISL_LA = $(top_builddir)/isl/libisl.la
endif

SUBDIRS = $(MAYBE_ISL) . doc test
DIST_SUBDIRS = $(MAYBE_ISL) . doc test

ACLOCAL_AMFLAGS = -I m4

FORCE:
isl/libisl.la: FORCE
	cd isl; $(MAKE) $(AM_MAKEFLAGS) libisl.la

bin_PROGRAMS = cloog
lib_LTLIBRARIES = @CLOOG_LIBRARIES@
EXTRA_LTLIBRARIES = libcloog-polylib.la libcloog-isl.la

if NEED_GET_MEMORY_FUNCTIONS
GET_MEMORY_FUNCTIONS=source/mp_get_memory_functions.c
endif

SOURCES_CORE = \
	$(GET_MEMORY_FUNCTIONS) \
	source/block.c \
	source/clast.c \
	source/state.c \
	source/int.c \
	source/loop.c \
	source/names.c \
	source/options.c \
	source/pprint.c \
	source/program.c \
	source/statement.c \
	source/version.c

# avoid -I$(top_builddir)/include/cloog ending up in DEFAULT_INCLUDES
# as that may result in polylib/matrix.h confusion
DEFAULT_INCLUDES = -I.
INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include
AM_CFLAGS = -Wall
libcloog_polylib_la_CPPFLAGS = @POLYLIB_CPPFLAGS@ -DCLOOG_POLYLIB
libcloog_polylib_la_LDFLAGS = -rpath $(libdir) @POLYLIB_LDFLAGS@
libcloog_polylib_la_LIBADD = @POLYLIB_LIBS@
libcloog_polylib_la_SOURCES = \
	$(SOURCES_CORE) \
	source/polylib/domain.c \
	source/polylib/matrix.c \
	source/polylib/backend.c
libcloog_isl_la_CPPFLAGS = @ISL_CPPFLAGS@ -DCLOOG_ISL
libcloog_isl_la_LDFLAGS = -rpath $(libdir) @ISL_LDFLAGS@
libcloog_isl_la_LIBADD = @ISL_LIBS@ $(ISL_LA)
libcloog_isl_la_SOURCES = \
	$(SOURCES_CORE) \
	source/isl/domain.c \
	source/isl/constraints.c \
	source/isl/backend.c
LDADD = libcloog-@BACKEND@.la
cloog_DEPENDENCIES = libcloog-@BACKEND@.la
cloog_SOURCES = source/cloog.c

pkginclude_HEADERS = \
	include/cloog/block.h \
	include/cloog/clast.h \
	include/cloog/cloog.h \
	include/cloog/state.h \
	include/cloog/domain.h \
	include/cloog/loop.h \
	include/cloog/constraints.h \
	include/cloog/names.h \
	include/cloog/options.h \
	include/cloog/pprint.h \
	include/cloog/program.h \
	include/cloog/statement.h \
	include/cloog/version.h
nodist_pkginclude_HEADERS = \
	include/cloog/int.h

pkgpolylibincludedir = $(pkgincludedir)/polylib
pkgpolylibinclude_HEADERS = \
	include/cloog/polylib/cloog.h \
	include/cloog/polylib/domain.h \
	include/cloog/polylib/matrix.h
nodist_pkgpolylibinclude_HEADERS = \
	include/cloog/polylib/backend.h

pkgislincludedir = $(pkgincludedir)/isl
pkgislinclude_HEADERS = \
	include/cloog/isl/backend.h \
	include/cloog/isl/cloog.h \
	include/cloog/isl/domain.h \
	include/cloog/isl/matrix.h

version.h: @GIT_INDEX@
	echo '#define CLOOG_HEAD "'`$(top_builddir)/genversion.sh`'"' > $@

EXTRA_DIST = autoconf/Doxyfile.in doc/images examples

dist-hook:
	rm -f $(distdir)/test/Makefile
	$(top_builddir)/genversion.sh > $(distdir)/CLOOG_HEAD
	(cd doc; make cloog.pdf) && cp doc/cloog.pdf $(distdir)/doc/

#/*****************************************************************************
# *                                   Rules                                   *
# *****************************************************************************/

check:
	$(MAKE) test -C test

valcheck:
	$(MAKE) valgrind -C test

total:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                     CLooG                     *"
	@echo "              *-----------------------------------------------*/"
	$(MAKE) uninstall
	$(MAKE) clean
	$(MAKE)
	$(MAKE) install

doc:
	@echo "             /*-----------------------------------------------*"
	@echo "              *        Generating CLooG's documentation       *"
	@echo "              *-----------------------------------------------*/"
	doxygen ./autoconf/Doxyfile
