#
#   /**-------------------------------------------------------------------**
#    **                              CLooG                                **
#    **-------------------------------------------------------------------**
#    **                            makefile                               **
#    **-------------------------------------------------------------------**
#    **                   First version: march 6th 2003                   **
#    **-------------------------------------------------------------------**/
#
# makefile for checking.

#/*****************************************************************************
# *               CLooG : the Chunky Loop Generator (experimental)            *
# *****************************************************************************
# *                                                                           *
# * Copyright (C) 2003 Cedric Bastoul                                         *
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
# * This file has been written with the help of the Bart Kienhuis's checking  *
# * makefile for PipLib.                                                      *
# *                                                                           *
# *****************************************************************************/

noinst_PROGRAMS = generate_test
INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include
LDADD = ../libcloog-@BACKEND@.la
generate_test_SOURCES = generate_test.c

FORCE:
../libcloog-@BACKEND@.la: FORCE
	cd ..; $(MAKE) $(AM_MAKEFLAGS) libcloog-@BACKEND@.la

ROOT = ..

FINITE_CLOOGTEST_C =                   \
	0D-1                           \
	0D-2                           \
	0D-3                           \
	1point-1                       \
	1point-2                       \
	4-param                        \
	basic-bounds-1                 \
	basic-bounds-2                 \
	basic-bounds-3                 \
	basic-bounds-4                 \
	basic-bounds-5                 \
	basic-bounds-6                 \
	block                          \
	byu98-1-2-3                    \
	cholesky                       \
	cholesky2                      \
	classen                        \
	constbound                     \
	darte                          \
	donotsimp                      \
	dot                            \
	dot2                           \
	double                         \
	durbin_e_s                     \
	emploi                         \
	esced                          \
	ex1                            \
	forwardsub-1-1-2               \
	forwardsub-2-1-2-3             \
	forwardsub-3-1-2               \
	gauss                          \
	gesced                         \
	gesced2                        \
	gesced3                        \
	guide                          \
	iftest                         \
	iftest2                        \
	largeur                        \
	levenshtein-1-2-3              \
	lineality-1-2                  \
	lineality-2-1-2                \
	logo                           \
	logopar                        \
	lu                             \
	lu2                            \
	lub                            \
	lux                            \
	min-1-1                        \
	min-2-1                        \
	min-3-1                        \
	min-4-1                        \
	mode                           \
	multi-mm-1                     \
	no_lindep                      \
	orc                            \
	rectangle                      \
	square+triangle-1-1-2-3        \
	swim                           \
	test                           \
	thomasset                      \
	tiling                         \
	uday_scalars		       \
	vivien                         \
	vivien2                        \
	wavefront                      \
	yosr                           \
	yosr2                          \
	./non_optimal/nul_complex1     \
	./non_optimal/usvd_e_t         \
	./non_optimal/youcef           \
	./reservoir/cholesky2          \
	./reservoir/fusion1            \
	./reservoir/fusion2            \
	./reservoir/jacobi2            \
	./reservoir/jacobi3            \
	./reservoir/lim-lam1           \
	./reservoir/lim-lam2           \
	./reservoir/lim-lam3           \
	./reservoir/lim-lam4           \
	./reservoir/lim-lam5           \
	./reservoir/lim-lam6           \
	./reservoir/liu-zhuge1         \
	./reservoir/loechner3          \
	./reservoir/loechner4          \
	./reservoir/loechner5          \
	./reservoir/mg-interp2         \
	./reservoir/mg-interp          \
	./reservoir/mg-psinv           \
	./reservoir/mg-resid           \
	./reservoir/mg-rprj3           \
	./reservoir/pingali1           \
	./reservoir/pingali2           \
	./reservoir/pingali3           \
	./reservoir/pingali4           \
	./reservoir/pingali5           \
	./reservoir/pingali6           \
	./reservoir/QR                 \
	./reservoir/tang-xue1          \
	./reservoir/two   

CLOOGTEST_C =                          \
	infinite                       \
	infinite2                      \
	infinite3                      \
	infinite4                      \
	$(FINITE_CLOOGTEST_C)

#	./reservoir/long               \
#	faber

CLOOGTEST_FORTRAN =                    \
	dartef                         \
	iftestf                        \
	levenshtein-1-2-3f             \
	yosrf

CLOOGTEST_STRIDED =                    \
	basic-bounds-2                 \
	nul_basic1                     \
	nul_basic2                     \
	nul_lcpc		       \
	multi-stride		       \
	multi-stride2		       \
	reservoir/bastoul3	       \
	reservoir/stride	       \
	reservoir/stride2

SPECIAL_TESTS = \
	vasilache \
	merge \
	equality \
	equality2 \
	stride \
	stride2
SPECIAL_OPTIONS = \
	"vasilache" "-f 8 -l 9" "\n" \
	"merge" "-f -1" "\n" \
	"equality" "-f -1 -l 2 -override" "\n" \
	"equality2" "-f -1 -l 4 -esp 1 -override" "\n" \
	"stride" "-f -1 -strides 1" "\n" \
	"stride2" "-f -1 -strides 1"

EXTRA_DIST = \
	$(CLOOGTEST_C:%=%.cloog) \
	$(CLOOGTEST_C:%=isl/%.c) \
	$(CLOOGTEST_C:%=polylib/%.c) \
	$(FINITE_CLOOGTEST_C:%=%.good.c) \
	$(CLOOGTEST_FORTRAN:%=%.cloog) \
	$(CLOOGTEST_FORTRAN:%=isl/%.f) \
	$(CLOOGTEST_FORTRAN:%=polylib/%.f) \
	$(CLOOGTEST_STRIDED:%=%.cloog) \
	$(CLOOGTEST_STRIDED:%=isl/%.c) \
	$(CLOOGTEST_STRIDED:%=polylib/%.c) \
	$(SPECIAL_TESTS:%=%.cloog) \
	$(SPECIAL_TESTS:%=isl/%.c) \
	$(SPECIAL_TESTS:%=polylib/%.c) \
	$(SPECIAL_TESTS:%=%.good.c)

generate:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                 Generate files                *"
	@echo "              *-----------------------------------------------*/"
	@for x in $(CLOOGTEST_C) ; do \
		echo "Generate $$x.c" ; \
		$(top_builddir)/cloog $(srcdir)/$$x.cloog -o $(srcdir)/@BACKEND@/$$x.c ; \
	done ; \
	for x in $(CLOOGTEST_FORTRAN) ; do \
		echo "Generate $$x.f" ; \
		$(top_builddir)/cloog $(srcdir)/$$x.cloog -o $(srcdir)/@BACKEND@/$$x.f ; \
	done ; \
	for x in $(CLOOGTEST_STRIDED) ; do \
		echo "Generate $$x.c" ; \
		$(top_builddir)/cloog $(srcdir)/$$x.cloog \
			-strides 1 -o $(srcdir)/@BACKEND@/$$x.c ; \
	done ; \
	(echo -e $(SPECIAL_OPTIONS)) | while read x options; do \
		echo "Generate file $$x ($$options)" ; \
		$(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog \
		    $$options > $(srcdir)/@BACKEND@/$$x.c ; \
	done

generate_good:
	@for x in $(FINITE_CLOOGTEST_C) ; do \
		echo "Generate $$x.good.c" ; \
		$(top_builddir)/cloog$(EXEEXT) -callable 1 $(srcdir)/$$x.cloog -o $(srcdir)/$$x.good.c ; \
	done; \
	(echo -e $(SPECIAL_OPTIONS)) | while read x options; do \
		echo "Generate $$x.good.c ($$options)" ; \
		$(top_builddir)/cloog$(EXEEXT) -callable 1 $$options $(srcdir)/$$x.cloog -o $(srcdir)/$$x.good.c ; \
	done

valgrind:
	@echo "             /*-----------------------------------------------*"
	@echo "              *            Checking files by Valgrind         *"
	@echo "              *-----------------------------------------------*/"
	for x in $(CLOOGTEST_C) ; do \
		echo "Valgrinding $$x.c" ; \
		valgrind $(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog ; \
	done ; \
	for x in $(CLOOGTEST_FORTRAN) ; do \
		echo "Valgrinding $$x.f" ; \
		valgrind $(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog ; \
	done

test_run: generate_test$(EXEEXT)
	@echo "             /*-----------------------------------------------*"
	@echo "              *                 Testing CLooG                 *"
	@echo "              *-----------------------------------------------*/"
	@failedtest=0 ; \
	echo "--- Testing C files... (correct output)" ; \
	for x in $(FINITE_CLOOGTEST_C) ; do \
		echo -n "Verify file $$x... " ; \
		echo -n "generating... " ; \
		$(top_builddir)/cloog$(EXEEXT) -q -callable 1 $(srcdir)/$$x.cloog > test_test.c ; \
		./generate_test$(EXEEXT) < $(srcdir)/$$x.cloog > test_main.c; \
		echo -n "compiling... " ; \
		$(COMPILE) -c test_test.c; \
		$(COMPILE) -Dtest=good -c $(srcdir)/$$x.good.c -o test_good.o; \
		$(LINK) test_main.c test_test.o test_good.o > /dev/null; \
		echo -n "comparing... " ; \
		./test_run$(EXEEXT); \
		result=$$?; \
		if test "$$result" -ne "0" ; then \
			failedtest=`expr $$failedtest + 1` ; \
		fi; \
		$(RM) $@; \
		echo; \
	done; \
	(echo -e $(SPECIAL_OPTIONS)) | while read x options; do \
		echo -n "Verify file $$x ($$options)... " ; \
		echo -n "generating... " ; \
		$(top_builddir)/cloog$(EXEEXT) -q -callable 1 $$options $(srcdir)/$$x.cloog > test_test.c ; \
		./generate_test$(EXEEXT) < $(srcdir)/$$x.cloog > test_main.c; \
		echo -n "compiling... " ; \
		$(COMPILE) -c test_test.c; \
		$(COMPILE) -Dtest=good -c $(srcdir)/$$x.good.c -o test_good.o; \
		$(LINK) test_main.c test_test.o test_good.o > /dev/null; \
		echo -n "comparing... " ; \
		./test_run$(EXEEXT); \
		result=$$?; \
		if test "$$result" -ne "0" ; then \
			failedtest=`expr $$failedtest + 1` ; \
		fi; \
		$(RM) $@; \
		echo; \
	done; \
	(echo -e $(SPECIAL_OPTIONS)) | while read x options; do \
		echo -n "Verify file $$x ($$options -sh 1)... " ; \
		echo -n "generating... " ; \
		$(top_builddir)/cloog$(EXEEXT) -q -callable 1 $$options -sh 1 $(srcdir)/$$x.cloog > test_test.c ; \
		./generate_test$(EXEEXT) < $(srcdir)/$$x.cloog > test_main.c; \
		echo -n "compiling... " ; \
		$(COMPILE) -c test_test.c; \
		$(COMPILE) -Dtest=good -c $(srcdir)/$$x.good.c -o test_good.o; \
		$(LINK) test_main.c test_test.o test_good.o > /dev/null; \
		echo -n "comparing... " ; \
		./test_run$(EXEEXT); \
		result=$$?; \
		if test "$$result" -ne "0" ; then \
			failedtest=`expr $$failedtest + 1` ; \
		fi; \
		$(RM) $@; \
		echo; \
	done; \
	if [ $$failedtest != 0 ] ; then \
		echo "$$failedtest tests failed" ; \
	else \
	        echo "CLooG works correctly :-) !" ; \
        fi
test:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                 Testing CLooG                 *"
	@echo "              *-----------------------------------------------*/"
	@failedtest=0 ; \
	echo "--- Testing C files... (expected output)" ; \
	for x in $(CLOOGTEST_C) ; do \
		echo "Verify file $$x" ; \
		$(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog > ttemp ; \
		diff -u -w --ignore-matching-lines='CLooG' ttemp $(srcdir)/@BACKEND@/$$x.c ; \
		result=$$?; \
		if [ "$$result" -ne "0" ] ; then \
			echo "!!! Error: $$x is not the same" ; \
			failedtest=`expr $$failedtest + 1` ; \
		else \
			echo "$$x passed" ; \
	        fi ; \
	done ; \
	echo "--- Testing FORTRAN files..." ; \
	for x in $(CLOOGTEST_FORTRAN) ; do \
		echo "Verify file $$x" ; \
		$(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog > ttemp ; \
		diff -u -w --ignore-matching-lines='CLooG' ttemp $(srcdir)/@BACKEND@/$$x.f ; \
		result=$$?; \
		if [ "$$result" -ne "0" ] ; then \
			echo "!!! Error: $$x is not the same" ; \
			failedtest=`expr $$failedtest + 1` ; \
		else \
			echo "$$x passed" ; \
	        fi ; \
	done ; \
	echo "--- Testing STRIDED files..." ; \
	for x in $(CLOOGTEST_STRIDED) ; do \
		echo "Verify file $$x" ; \
		$(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog \
		-strides 1 > ttemp ; \
		diff -u -w --ignore-matching-lines='CLooG' ttemp $(srcdir)/@BACKEND@/$$x.c ; \
		result=$$?; \
		if [ "$$result" -ne "0" ] ; then \
			echo "!!! Error: $$x is not the same" ; \
			failedtest=`expr $$failedtest + 1` ; \
		else \
			echo "$$x passed" ; \
	        fi ; \
	done ; \
	echo "--- Testing SPECIAL OPTIONS files..." ; \
	(echo -e $(SPECIAL_OPTIONS)) | (while read x options; do \
	        echo "Verify file $$x ($$options)" ; \
		$(top_builddir)/cloog$(EXEEXT) $(srcdir)/$$x.cloog \
			$$options > ttemp ; \
		diff -u -w --ignore-matching-lines='CLooG' ttemp $(srcdir)/@BACKEND@/$$x.c ; \
		result=$$?; \
		if [ "$$result" -ne "0" ] ; then \
			echo "!!! Error: $$x is not the same" ; \
			failedtest=`expr $$failedtest + 1` ; \
		else \
			echo "$$x passed" ; \
	        fi ; \
	done; exit $$failedtest; ); \
	failedtest=$$?; \
	rm ttemp ; \
	if [ $$failedtest != 0 ] ; then \
		echo "$$failedtest tests failed" ; \
	else \
	        echo "CLooG works correctly :-) !" ; \
        fi
