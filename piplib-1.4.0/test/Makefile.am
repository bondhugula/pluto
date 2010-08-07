#
# Copyright 2002 (c) Leiden University, LIACS,
# All Rights Reserved.
#
# Permission is granted to copy, use, and distribute
# for any commercial or noncommercial purpose under the terms
# of the GNU Lesser General Public license, version 2.1, February 1999
# (see file : ../LICENSE).
#
# @author: Bart Kienhuis
# $Id: Makefile.am,v 1.1.1.1 2010/01/25 07:06:37 uday Exp $
#
# ndCed: useful line for names with their extensions:
# diff -w xyz `basename $$x .dat`.ll ;

# pbs with boulet bouleti and dirk 
#	boulet  \
#	bouleti \
#	dirk    \
#

PIPTEST =       \
	crescat \
	discr   \
	equus   \
	invert  \
	linear  \
	lineri  \
	loz     \
	max     \
	maxb    \
	pairi   \
	petit   \
	rairo   \
	rairoi  \
	test    \
	test2   \
	test2i  \
	test3   \
	test3i  \
	test4   \
	test4i  \
	test5   \
	test5i  \
	test6   \
	test6i  \
	test7   \
	test7i  \
	test8i  \
	test9i  \
	test10i \
	test11  \
	test11i \
	test12  \
	test12i

EXTRA_DIST = \
	$(PIPTEST:%=%.dat) \
	$(PIPTEST:%=%.ll)

generate:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                 Generate files                *"
	@echo "              *-----------------------------------------------*/"
	for x in $(PIPTEST) ; do \
		echo "Generate $(srcdir)/$$x.ll" ; \
		$(top_builddir)/pip64$(EXEEXT) $(srcdir)/$$x.dat $(srcdir)/$$x.ll ; \
	done

test:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                  Testing PIP                  *"
	@echo "              *-----------------------------------------------*/"
	@failedtest=0; \
	for x in $(PIPTEST) ; do \
		echo "Verify file $$x" ; \
		$(top_builddir)/$(PIP) -s $(srcdir)/$$x.dat > xyz; \
		diff -w xyz $(srcdir)/$$x.ll ; \
		result=$$?; \
		if [ "$$result" -ne "0" ]; then \
			echo "Error: $$x is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
		else \
			echo "$$x passed"; \
	        fi; \
	done ; \
	rm xyz ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
	else \
	        echo "PIP works correctly :-)"; \
        fi
