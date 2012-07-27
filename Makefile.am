#
# Makefile.am
#
# Copyright (C) 2007 Uday Bondhugula 
# 
#

RM = rm -f
LN = ln -s

SUBDIRS = piplib isl cloog-isl clan-0.6.0 candl-0.4.0 src

ACLOCAL_AMFLAGS = -I m4

dist-hook:
	rm -rf `find $(distdir)/doc -name CVS`


bin_SCRIPTS = polycc vloog ploog plorc plutune getversion.sh
CLEANFILES: $(bin_SCRIPTS) parsetab.py
EXTRA_DIST = polycc.sh.in  examples test


pkginclude_HEADERS = include/pluto/libpluto.h

polycc: polycc.sh
	rm -f polycc
	echo "#! " $(BASH) > polycc
	cat $(srcdir)/polycc.sh >> polycc
	chmod ugo+x polycc

.PHONY: bin binit

install:
	@echo "No install target"
	@echo "Please run polycc in the top-level directory"

bindir:
	@if [ ! -d bin ];then \
		echo mkdir bin; \
		mkdir bin; \
		fi

binit:
	cp -f src/pluto bin/
	cp -r polycc bin/

pclean:
	$(MAKE) -C src clean

test: .PHO

.PHO: 
	./test.sh
	./test.sh --tile --parallel --unroll
	./test.sh --maxfuse
	./test.sh --nofuse
