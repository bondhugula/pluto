#
# Makefile.am
#
# Copyright (C) 2007 Uday Bondhugula 
# 
#

RM = rm -f
LN = ln -s

if EXTERNAL_ISL
MAY_ISL =
# ISL_INCLUDE and ISL_LIBADD were already defined
else
MAY_ISL = isl
ISL_INCLUDE += -I$(top_srcdir)/isl/include -I$(top_builddir)/isl/include
ISL_LIBADD += $(top_builddir)/isl/libisl.la
endif


SUBDIRS = piplib polylib $(MAY_ISL) openscop cloog-isl clan candl pet src

ACLOCAL_AMFLAGS = -I m4

dist-hook:
	rm -rf `find $(distdir)/doc -name CVS`


bin_SCRIPTS = polycc vloog ploog plorc plutune getversion.sh
CLEANFILES: $(bin_SCRIPTS) parsetab.py
EXTRA_DIST = polycc.sh.in  examples test

pkginclude_HEADERS = include/pluto/pluto.h include/pluto/matrix.h

polycc: polycc.sh
	rm -f polycc
	echo "#! " $(BASH) > polycc
	cat polycc.sh >> polycc
	chmod ugo+x polycc

.PHONY: bin binit

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

test_libpluto: test/test_libpluto.c
	libtool --mode=link --tag=CC $(CC) $(CLFAGS) $<  -Lsrc/.libs/ -I $(top_srcdir)/include $(ISL_INCLUDE) $(ISL_LIBADD) -I openscop/include -lpluto -lgomp -o test_libpluto

unit_tests: test/unit_tests.c
	libtool --mode=link --tag=CC $(CC) $(CLFAGS) $< -Lsrc/.libs/ -I $(top_srcdir)/include -I $(top_srcdir)/src -lpluto -o unit_tests

src/pluto: force
	$(MAKE) $(MFLAGS) -C src pluto

force:
	true

test: PHONY
	
PHONY: polycc src/pluto test_libpluto unit_tests
	./test.sh

test-alt: polycc src/pluto
	./test-alt.sh --silent --notile --noparallel
	./test-alt.sh --silent --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --nodiamond-tile --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --lastwriter --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --lastwriter --nodiamond-tile --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --lastwriter
	./test-alt.sh --silent --nofuse --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --maxfuse --cloogsh --nocloogbacktrack
	./test-alt.sh --silent --nofuse
	./test-alt.sh --silent --iss
