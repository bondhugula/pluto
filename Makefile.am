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

SUBDIRS = piplib polylib $(MAY_ISL) openscop cloog-isl clan candl pet lib tool

ACLOCAL_AMFLAGS = -I m4

dist-hook:
	rm -rf `find $(distdir)/doc -name CVS`

bin_SCRIPTS = polycc plutune getversion.sh
CLEANFILES: $(bin_SCRIPTS) parsetab.py
EXTRA_DIST = polycc.sh.in  examples test

pkginclude_HEADERS = include/pluto/pluto.h include/pluto/matrix.h

polycc: polycc.sh
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
	cp -f tool/pluto bin/
	cp -r polycc bin/

pclean:
	$(MAKE) -C lib clean
	$(MAKE) -C tool clean

test_libpluto: test/test_libpluto.c
	libtool --mode=link --tag=CC $(CC) $(CFLAGS) $<  -Llib/.libs/ -I $(top_srcdir)/include $(ISL_INCLUDE) $(ISL_LIBADD) -I openscop/include -lpluto -lgomp -o test_libpluto

unit_tests: test/unit_tests.c
	libtool --mode=link --tag=CC $(CC) $(CFLAGS) $< -Llib/.libs/ -I $(top_srcdir)/include $(ISL_INCLUDE) $(ISL_LIBADD) -I$(top_srcdir)/lib -lpluto -o unit_tests

force:
	true

libpluto: force
	$(MAKE) $(MFLAGS) -C lib libpluto.la

pluto: libpluto
	$(MAKE) $(MFLAGS) -C tool pluto

test: PHONY

PHONY: pluto polycc test_libpluto unit_tests
	./test.sh

test-alt: polycc pluto
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
