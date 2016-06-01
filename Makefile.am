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
ISL_INCLUDE += -I$(top_srcdir)/isl/include
ISL_LIBADD += $(top_srcdir)/isl/libisl.la
endif


SUBDIRS = piplib polylib $(MAY_ISL) openscop cloog-isl clan candl pet src

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

test_libpluto: src/test_libpluto.c
	libtool --mode=link gcc -g src/test_libpluto.c -Lsrc/.libs/ $(ISL_INCLUDE) $(ISL_LIBADD) -I include -I openscop/include -lpluto -lgomp -o test_libpluto

test: .PHO

.PHO: 
	./test.sh --silent
	./test.sh --silent --tile --parallel --cloogsh --nocloogbacktrack
	./test.sh --silent --lastwriter --tile --parallel --cloogsh --nocloogbacktrack
	./test.sh --silent --partlbtile --cloogsh --nocloogbacktrack
	./test.sh --silent --lastwriter
	./test.sh --silent --nofuse --tile --parallel --cloogsh --nocloogbacktrack
	./test.sh --silent --lastwriter --partlbtile --cloogsh --nocloogbacktrack
	./test.sh --silent --maxfuse --tile --parallel --cloogsh --nocloogbacktrack
	./test.sh --silent --nofuse
	./test.sh --silent --iss
