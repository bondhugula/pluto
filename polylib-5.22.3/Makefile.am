SUBDIRS = include 

CPPFLAGS = -I$(srcdir)/include -Iinclude @CPPFLAGS@

ACLOCAL_AMFLAGS = -I m4

# poly_src = $(srcdir)/applications
poly_arith = $(srcdir)/source/arith

CFILES= \
	source/arith/errors.c \
	source/kernel/errormsg.c \
	source/kernel/vector.c \
	source/kernel/matrix.c \
	source/kernel/polyhedron.c \
	source/kernel/polyparam.c \
	source/kernel/param.c \
	source/kernel/alpha.c \
	source/kernel/matrix_addon.c \
	source/kernel/matrix_permutations.c \
	source/kernel/compress_parms.c \
	source/kernel/SolveDio.c \
	source/kernel/Lattice.c \
	source/kernel/Matop.c \
	source/kernel/NormalForms.c

sources = $(CFILES) \
	$(poly_arith)/assert.h \
	$(poly_arith)/arithmetique.h \
	$(poly_arith)/arithmetic_errors.h

EXTRA_DIST = doc source/oldpolytest.c

lib_LTLIBRARIES = @polylibs@
libpolylib32_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpolylib32_la_SOURCES = $(sources)
libpolylib32_la_CPPFLAGS = -DPOLYLIB_BITS=32
libpolylib32_la_LIBADD = @LTLIBOBJS@
libpolylib64_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpolylib64_la_SOURCES = $(sources)
libpolylib64_la_CPPFLAGS = -DPOLYLIB_BITS=64
libpolylib64_la_LIBADD = @LTLIBOBJS@
libpolylib128_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpolylib128_la_SOURCES = $(sources)
libpolylib128_la_CPPFLAGS = -DPOLYLIB_BITS=128
libpolylib128_la_LIBADD = @LTLIBOBJS@
libpolylibgmp_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpolylibgmp_la_SOURCES = $(sources)
libpolylibgmp_la_CPPFLAGS = -DGNUMP
libpolylibgmp_la_LIBADD = @LTLIBOBJS@


EXTRA_LTLIBRARIES = libpolylib32.la libpolylib64.la libpolylib128.la \
		    libpolylibgmp.la

multi_noinst_programs = 

noinst_PROGRAMS = \
	@multi_noinst_programs@

multi_bin_programs = 

bin_PROGRAMS = \
	@multi_bin_programs@

EXTRA_PROGRAMS = 

LDADD = @polylib@

check: tests

# UB: Removed
#tests: all
#(cd Test; $(MAKE) tests)
	
# UB: Removed
#longtest: all
#(cd Test; $(MAKE) long_tests)

# Make the 'javadoc' style documentation
document:
	-mkdir doc
	$(DOXYGEN) polylib.doxygen

dist-hook: document
	rm -rf `find $(distdir)/doc -name CVS`
	cp -r doc/codeDoc $(distdir)/doc

install-data-local:
	@test -z "$(pkgconfig_libdir)" || $(mkdir_p) "$(DESTDIR)$(pkgconfig_libdir)"
	$(INSTALL_DATA) $(pkgconfig_libfile) "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"

uninstall-local:
	rm -f "$(DESTDIR)$(pkgconfig_libdir)/$(pkgconfig_libfile)"
