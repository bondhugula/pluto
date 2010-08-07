#
#  /**-------------------------------------------------------------------**
#   **                            PIP/PipLib                             **
#   **-------------------------------------------------------------------**
#   **                           makefile.in                             **
#   **-------------------------------------------------------------------**
#   **                   First version: july 28th 2001                   **
#   **-------------------------------------------------------------------**/
#
# Makefile.am (or Makefile if generated) of PIP/PipLib. Makefile.am is not a
# makefile, you must run the 'configure' shellscript to generate the Makefile
# thanks to this file. 

SUBDIRS = . example doc test

ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = @pips@
EXTRA_PROGRAMS = pip32 pip64 pipMP
lib_LTLIBRARIES = @piplibs@
EXTRA_LTLIBRARIES = libpiplib32.la libpiplib64.la libpiplibMP.la

LIBRARY_FILES = \
	source/integrer.c \
	source/piplib.c \
	source/sol.c \
	source/tab.c \
	source/traiter.c \
	source/funcall.h \
	source/pip.h \
	source/sol.h \
	source/tab.h \
	source/type.h

pkginclude_HEADERS = \
	include/piplib/piplib.h \
	include/piplib/piplib32.h \
	include/piplib/piplib64.h \
	include/piplib/piplibMP.h

libpiplib32_la_SOURCES = $(LIBRARY_FILES)
libpiplib32_la_CPPFLAGS = @DFLAGS32@ -I$(srcdir)/include
libpiplib32_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpiplib64_la_SOURCES = $(LIBRARY_FILES)
libpiplib64_la_CPPFLAGS = @DFLAGS64@ -I$(srcdir)/include
libpiplib64_la_LDFLAGS = -rpath $(libdir) -version-info @versioninfo@
libpiplibMP_la_SOURCES = $(LIBRARY_FILES)
libpiplibMP_la_CPPFLAGS = @MP_DFLAGS@ -I$(srcdir)/include
libpiplibMP_la_LDFLAGS = @EXTRA_LIBS@ -rpath $(libdir) -version-info @versioninfo@

pip32_SOURCES = source/maind.c
pip32_CPPFLAGS = @DFLAGS32@ -I$(srcdir)/include
pip32_LDADD = libpiplib32.la
pip64_SOURCES = source/maind.c
pip64_CPPFLAGS = @DFLAGS64@ -I$(srcdir)/include
pip64_LDADD = libpiplib64.la
pipMP_SOURCES = source/maind.c
pipMP_CPPFLAGS = @MP_DFLAGS@ -I$(srcdir)/include
pipMP_LDADD = libpiplibMP.la @EXTRA_LIBS@

version.h: @GIT_INDEX@
	echo '#define PIPLIB_HEAD "'`$(top_builddir)/genversion.sh`'"' > $@

EXTRA_DIST = \
	autogen.sh

dist-hook:
	(cd doc; make piplib.pdf) && cp doc/piplib.pdf $(distdir)/doc/
	$(top_builddir)/genversion.sh > $(distdir)/PIPLIB_HEAD
	
check:
	@for pip in @pips@; do \
	    (cd test; $(MAKE) test PIP=$$pip); \
	done
	@for example in @examples@; do \
	    (cd example; $(MAKE) test  EXAMPLE=$$example); \
	done
	
total:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                   PIP/PipLib                  *"
	@echo "              *-----------------------------------------------*/"
	$(MAKE) uninstall
	$(MAKE) clean
	$(MAKE)
	$(MAKE) install
