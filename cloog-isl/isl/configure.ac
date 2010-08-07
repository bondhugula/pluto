AC_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE(isl, 0.00)

AC_PROG_CC

AC_PROG_LIBTOOL

AX_CREATE_STDINT_H(include/isl_stdint.h)

AC_ARG_WITH(gmp_prefix,
	[AS_HELP_STRING([--with-gmp-prefix=DIR],
			[Location of GMP installation])])

AC_SUBST(GMP_CPPFLAGS)
AC_SUBST(GMP_LDFLAGS)
if test "x$with_gmp_prefix" != "x"; then
	isl_configure_args="$isl_configure_args --with-gmp=$with_gmp_prefix"
	GMP_CPPFLAGS="-I$with_gmp_prefix/include"
	GMP_LDFLAGS="-L$with_gmp_prefix/lib"
fi
SAVE_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$GMP_CPPFLAGS $CPPFLAGS"
need_get_memory_functions=false
AC_CHECK_DECLS(mp_get_memory_functions,[],[
	need_get_memory_functions=true
],[#include <gmp.h>])
CPPFLAGS="$SAVE_CPPFLAGS"
AM_CONDITIONAL(NEED_GET_MEMORY_FUNCTIONS, test x$need_get_memory_functions = xtrue)

AC_DEFUN([ISL_SUBMODULE],[
	AC_ARG_WITH($1_prefix,
		[AS_HELP_STRING([--with-$1-prefix=DIR],
				[Location of system $1])])
	AC_ARG_WITH($1_builddir,
		    [AS_HELP_STRING([--with-$1-builddir=DIR],
				    [Location of $1 builddir])])
	if test "x$with_$1_prefix" != "x"; then
		if test "x$with_$1" != "x" -a "x$with_$1" != "xsystem"; then
			AC_MSG_ERROR([Setting $with_$1_prefix implies use of system $1])
		fi
		with_$1="system"
	fi
	if test "x$with_$1_builddir" != "x"; then
		if test "x$with_$1" != "x" -a "x$with_$1" != "xbuild"; then
			AC_MSG_ERROR([Setting $with_$1_builddir implies use of build $1])
		fi
		with_$1="build"
		$1_srcdir=`echo @abs_srcdir@ | $with_$1_builddir/config.status --file=-`
		AC_MSG_NOTICE($1 sources in $$1_srcdir)
	fi
	case "$with_$1" in
	no|bundled|system|build)
		;;
	*)
		if test -d $srcdir/.git -a \
			-d $srcdir/$1 -a \
			! -d $srcdir/$1/.git; then
			AC_MSG_WARN(
[git repo detected, but submodule $1 not initialized])
			AC_MSG_WARN([You may want to run])
			AC_MSG_WARN([	git submodule init])
			AC_MSG_WARN([	git submodule update])
			AC_MSG_WARN([	sh autogen.sh])
		fi
		if test -f $srcdir/$1/configure; then
			with_$1="bundled"
		fi
		;;
	esac
])

AC_ARG_WITH(polylib,
	[AS_HELP_STRING([--with-polylib=build|system|no],
			[Which PolyLib to use])])
ISL_SUBMODULE(polylib)
if test "x$with_polylib" = "x"; then
	with_polylib="no"
fi
AC_MSG_CHECKING([which polylib to use])
AC_MSG_RESULT($with_polylib)

have_polylib=false
AC_SUBST(POLYLIB_CPPFLAGS)
AC_SUBST(POLYLIB_LDFLAGS)
AC_SUBST(POLYLIB_LIBS)
case "$with_polylib" in
	build)
		polylibs=`echo @polylibs@ | $with_polylib_builddir/config.status --file=-`
		AC_MSG_NOTICE(Configured polylibs: $polylibs)
		isl_cv_polylib=missing
		for bits in $polylibs; do
			if test "$bits" = "libpolylibgmp.la"; then
				isl_cv_polylib=ok
			fi
		done
		if test "$isl_cv_polylib" = "missing"; then
			AC_MSG_ERROR(no gmp polylib configured)
		fi
		POLYLIB_CPPFLAGS="-I$with_polylib_builddir/include -I$polylib_srcdir/include"
		POLYLIB_LIBS="$with_polylib_builddir/libpolylibgmp.la"
	;;
	system)
		POLYLIB_LIBS="-lpolylibgmp"
		if test "x$with_polylib_prefix" != "x"; then
			POLYLIB_CPPFLAGS="-I$with_polylib_prefix/include"
			POLYLIB_LDFLAGS="-L$with_polylib_prefix/lib"
		fi
		SAVE_CPPFLAGS="$CPPFLAGS"
		SAVE_LDFLAGS="$LDFLAGS"
		CPPFLAGS="$POLYLIB_CPPFLAGS $CPPFLAGS"
		LDFLAGS="$POLYLIB_LDFLAGS $LDFLAGS"
		AC_CHECK_LIB(polylibgmp, PolyhedronTSort,[ true ],[
			AC_MSG_ERROR(Need polylib)
		])
		CPPFLAGS="$SAVE_CPPFLAGS"
		LDFLAGS="$SAVE_LDFLAGS"
	;;
	no)
	;;
	*)
		AC_MSG_ERROR(unsupported)
	;;
esac
if test "$with_polylib" != "no"; then
	AC_DEFINE(ISL_POLYLIB,,polylib is available)
	have_polylib=true
fi
AM_CONDITIONAL(HAVE_POLYLIB, test x$have_polylib = xtrue)

AC_ARG_WITH(piplib,
	[AS_HELP_STRING([--with-piplib=build|bundled|system|no],
			[Which piplib to use])])
ISL_SUBMODULE(piplib)
if test "x$with_piplib" = "x"; then
	with_piplib="system"
fi
AC_MSG_CHECKING([which piplib to use])
AC_MSG_RESULT($with_piplib)

have_piplib=false
AC_SUBST(PIPLIB_CPPFLAGS)
AC_SUBST(PIPLIB_LDFLAGS)
AC_SUBST(PIPLIB_LIBS)
case "$with_piplib" in
	bundled)
		PIPLIB_CPPFLAGS="-I$srcdir/piplib/include"
		isl_configure_args="$isl_configure_args --with-bits=gmp"
	;;
	build)
		PIPLIB_CPPFLAGS="-I$piplib_srcdir/include"
		PIPLIB_LIBS="$with_piplib_builddir/libpiplibMP.la"
	;;
	system)
		PIPLIB_LIBS="-lpiplibMP"
		if test "x$with_piplib_prefix" != "x"; then
			PIPLIB_CPPFLAGS="-I$with_piplib_prefix/include"
			PIPLIB_LDFLAGS="-L$with_piplib_prefix/lib"
		fi
		SAVE_CPPFLAGS="$CPPFLAGS"
		SAVE_LDFLAGS="$LDFLAGS"
		CPPFLAGS="$PIPLIB_CPPFLAGS $CPPFLAGS"
		LDFLAGS="$PIPLIB_LDFLAGS $LDFLAGS"
		AC_CHECK_LIB(piplibMP, pip_solve,[
			AC_CHECK_MEMBER(PipOptions.Urs_parms, [], [
				AC_MSG_ERROR([Piplib too old; please install version 1.3.6 or newer])
			],[#include <piplib/piplibMP.h>])
		],[
			AC_MSG_ERROR([Piplib not found])
		])
		CPPFLAGS="$SAVE_CPPFLAGS"
		LDFLAGS="$SAVE_LDFLAGS"
	;;
	no)
	;;
	*)
		AC_MSG_ERROR(unsupported)
	;;
esac
if test "$with_piplib" != "no"; then
	AC_DEFINE(ISL_PIPLIB,,piplib is available)
	have_piplib=true
fi
AM_CONDITIONAL(HAVE_PIPLIB, test x$have_piplib = xtrue)
AM_CONDITIONAL(BUNDLED_PIPLIB, test $with_piplib = bundled)

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_HEADERS(include/isl_ctx.h)
AC_CONFIG_FILES(Makefile)
if test $with_piplib = bundled; then
	AC_CONFIG_SUBDIRS(piplib)
fi
AC_CONFIG_COMMANDS_POST([
	dnl pass on arguments to subdir configures, but don't
	ac_configure_args="$ac_configure_args $isl_configure_args"
	ac_configure_args="$ac_configure_args $isl_configure_args"
])
AC_OUTPUT
