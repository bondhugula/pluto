AC_DEFUN([CANDL_ARG_LIBS_DEPENDENCIES],
[
dnl Add $prefix to the library path (convenience).
  if test -e ${prefix}/include; then
    CPPFLAGS="${CPPFLAGS} -I${prefix}/include"
  fi;
  if test -e ${prefix}/lib; then
    LDFLAGS="${LDFLAGS} -L${prefix}/lib"
  fi;
dnl Offer --with-piplib.
  AC_ARG_WITH(piplib,
	      AC_HELP_STRING([--with-piplib=DIR],
              	             [DIR Location of PIPLib package]),
              [with_piplib=$withval;
	       CPPFLAGS="${CPPFLAGS} -I$withval/include";
	       LDFLAGS="${LDFLAGS} -L$withval/.libs"
	      ],
              [with_piplib=yes])
dnl Check for piplib existence.
dnl UB: Do not check ofr piplib existence
dnl It wouldn't have been built by now
dnl  AS_IF([test "x$with_piplib" != xno],
dnl	[AC_CHECK_LIB([piplib$BITS], [pip_solve],
dnl	 [LIBS="-lpiplib$BITS $LIBS";
dnl	 AC_DEFINE([HAVE_LIBPIPLIB], [1], [Define if you have libpiplib$BITS])
dnl       ],
dnl     [if test "x$with_piplib" != xcheck; then
dnl           AC_MSG_FAILURE([--with-piplib was given, but test for piplib failed])
dnl          fi
dnl         ])
dnl	])
LIBS="-lpiplib$BITS $LIBS"

dnl
dnl
dnl Offer --with-scoplib.
  AC_ARG_WITH(scoplib,
	      AC_HELP_STRING([--with-scoplib=DIR],
              	             [DIR Location of ScopLib package]),
              [with_scoplib=$withval;
	       CPPFLAGS="${CPPFLAGS} -I$withval/include";
	       LDFLAGS="${LDFLAGS} -L$withval/source/.libs"
	      ],
              [with_scoplib=check])
dnl Check for scoplib existence.
dnl UB: Do not check for its existence
dnl  AS_IF([test "x$with_scoplib" != xno],
dnl	[AC_CHECK_LIB([scoplib], [scoplib_scop_read],
dnl	 [LIBS="-lscoplib $LIBS";
dnl	 DEFINE_HAS_SCOPLIB_LIB="# define CANDL_SUPPORTS_SCOPLIB"
dnl	 ],
dnl       [DEFINE_HAS_SCOPLIB_LIB=""
dnl  	  if test "x$with_scoplib" != xcheck; then
dnl           AC_MSG_FAILURE([Test for ScopLib failed. Use --with-scoplib to specify libscoplib path.])
dnl          fi
dnl         ])
dnl	])
LIBS="-lscoplib $LIBS"
DEFINE_HAS_SCOPLIB_LIB="# define CANDL_SUPPORTS_SCOPLIB"

])

