dnl
dnl configure.in
dnl 
dnl PLUTO - parallelism and locality transformer (experimental)
dnl 
dnl Copyright (C) 2007-2008 Uday Bondhugula
dnl
dnl This program is free software; you can redistribute it and/or
dnl modify it under the terms of the GNU General Public License
dnl as published by the Free Software Foundation; either version 2
dnl of the License, or (at your option) any later version.
dnl
dnl The complete GNU General Public Licence Notice can be found as the
dnl `COPYING' file in the top-level directory.
dnl

AC_PREREQ(2.13)

AC_INIT(pluto, 0.12.0, [udayb@iisc.ac.in])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE()
AM_SILENT_RULES([yes])
AM_CONFIG_HEADER(config.h)

dnl Determine a C compiler to use
AC_PROG_CC
dnl Determine a C++ compiler to use
AC_PROG_CXX

AX_CXX_COMPILE_STDCXX(11, noext, mandatory)

AC_PROG_LIBTOOL
AC_PROG_SED

dnl clan, candl, cloog are all configured with 64-bit data type support
dnl NOTE: configuring these 64-bit data type support has nothing to do with arch
dnl being 32-bit/64-bit
CXXFLAGS="$CXXFLAGS -Wall -std=c++11"
if test CXX = icc; then
    LDFLAGS="$LDFLAGS -lirc -limf"
fi

debug="false"
dnl debugging flag
AC_ARG_ENABLE(debug,
              [  --enable-debug  enable compilation with debugging info ],
              debug="true")

AM_CONDITIONAL(PLUTO_DEBUG, test x$debug = xtrue)

if test $debug = "true"; then
    CXXFLAGS="$CXXFLAGS -g"
fi

AC_CHECK_PROG(LD,ld,ld)
AC_PROG_INSTALL
AC_PROG_LN_S


dnl check for some compiler characteristics
dnl MC_CHECK_GUIDING_DECLS
AC_SUBST(GUIDING_DECLS)

dnl check for typedefs, structures, compiler features
AC_C_CONST

AC_CHECK_HEADERS(math.h stdlib.h stdio.h assert.h string.h)

dnl set SOURCE_DIR to current directory
SOURCE_DIR=`pwd`
AC_SUBST(SOURCE_DIR)

AC_CHECK_PROG(CP,cp,cp)
AC_CHECK_PROG(RM,rm,rm -f)
AC_CHECK_PROG(MV,mv,mv)
AC_CHECK_PROG(CAT,cat,cat)
AC_CHECK_PROG(AR,ar,ar)
AC_CHECK_PROG(STRIP,strip,strip)

AC_CHECK_PROG(FILECHECK,FileCheck,yes,no,path = '/usr/lib64/llvm/')
AS_IF([test x"$FILECHECK" != x"yes"], [AC_MSG_ERROR([Please install LLVM FileCheck before configuring.])])

glpk="false"
dnl glpk flag
AC_ARG_ENABLE(glpk,
              [  --enable-glpk  enable ILP solving with GLPK], glpk="true")

gurobi="false"
dnl check for gurobi
AC_ARG_ENABLE(gurobi,
              [  --enable-gurobi  enable ILP solving with GUROBI], gurobi="true")

if test $glpk = "true"; then
    dnl check for glpk
    AC_CHECK_HEADER([glpk.h],[AC_DEFINE([GLPK],[1],[GLPK solver for solving ILP])]
                                CFLAGS="$CFLAGS -DGLPK" LDFLAGS="$LDFLAGS -lglpk"
                                , [AC_ERROR([glpk.h header not found])])
    AC_CHECK_LIB([glpk], [glp_create_prob], []
                 , [AC_ERROR([GLPK library not found])])
fi

dnl Offer --with-gurobi-prefix.
AC_ARG_WITH(gurobi-prefix,
            AC_HELP_STRING([--with-gurobi-prefix=DIR],
                           [DIR location of GUROBI package]),
                           [CPPFLAGS="${CPPFLAGS} -I$withval/include";
                            LDFLAGS="${LDFLAGS} -L$withval/lib";
                            # MAJOR_VERSION=`cat $withval/include/gurobi_c.h| grep VERSION_MAJOR`;
                            # echo "major version = ${MAJOR_VERSION}";
                            # LIBS="-lgurobi75 $LIBS";
                            GUROBI_PREFIX="$withval";
                            ])

if test $gurobi = "true"; then
    dnl check for gurobi
    AC_CHECK_HEADER([gurobi_c.h],[AC_DEFINE([GUROBI],[1],[GUROBI solver for solving ILP])]
                                   major_version=`cat $GUROBI_PREFIX/include/gurobi_c.h| grep VERSION_MAJOR| awk '{print $3}'`
                                   minor_version=`cat $GUROBI_PREFIX/include/gurobi_c.h| grep VERSION_MINOR| awk '{print $3}'`
                                   grb_version=gurobi$major_version$minor_version
                                   CFLAGS="$CFLAGS -DGUROBI" LDFLAGS="$LDFLAGS -l$grb_version"
                                   , [AC_ERROR([gurobi_c.h header not found])])
    AC_CHECK_LIB([$grb_version], [GRBnewmodel], []
                 , [AC_ERROR([Gurobi library not found])])
fi

dnl Offer --with-gmp-prefix.
AC_ARG_WITH(gmp-prefix,
            AC_HELP_STRING([--with-gmp-prefix=DIR],
                           [DIR Location of GMP package]),
                           [CPPFLAGS="${CPPFLAGS} -I$withval/include";
                            LDFLAGS="${LDFLAGS} -L$withval/lib";
                            LIBS="-lgmp $LIBS";
                            GMP_PREFIX="$withval";
                            ])

AC_CHECK_HEADER(gmp.h,
                [AC_CHECK_LIB(gmp,main,
                 [poly_cv_gmpfatal="no"],
                 [poly_cv_gmpfatal="yes"])],
                 [poly_cv_gmpfatal="yes"])

dnl Request arch flag
ARCHGIVEN="false"
AC_ARG_WITH(gcc-arch,
			AC_HELP_STRING([--with-gcc-arch=ARCH],
						   [ARCH Architecture to compile pluto for (march, mcpu etc.]),
						   [ARCHGIVEN="true";
						   ARCH="$withval";])

archopts=""
if test $ARCHGIVEN = "true"; then
    archopts="--with-gcc-arch=$ARCH"
fi

if test "$poly_cv_gmpfatal" = "yes"; then
    AC_MSG_ERROR([GMP not found])
    else
        AC_CHECK_DECLS(mp_get_memory_functions,[],[
                                                AC_LIBOBJ(mp_get_memory_functions)
                                                ],[#include <gmp.h>
                                                   ])
        fi



dnl Configuring PipLib (long long int is 64 bits)
echo ""
echo "=========================="
echo "Configuring PipLib"
echo "=========================="
configureopts="--with-gmp --with-bits=64 \
--with-gmp-prefix=$GMP_PREFIX \
--prefix=$prefix"
(cd piplib/
 $RM config.cache;
 ./configure ${configureopts} ${archopts}
 )


dnl Configuring Openscop
echo ""
echo "=========================="
echo "Configuring Openscop"
echo "=========================="
configureopts="--with-gmp=system \
--with-gmp-prefix=$GMP_PREFIX \
--prefix=$prefix"
(cd openscop/
 ./configure ${configureopts} ${archopts}
 )


dnl Configuring Clan-osl
echo ""
echo "=========================="
echo "Configuring Clan"
echo "=========================="
configureopts="--prefix=$prefix \
--with-osl=build \
--with-osl-builddir=$SOURCE_DIR/openscop"
(cd clan/
 ./configure ${configureopts} ${archopts}
 )

dnl check if clan configure failed
AC_CHECK_FILE([clan/Makefile], [], [AC_MSG_ERROR([configure in clan/ failed])])


dnl Configuring Candl-osl
dnl can't enable the gmp version without modifying
dnl the precision the is used elsewhere in pluto
echo ""
echo "=========================="
echo "Configuring Candl"
echo "=========================="
configureopts="--enable-llint-version \
--prefix=$prefix \
--with-piplib=build \
--with-piplib-builddir=$SOURCE_DIR/piplib \
--with-osl=build \
--with-osl-builddir=$SOURCE_DIR/openscop"
(cd candl/
 ./configure ${configureopts} ${archopts}
 )


dnl Configuring polylib
echo ""
echo "======================"
echo "Configuring polylib"
echo "======================"
configureopts="--enable-longlongint-lib \
--prefix=$prefix"
(cd polylib/
 ./configure ${configureopts}
)


AC_SUBST(ISL_INCLUDE)
AC_SUBST(ISL_LDFLAGS)
AC_SUBST(ISL_LIBADD)

dnl Provide --with-isl-prefix
AC_ARG_WITH(isl-prefix,
            AC_HELP_STRING([--with-isl-prefix=DIR],
              	             [DIR location where isl is installed]),
              [isl_prefix=$withval;
               ISL_INCLUDE="-I$withval/include";
               ISL_LDFLAGS="-L$withval/lib";
               ISL_LIBADD="$withval/lib/libisl.la";
               external_isl=true
               ],
              [external_isl=false])
AM_CONDITIONAL(EXTERNAL_ISL, test x$external_isl = xtrue)

if test x$external_isl = xfalse; then
dnl Configuring isl
dnl the current isl submodule still supports the piplib
dnl options, but the trunk of isl doesn't
echo ""
echo "======================"
echo "Configuring isl"
echo "======================"
configureopts="--with-gmp=system \
--with-gmp-prefix=$GMP_PREFIX \
--with-gmp-exec-prefix=$GMP_PREFIX \
--prefix=$prefix"
(cd isl/
 ./configure ${configureopts} ${archopts}
)
fi

dnl Offer --with-clang-prefix
AC_ARG_WITH(clang-prefix,
    AC_HELP_STRING([--with-clang-prefix=DIR], [DIR location of clang install]),
    [
        CLANG_PREFIX="$withval";
        external_clang=true
    ],
    [external_clang=false]
)

AM_CONDITIONAL(EXTERNAL_CLANG, test x$external_clang = xtrue)

dnl Configure PET
echo -e "\n=========================="
echo "Configuring PET"
echo "=========================="

if test x$external_clang = xtrue; then
    configureopts="CXXFLAGS=-fno-rtti --with-clang-prefix=$CLANG_PREFIX --with-clang-exec-prefix=$CLANG_PREFIX --with-isl-builddir=../isl --prefix=$prefix"
else
    configureopts="CXXFLAGS=-fno-rtti --with-isl-builddir=../isl --prefix=$prefix"
fi

(cd pet && ./configure ${configureopts})

dnl check if pet configure had failed
AC_CHECK_FILE([pet/Makefile], [], [AC_MSG_ERROR([configure in pet/ failed])])

dnl Configuring Cloog 0.14.1 (long long int)
echo ""
echo "=========================="
echo "Configuring Cloog-isl"
echo "=========================="
if test x$external_isl = xfalse; then
configureopts="--with-isl=build \
--with-isl-builddir=$SOURCE_DIR/isl \
--with-gmp=system \
--with-gmp-prefix=$GMP_PREFIX \
--with-osl=build \
--with-osl-builddir=$SOURCE_DIR/openscop \
--prefix=$prefix"
else
configureopts="--with-isl=system \
--with-isl-prefix=$isl_prefix \
--with-gmp=system \
--with-gmp-prefix=$GMP_PREFIX \
--with-osl=build \
--with-osl-builddir=$SOURCE_DIR/openscop \
--prefix=$prefix"
fi
(cd cloog-isl/
 ./configure ${configureopts} ${archopts}
 )

AC_PATH_PROGS(BASH, bash)

AC_CONFIG_FILES([getversion.sh], [chmod +x ./getversion.sh])
AC_CONFIG_COMMANDS([version.h],
                   [echo '#define PLUTO_VERSION "'`./getversion.sh`'"' > src/version.h])

AC_OUTPUT(Makefile \
          src/Makefile \
          polycc.sh
          )

echo " /*-----------------------------------------------*"
echo "  *           PLUTO configuration is OK           *"
echo "  *-----------------------------------------------*/"
echo "Your system is ready to compile PLUTO"
echo "Run 'make' next to compile"
echo "Then, run 'make test' to check if everything is working correctly"
echo "Use polycc (see README)"
