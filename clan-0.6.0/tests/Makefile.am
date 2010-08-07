#
#   /**------- <| --------------------------------------------------------**
#    **         A                     Clan                                **
#    **---     /.\   -----------------------------------------------------**
#    **   <|  [""M#                makefile.am                            **
#    **-   A   | #   -----------------------------------------------------**
#    **   /.\ [""M#         First version: 30/04/2008                     **
#    **- [""M# | #  U"U#U  -----------------------------------------------**
#         | #  | #  \ .:/
#         | #  | #___| #
# ******  | "--'     .-"  *****************************************************
# *     |"-"-"-"-"-#-#-##   Clan : the Chunky Loop Analyser (experimental)    *
# ****  |     # ## ######  ****************************************************
# *      \       .::::'/                                                      *
# *       \      ::::'/     Copyright (C) 2008 Cedric Bastoul                 *
# *     :8a|    # # ##                                                        *
# *     ::88a      ###      This is free software; you can redistribute it    *
# *    ::::888a  8a ##::.   and/or modify it under the terms of the GNU       *
# *  ::::::::888a88a[]:::   Lesser General Public License as published by     *
# *::8:::::::::SUNDOGa8a::. the Free Software Foundation, either version 3 of *
# *::::::::8::::888:Y8888::                the License, or (at your option)   *
# *::::':::88::::888::Y88a::::::::::::...  any later version.                 *
# *::'::..    .   .....   ..   ...  .                                         *
# * This software is distributed in the hope that it will be useful, but      *
# * WITHOUT ANY WARRANTY; without even the implied warranty of                *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General *
# * Public License  for more details.	                                      *
# *                                                                           *
# * You should have received a copy of the GNU Lesser General Public          *
# * License along with software; if not, write to the Free Software           *
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA  *
# *                                                                           *
# * Clan, the Chunky Loop Analyser                                            *
# * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                        *
# *                                                                           *
# *****************************************************************************/
#
# Makefile.am (or makefile if generated) of Clan, the Chunky Loop Analyser.
# Makefile.am is not a makefile, you must run the 'autogen.sh' THEN the
# configure shellscript to generate the Makefile thanks to this file.


#############################################################################
SUBDIRS 		=

#############################################################################

MAINTAINERCLEANFILES    = Makefile.in

TESTS_ENVIRONMENT       = top_builddir=$(top_builddir) SRCDIR=$(srcdir)

#############################################################################

EXTRA_DIST              = checker.sh

check_SCRIPTS		= check_unitary.sh check_complex.sh

## Replace by program names when needed
# check_PROGRAMS = ...
TESTS			= $(check_SCRIPTS)

TEST_FILES             	=			\
	$(top_srcdir)/tests/ax-do.c		\
	$(top_srcdir)/tests/gauss.c		\
	$(top_srcdir)/tests/gemver.c		\
	$(top_srcdir)/tests/iir_4_64.c		\
	$(top_srcdir)/tests/matmul.c		\
	$(top_srcdir)/tests/swim.c

UNITARY_TEST_FILES             	=					\
	$(top_srcdir)/tests/unitary/classicloop.c			\
	$(top_srcdir)/tests/unitary/simpleloop.c			\
	$(top_srcdir)/tests/unitary/basicnest.c				\
	$(top_srcdir)/tests/unitary/empty.c				\
	$(top_srcdir)/tests/unitary/simplecond.c			\
	$(top_srcdir)/tests/unitary/twoconds.c				\
	$(top_srcdir)/tests/unitary/twocondincond.c			\
	$(top_srcdir)/tests/unitary/complexcond.c			\
	$(top_srcdir)/tests/unitary/multiparenthcond.c			\
	$(top_srcdir)/tests/unitary/maxloopsimple.c			\
	$(top_srcdir)/tests/unitary/maxloopcomplex.c			\
	$(top_srcdir)/tests/unitary/minloopsimple.c			\
	$(top_srcdir)/tests/unitary/minloopcomplex.c			\
	$(top_srcdir)/tests/unitary/minmax.c				\
	$(top_srcdir)/tests/unitary/negassign.c				\
	$(top_srcdir)/tests/unitary/binop.c				\
	$(top_srcdir)/tests/unitary/complexop.c				\
	$(top_srcdir)/tests/unitary/complexprog.c			\
	$(top_srcdir)/tests/unitary/noloop.c				\
	$(top_srcdir)/tests/unitary/noloop2.c				\
	$(top_srcdir)/tests/unitary/ternary.c				\
	$(top_srcdir)/tests/unitary/functioncall.c			\
	$(top_srcdir)/tests/unitary/functioncall2.c			\
	$(top_srcdir)/tests/unitary/redefiter.c				\
	$(top_srcdir)/tests/unitary/comments.c				\
	$(top_srcdir)/tests/unitary/underscore.c			\
	$(top_srcdir)/tests/unitary/parenthesis.c			\
	$(top_srcdir)/tests/unitary/intaddit.c				\
	$(top_srcdir)/tests/unitary/assigniter.c			\
	$(top_srcdir)/tests/unitary/assigniter2.c			\
	$(top_srcdir)/tests/unitary/increment.c				\
	$(top_srcdir)/tests/unitary/field.c				\
	$(top_srcdir)/tests/unitary/object.c				\
	$(top_srcdir)/tests/unitary/nestedcurly.c


TEST_FILES_REF    =		 $(TEST_FILES:.c=.c.scop)
UNITARY_TEST_FILES_REF		= $(UNITARY_TEST_FILES:.c=.c.scop)

TEST_FILES_CLEAN		= $(TEST_FILES:.c=.c.test)
UNITARY_TEST_FILES_CLEAN	= $(UNITARY_TEST_FILES:.c=.c.test)


CLEANFILES			=		\
	$(TEST_FILES_CLEAN)			\
	$(UNITARY_TEST_FILES_CLEAN)


EXTRA_DIST	+=				\
	$(TEST_FILES)				\
	$(TEST_FILES_REF)			\
	$(UNITARY_TEST_FILES)			\
	$(UNITARY_TEST_FILES_REF)		\
	$(check_SCRIPTS)


TESTS_ENVIRONMENT +=					\
	CHECKER="$(top_srcdir)/tests/checker.sh"	\
	TEST_FILES="$(TEST_FILES)"			\
	UNITARY_TEST_FILES="$(UNITARY_TEST_FILES)"
