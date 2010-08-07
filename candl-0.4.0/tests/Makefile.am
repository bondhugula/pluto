#
#   /**-------------------------------------------------------------------**
#    **                              CAnDL                                **
#    **-------------------------------------------------------------------**
#    **                           Makefile.am                             **
#    **-------------------------------------------------------------------**
#    **                 First version: september 8th 2003                 **
#    **-------------------------------------------------------------------**/
#
#/*****************************************************************************
# *   CAnDL : the Chunky Analyser for Dependences in Loops (experimental)     *
# *****************************************************************************
# *                                                                           *
# * Copyright (C) 2003-2008 Cedric Bastoul                                    *
# *                                                                           *
# * This is free software; you can redistribute it and/or modify it under the *
# * terms of the GNU Lesser General Public License as published by the Free   *
# * Software Foundation; either version 3 of the License, or (at your option) *
# * any later version.							      *
# *                                                                           *
# * This software is distributed in the hope that it will be useful, but      *
# * WITHOUT ANY WARRANTY; without even the implied warranty of                *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General *
# * Public License for more details.                                          *
# *                                                                           *
# * You should have received a copy of the GNU Lesser General Public License  *
# * along with software; if not, write to the Free Software Foundation, Inc., *
# * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                    *
# *                                                                           *
# * CAnDL, the Chunky Dependence Analyser                                     *
# * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                        *
# *                                                                           *
# *****************************************************************************/


#############################################################################
SUBDIRS 		=

#############################################################################

MAINTAINERCLEANFILES    = Makefile.in

TESTS_ENVIRONMENT       = top_builddir=$(top_builddir) SRCDIR=$(srcdir)

#############################################################################

EXTRA_DIST              = checker.sh

check_SCRIPTS		= check_suite.sh

## Replace by program names when needed
# check_PROGRAMS = ...
TESTS			= $(check_SCRIPTS)

TEST_FILES             	=			\
	ax-do.scop				\
	gemver.scop				\
	matmul.scop				\
	swim.scop				\
	simple1.scop				\
	lu.scop					\
	jacobi1d.scop


TEST_FILES_DEP			= $(TEST_FILES:.scop=.dep)
TEST_FILES_CANDL		= $(TEST_FILES:.scop=.candl)
TEST_FILES_STRUCT		= $(TEST_FILES:.scop=.struct)
TEST_FILES_OPT			= $(TEST_FILES:.scop=.opt.scop)
TEST_FILES_SCAL			= $(TEST_FILES:.scop=.scalstruct)
TEST_FILES_DEPSCOP		= $(TEST_FILES:.scop=.depscop)

TEST_FILES_CLEAN_CANDL		= $(TEST_FILES:.scop=.candltest)
TEST_FILES_CLEAN_OPTSCOP	= $(TEST_FILES:.scop=.optscoptest)
TEST_FILES_CLEAN_STRUCT		= $(TEST_FILES:.scop=.structest)
TEST_FILES_CLEAN_SCAL		= $(TEST_FILES:.scop=.scaltest)
TEST_FILES_CLEAN_OUTSCOP	= $(TEST_FILES:.scop=.outscop)

BASE_TEST_FILES			= $(TEST_FILES:.scop=)

CLEANFILES			=		\
	$(TEST_FILES_CLEAN_CANDL)		\
	$(TEST_FILES_CLEAN_OPTSCOP)		\
	$(TEST_FILES_CLEAN_SCAL)		\
	$(TEST_FILES_CLEAN_STRUCT)		\
	$(TEST_FILES_CLEAN_OUTSCOP)


EXTRA_DIST	+=				\
	$(TEST_FILES)				\
	$(TEST_FILES_DEP)			\
	$(TEST_FILES_CANDL)			\
	$(TEST_FILES_STRUCT)			\
	$(TEST_FILES_OPT)			\
	$(TEST_FILES_SCAL)			\
	$(TEST_FILES_DEPSCOP)			\
	$(check_SCRIPTS)


TESTS_ENVIRONMENT +=				\
	TEST_FILES="$(BASE_TEST_FILES)"
