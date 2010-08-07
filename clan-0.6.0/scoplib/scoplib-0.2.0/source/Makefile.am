#
#   /**------- <| --------------------------------------------------------**
#    **         A                  Clan/Scop                              **
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
# *::8:::::::::SUNDOGa8a::. the Free Software Foundation, either version 2.1  *
# *::::::::8::::888:Y8888::                of the License, or (at your option)*
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
SUBDIRS			=

#############################################################################
MAINTAINERCLEANFILES	= Makefile.in

INCLUDES		= -I$(top_builddir) -I$(top_srcdir)	\
			  -I$(top_builddir)/include		\
		          -I$(top_srcdir)/include

#############################################################################

lib_LTLIBRARIES 	= libscoplib.la


libscoplib_la_SOURCES 	=			\
	scop.c					\
	statement.c				\
	matrix.c				\
	vector.c

AM_CFLAGS = -Wall -fomit-frame-pointer -g
