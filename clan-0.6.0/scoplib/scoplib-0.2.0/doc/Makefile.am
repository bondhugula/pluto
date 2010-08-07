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

MAINTAINERCLEANFILES	= Makefile.in

EXTRA_DIST		=			\
	htmldoc.tar.gz				\
	clan.pdf				\
	S4.cloog				\
	ascii-art.txt				\
	images/basic.eps			\
	images/basic.fig			\
	images/basic.jpg			\
	images/basic.pdf			\
	images/basic.txt			\
	images/basic1.eps			\
	images/basic1.fig			\
	images/basic1.jpg			\
	images/basic1.pdf			\
	images/basic1.txt			\
	images/basic2.eps			\
	images/basic2.fig			\
	images/basic2.jpg			\
	images/basic2.pdf			\
	images/basic2.txt			\
	images/tree.eps				\
	images/tree.fig				\
	images/tree.jpg				\
	images/tree.pdf				\
	images/tree.txt


info_TEXINFOS		= clan.texi

DOCDIR			= htmldoc
DOXYFILE		= Doxyfile
DOXYGEN			= @DOXYGEN@

# The following requires a fixed version of the Emacs 19.30 etags.
ETAGS_ARGS		= --lang=none \
 --regex='/^@node[ \t]+\([^,]+\)/\1/' $(srcdir)/doc/clan.texi
TAGS_DEPENDENCIES	= clan.texi $(DOXYFILE).in

htmldoc.tar.gz:
	mkdir -p $(top_srcdir)/doc/htmldoc
	$(DOXYGEN) $(DOXYFILE)
	tar czf htmldoc.tar.gz htmldoc
