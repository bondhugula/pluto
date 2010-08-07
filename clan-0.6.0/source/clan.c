
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                  clan.c                               **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 30/04/2008                     **
    **- [""M# | #  U"U#U  -----------------------------------------------**
         | #  | #  \ .:/
         | #  | #___| #
 ******  | "--'     .-"  ******************************************************
 *     |"-"-"-"-"-#-#-##   Clan : the Chunky Loop Analyzer (experimental)     *
 ****  |     # ## ######  *****************************************************
 *      \       .::::'/                                                       *
 *       \      ::::'/     Copyright (C) 2008 Cedric Bastoul                  *
 *     :8a|    # # ##                                                         *
 *     ::88a      ###      This is free software; you can redistribute it     *
 *    ::::888a  8a ##::.   and/or modify it under the terms of the GNU Lesser *
 *  ::::::::888a88a[]:::   General Public License as published by the Free    *
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 3 of the       *
 *::::::::8::::888:Y8888:: License, or (at your option) any later version.    *
 *::::':::88::::888::Y88a::::::::::::...                                      *
 *::'::..    .   .....   ..   ...  .                                          *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with software; if not, write to the Free Software Foundation, Inc.,  *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Clan, the Chunky Loop Analyzer                                             *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <clan/scop.h>
#include <clan/options.h>


int main(int argc, char * argv[])
{
  scoplib_scop_p scop;
  clan_options_p options;
  FILE * input;
  FILE * output;

  /* Options and input/output file setting. */
  options = clan_options_read(argc,argv,&input,&output);

  /* Extraction of the polyhedral representation of the SCoP from the input. */
  if (options->inputscop)
    /* Input is a .scop file. */
    scop = scoplib_scop_read(input);
  else
    /* Input is a source code. */
    scop = clan_scop_extract(input,options);

  /* Printing of the internal data structure of the SCoP if asked. */
  if (options->structure)
    scoplib_scop_print(stdout,scop);

  /* Generation of the .scop output file. */
  int sopt = 0;
  if (options->castle)
    sopt |= SCOPLIB_SCOP_PRINT_CASTLE;
  if (options->arraystag)
    sopt |= SCOPLIB_SCOP_PRINT_ARRAYSTAG;
  scoplib_scop_print_dot_scop_options(output,scop,sopt);

  /* Save the planet. */
  clan_options_free(options);
  scoplib_scop_free(scop);

  return 0;
}
