
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                statement.c                            **
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


# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <ctype.h>
# include <clan/statement.h>


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

/**
 * clan_statement_compact function:
 * This function scans the statement list to put the right number of columns
 * to every matrix (during construction we used CLAN_MAX_DEPTH and
 * CLAN_MAX_PARAMETERS to define matrix and vector sizes).
 * \param statement     The first statement to scan to compact matrices.
 * \param nb_parameters The true number of parameters in the SCoP.
 **
 * - 02/05/2008: first version.
 */
void
clan_statement_compact(scoplib_statement_p statement, int nb_parameters)
{
  int nb_iterators;
  scoplib_matrix_list_p tmp;

  while (statement != NULL)
  {
    nb_iterators = statement->nb_iterators;
    tmp = statement->domain;
    while (tmp)
      {
	clan_matrix_compact(tmp->elt,nb_iterators,nb_parameters);
	tmp = tmp->next;
      }
    clan_matrix_compact(statement->schedule,nb_iterators,nb_parameters);
    clan_matrix_compact(statement->read,nb_iterators,nb_parameters);
    clan_matrix_compact(statement->write,nb_iterators,nb_parameters);
    statement = statement->next;
  }
}
