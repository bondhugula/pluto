
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 vector.c                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 01/05/2008                     **
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
# include <ctype.h>
# include <clan/vector.h>


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * clan_vector_term function:
 * This function generates the vector representation of a term. It allocates
 * a vector with maximal size and put the term value at the right place
 * depending if the term is a constant, an iterator coefficient or a
 * parameter coefficient (see the structure of a PolyLib row if unsure!).
 * \param symbol      The first node of the list of symbols.
 * \param coefficient The constant or coefficient.
 * \param identifier  Identifier of iterator or parameter (NULL for constant).
 **
 * - 01/05/2008: first version.
 */
scoplib_vector_p
clan_vector_term(clan_symbol_p symbol, int coefficient, char * identifier)
{
  int rank, size;
  scoplib_vector_p vector;

  size = CLAN_MAX_DEPTH + CLAN_MAX_PARAMETERS + 2 ;
  vector = scoplib_vector_malloc(size);

  if (identifier == NULL)
    SCOPVAL_set_si(vector->p[size - 1],coefficient);
  else
  {
    rank = clan_symbol_get_rank(symbol,identifier);

    if (clan_symbol_get_type(symbol,identifier) == SCOPLIB_TYPE_ITERATOR)
      SCOPVAL_set_si(vector->p[rank],coefficient);
    else
      SCOPVAL_set_si(vector->p[CLAN_MAX_DEPTH + rank],coefficient);
  }
  return vector;
}
