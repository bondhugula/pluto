
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 symbol.h                              **
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


#ifndef CLAN_SYMBOL_H
# define CLAN_SYMBOL_H

# include <stdio.h>
# include <scoplib/scop.h>
# include <scoplib/macros.h>

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The clan_symbol_t structure is a node of the symbol table of the parser.
 */
struct clan_symbol
{
  char * identifier;         /**< Symbol identifier */
  int type;                  /**< Symbol type (variable, iterator...) */
  int rank;                  /**< Depth for iterators, number for others */
  struct clan_symbol * next; /**< Next symbol in the symbol table */
};
typedef struct clan_symbol   clan_symbol_t;
typedef struct clan_symbol * clan_symbol_p;


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void          clan_symbol_print_structure(FILE *, clan_symbol_p, int);
void          clan_symbol_print(FILE *, clan_symbol_p);


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
clan_symbol_p clan_symbol_malloc();
void          clan_symbol_free(clan_symbol_p);


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
clan_symbol_p clan_symbol_lookup(clan_symbol_p, char *);
clan_symbol_p clan_symbol_add(clan_symbol_p *, char *, int, int);
void          clan_symbol_remove(clan_symbol_p*, clan_symbol_p);
int           clan_symbol_get_rank(clan_symbol_p, char *);
int           clan_symbol_get_type(clan_symbol_p, char *);
char **       clan_symbol_iterators(clan_symbol_p *, int);
char **       clan_symbol_id_array(clan_symbol_p, int, int *);


# if defined(__cplusplus)
  }
# endif
#endif /* define CLAN_SYMBOL_H */
