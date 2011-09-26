
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan/Scoplib                        **
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
 * Written by Cedric Bastoul,            Cedric.Bastoul@inria.fr              *
 *
 * Symbol support
 *            Prasanth Chatharasi,       prasanth@iith.ac.in                  *
 *                                                                            *
 ******************************************************************************/


#ifndef SCOPLIB_SYMBOL_H
# define SCOPLIB_SYMBOL_H

# include <stdio.h>
# include <scoplib/macros.h>
# include <scoplib/matrix.h>
# include <scoplib/vector.h>

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The scoplib_symbol_t structure is a node of the symbol table of the parser.
 */
struct scoplib_symbol
{
  char *identifier;             /**< Symbol identifier */
  char *data_type;              /**<Symbol Datatype(int , float......)*/  
  int type;                     /**< Symbol type (variable, iterator...) */
  int flag;                     /**<Symbol Flag, 1- Newly Introduced in the code
                                    while optimization,0 - Already present in 
                                    original code */
  int num_of_dimensions;        /**<Symbol Dimensions ex: a[10][2] -- 2 */                                    
  scoplib_matrix_p dimensions_bounds;      /**<Symbol Dimensions bounds */                                    
  struct scoplib_symbol * next; /**< Next symbol in the symbol table */
};
typedef struct scoplib_symbol   scoplib_symbol_t;
typedef struct scoplib_symbol * scoplib_symbol_p;

/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void          scoplib_symbol_print_structure(FILE *, scoplib_symbol_p, int);
void          scoplib_symbol_print(FILE *, scoplib_symbol_p);
// Modified to print the dimensions_bounds properly
void          scoplib_symbol_print_dot_scop(FILE* ,scoplib_symbol_p,int nb_parameters,
                              char** parameters,int nb_arrays,char** arrays );


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
scoplib_symbol_p scoplib_symbol_malloc();
void          scoplib_symbol_free(scoplib_symbol_p);


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
scoplib_symbol_p scoplib_symbol_lookup(scoplib_symbol_p, char *);
scoplib_symbol_p scoplib_symbol_add(scoplib_symbol_p *, char *, int,scoplib_matrix_p,int);
scoplib_symbol_p scoplib_symbol_add_with_structure(scoplib_symbol_p *,scoplib_symbol_p);
scoplib_symbol_p scoplib_symbol_copy(scoplib_symbol_p );
void          scoplib_symbol_remove(scoplib_symbol_p*, scoplib_symbol_p);
int           scoplib_symbol_get_type(scoplib_symbol_p, char *);
void          scoplib_symbol_pad(scoplib_symbol_p,int );
void          scoplib_symbol_table_compact(scoplib_symbol_p , int,int);
char*         scoplib_symbol_table_get_bound(scoplib_symbol_p, int,  char** ,int );

# if defined(__cplusplus)
  }
# endif
#endif /* define SCOPLIB_SYMBOL_H */
