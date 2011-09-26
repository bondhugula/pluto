
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 access.h                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 14/06/2011                     **
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
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 2.1 of the     *
 *::::::::8::::888:Y8888:: License, or (at your option) any later version.    *
 *::::':::88::::888::Y88a::::::::::::...                                      *
 *::'::..    .   .....   ..   ...  .                                          *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							                                            *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with software; if not, write to the Free Software Foundation, Inc.,  *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Clan, the Chunky Loop Analyzer                                             *
 * Written by Cedric Bastoul,            Cedric.Bastoul@inria.fr              *
 *
 * Access support
 *            Prasanth Chatharasi,       prasanth@iith.ac.in                  *
 ******************************************************************************/
 
#ifndef SCOPLIB_ACCESS_H
# define SCOPLIB_ACCESS_H

# include <stdio.h>
# include <scoplib/scop.h>


# if defined(__cplusplus)
extern "C"
  {
# endif

/**
 * The scoplib_access_t structure stores the information regarding the access by
 * the symbol in a statement.
 */
 
struct scoplib_access
{
  scoplib_symbol_p symbol;    /* The information regarding symbol used in access*/
  scoplib_matrix_p matrix;    /* Accesses information by symbol */
};
typedef struct scoplib_access   scoplib_access_t;
typedef struct scoplib_access * scoplib_access_p;


/**
 * The scoplib_access_list_t structure describes a (chained) list of
 * accesses. It is used to store the list of accesess for the
 * read / writes  in the statements.
 */
 
struct scoplib_access_list
{
  scoplib_access_p elt;            /* An element in the list */
  struct scoplib_access_list* next;/*  Pointer to the next element in the list */
};
typedef struct scoplib_access_list	  scoplib_access_list_t;
typedef struct scoplib_access_list *	scoplib_access_list_p;

/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void scoplib_access_print(FILE * , scoplib_access_p );
void scoplib_access_list_print(FILE *,scoplib_access_list_p);
void scoplib_access_print_dot_scop(FILE * ,scoplib_scop_p, scoplib_access_p );
void scoplib_access_list_print_dot_scop(FILE *,scoplib_scop_p,scoplib_access_list_p);

/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
scoplib_access_p	scoplib_access_malloc(unsigned, unsigned);
void			scoplib_access_free(scoplib_access_p);

scoplib_access_list_p	scoplib_access_list_malloc();
void			scoplib_access_list_free(scoplib_access_list_p);

/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
 
scoplib_access_p  scoplib_access_copy(scoplib_access_p ) ;

scoplib_access_list_p  scoplib_access_list_copy_access_list(scoplib_access_list_p ) ;

void scoplib_access_list_add_access(scoplib_access_list_p* ,scoplib_access_p );

int scoplib_access_list_lookup_symbol(scoplib_access_list_p ,scoplib_symbol_p );


/* Converting the matrix to access format */
scoplib_access_list_p scoplib_access_matrix_to_access_format(scoplib_scop_p,
                                                scoplib_matrix_p ) ;
/* Getting the read_access from the statement */                                                
scoplib_access_list_p scoplib_access_get_read_access_list(scoplib_scop_p,scoplib_statement_p);
/* Getting the write_access from the statement */
scoplib_access_list_p scoplib_access_get_write_access_list(scoplib_scop_p,scoplib_statement_p);

# if defined(__cplusplus)
  }
# endif
#endif /* define SCOPLIB_ACCESS_H */

 
