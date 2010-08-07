
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                statement.h                            **
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
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 2.1 of the     *
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


#ifndef SCOPLIB_STATEMENT_H
# define SCOPLIB_STATEMENT_H

# include <stdio.h>
# include <scoplib/macros.h>
# include <scoplib/matrix.h>

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The scoplib_statement_t structure stores the useful informations of a given
 * statement to process it within a polyhedral framework.
 */
struct scoplib_statement
{
  scoplib_matrix_list_p domain;   /**< Iteration domain of the statement */
  scoplib_matrix_p schedule;      /**< Scheduling function for the statement */
  scoplib_matrix_p read;          /**< Array read access informations */
  scoplib_matrix_p write;         /**< Array write access informations */
  int nb_iterators;               /**< Original depth of the statement */
  char ** iterators;              /**< Array of (nb_iterators) iterator names */
  char * body;                    /**< Original statement body */


  /** Support for non-static code analysis (See Benabderrahmane's
      Research Report #6814). */
  int nb_exit_predicates;
  char ** exit_predicates;	/**< Array of exit predicats of all
				   while loops of the statement  */
  int nb_control_predicates;
  char ** control_predicates;	/**< Array of control predicats of all
				   irregular if of a statement  */


  struct scoplib_statement * next;/**< Next statement in the linked list */
};
typedef struct scoplib_statement   scoplib_statement_t;
typedef struct scoplib_statement * scoplib_statement_p;


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void             scoplib_statement_print_structure(FILE *, scoplib_statement_p,
						   int);
void             scoplib_statement_print(FILE *, scoplib_statement_p);
void             scoplib_statement_print_dot_scop(FILE *, scoplib_statement_p,
						  int, char **, int, char **);


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/
scoplib_statement_p scoplib_statement_read (FILE*, int, char***, int*);


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
scoplib_statement_p scoplib_statement_malloc();
void                scoplib_statement_free(scoplib_statement_p);


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
void             scoplib_statement_add(scoplib_statement_p *, scoplib_statement_p);
void             scoplib_statement_compact(scoplib_statement_p, int);
int              scoplib_statement_number(scoplib_statement_p);


# if defined(__cplusplus)
  }
# endif
#endif /* define SCOPLIB_STATEMENT_H */
