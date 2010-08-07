
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 matrix.h                              **
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


#ifndef SCOPLIB_MATRIX_H
# define SCOPLIB_MATRIX_H

# include <stdio.h>
# include <scoplib/macros.h>
# include <scoplib/vector.h>


# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The scoplib_matrix_t structure stores a matrix information in the PolyLib
 * format (the first entry of each row has a specific meaning). When a row
 * describes a linear constraint, a 0 means it is an equality == 0, a 1 means
 * an inequality >= 0. When a row describes an array access, a number different
 * than 0 is the array identifier (the remainder of the row describes the
 * access function of the first dimension of this array), otherwise it means
 * the row describes access functions for next array dimensions.
 */
struct scoplib_matrix
{
  unsigned NbRows;       /**< The number of rows */
  unsigned NbColumns;	 /**< The number of columns */
  scoplib_int_t ** p;    /**< An array of pointers to the beginning
			    of each row */
  scoplib_int_t * p_Init;/**< The matrix is stored here, contiguously
			    in memory */
  int p_Init_size;       /**< Needed to free the memory allocated by
			    mpz_init. */
};
typedef struct scoplib_matrix   scoplib_matrix_t;
typedef struct scoplib_matrix * scoplib_matrix_p;


/**
 * The scoplib_matrix_list_t structure describes a (chained) list of
 * matrices. It is used to store the list of matrices for the
 * iteration domain of a statement (possibly being a union of
 * convex domains).
 *
 */
struct scoplib_matrix_list
{
  scoplib_matrix_p elt;             /**< An element of the list. */
  struct scoplib_matrix_list* next; /**< Pointer to the next element
				       of the list.*/
};
typedef struct scoplib_matrix_list	scoplib_matrix_list_t;
typedef struct scoplib_matrix_list *	scoplib_matrix_list_p;


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void          scoplib_matrix_print_structure(FILE *, scoplib_matrix_p, int);
void          scoplib_matrix_print(FILE *, scoplib_matrix_p);
void          scoplib_matrix_print_dot_scop(FILE *, scoplib_matrix_p, int,
					    int, char **, int, char **,
					    int, char **);

void          scoplib_matrix_list_print_structure(FILE *,
						  scoplib_matrix_list_p, int);
void          scoplib_matrix_list_print(FILE *, scoplib_matrix_list_p);
void          scoplib_matrix_list_print_dot_scop(FILE *, scoplib_matrix_list_p,
						 int, int, char **, int,
						 char **, int, char **);


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/
scoplib_matrix_p	scoplib_matrix_read(FILE *);
scoplib_matrix_list_p	scoplib_matrix_list_read(FILE *);
scoplib_matrix_p	scoplib_matrix_read_arrays(FILE *, char ***, int *);


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
scoplib_matrix_p	scoplib_matrix_malloc(unsigned, unsigned);
void			scoplib_matrix_free_inside(scoplib_matrix_p);
void			scoplib_matrix_free(scoplib_matrix_p);

scoplib_matrix_list_p	scoplib_matrix_list_malloc();
void			scoplib_matrix_list_free(scoplib_matrix_list_p);


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
scoplib_matrix_p scoplib_matrix_ncopy(scoplib_matrix_p, int);
scoplib_matrix_p scoplib_matrix_copy(scoplib_matrix_p);
void	scoplib_matrix_replace_vector(scoplib_matrix_p, scoplib_vector_p, int);
void    scoplib_matrix_insert_vector(scoplib_matrix_p, scoplib_vector_p, int);
void	scoplib_matrix_add_vector(scoplib_matrix_p, scoplib_vector_p, int);
void	scoplib_matrix_sub_vector(scoplib_matrix_p, scoplib_vector_p, int);
scoplib_matrix_p scoplib_matrix_from_vector(scoplib_vector_p);
void    scoplib_matrix_replace_matrix(scoplib_matrix_p, scoplib_matrix_p, int);
void    scoplib_matrix_insert_matrix(scoplib_matrix_p, scoplib_matrix_p, int);
scoplib_matrix_p scoplib_matrix_concat(scoplib_matrix_p, scoplib_matrix_p);
int	scoplib_matrix_equal(scoplib_matrix_p, scoplib_matrix_p);    
    
# if defined(__cplusplus)
  }
# endif
#endif /* define SCOPLIB_MATRIX_H */
