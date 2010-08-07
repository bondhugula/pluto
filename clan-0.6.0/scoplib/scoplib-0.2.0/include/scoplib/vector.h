
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 vector.h                              **
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


#ifndef SCOPLIB_VECTOR_H
# define SCOPLIB_VECTOR_H

# include <stdio.h>
# include <scoplib/macros.h>


# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The scoplib_vector_t structure stores a vector information in the PolyLib
 * format (the first entry has a specific meaning). When a vector
 * describes a linear constraint, a 0 means it is an equality == 0, a 1 means
 * an inequality >= 0. When the vector describes an array access, a number
 * different than 0 is the array identifier.
 */
struct scoplib_vector
{
  unsigned Size;  /**< The number of vector entries */
  scoplib_int_t * p; /**< An array of values */
};
typedef struct scoplib_vector   scoplib_vector_t;
typedef struct scoplib_vector * scoplib_vector_p;


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void          scoplib_vector_print_structure(FILE *, scoplib_vector_p, int);
void          scoplib_vector_print(FILE *, scoplib_vector_p);


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/
scoplib_vector_p scoplib_vector_malloc(unsigned);
void		 scoplib_vector_free(scoplib_vector_p);


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
scoplib_vector_p scoplib_vector_add_scalar(scoplib_vector_p, int);
scoplib_vector_p scoplib_vector_add(scoplib_vector_p, scoplib_vector_p);
scoplib_vector_p scoplib_vector_sub(scoplib_vector_p, scoplib_vector_p);
void		 scoplib_vector_tag_inequality(scoplib_vector_p);
void		 scoplib_vector_tag_equality(scoplib_vector_p);

# if defined(__cplusplus)
  }
# endif
#endif /* define SCOPLIB_VECTOR_H */
