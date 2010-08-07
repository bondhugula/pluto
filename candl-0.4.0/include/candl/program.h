
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                   program.h                             **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 9th 2003               **
    **--- |"-.-"| -------------------------------------------------------**
          |     |
          |     |
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2003-2008 Cedric Bastoul                                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * CAnDL, the Chunky Dependence Analyzer                                      *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


#ifndef CANDL_PROGRAM_H
# define CANDL_PROGRAM_H

# include <stdio.h>
# include <stdlib.h>
# include <candl/matrix.h>
# include <candl/statement.h>

# ifdef CANDL_SUPPORTS_SCOPLIB
#  include <scoplib/scop.h>
# endif


# if defined(__cplusplus)
extern "C"
  {
# endif

/**
 * candl_program_t structure:
 * this structure contains all the informations about a program.
 */
struct candl_program
{
  CandlMatrix * context;            /**< The context of the program. */
  int nb_statements;                /**< The number of statements. */
  CandlStatement ** statement;      /**< Array of nb_statements pointers on
                                      *   the statements of the program.
                                      */
  CandlMatrix    ** transformation; /**< Array of nb_statements pointers on
                                      *   the transformation candidate (one
                                      *   function per statement). If NULL:
                                      *   no tranformation candidate.
                                      */
  int* scalars_privatizable;
};
typedef struct candl_program CandlProgram;
typedef struct candl_program candl_program_t;
typedef struct candl_program * candl_program_p;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void candl_program_print_structure(FILE *, candl_program_p, int);
void candl_program_print(FILE *, candl_program_p);
void candl_program_print_candl_file(FILE *, candl_program_p);

/******************************************************************************
 *                         Memory alloc/dealloc function                      *
 ******************************************************************************/
candl_program_p candl_program_malloc();
void candl_program_free(candl_program_p);


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/
candl_program_p candl_program_read(FILE *);
/* This function is compiled if candl was configured with CLAN support. */
# ifdef CANDL_SUPPORTS_SCOPLIB
candl_program_p candl_program_read_scop(FILE *);
# endif

/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
/* This function is compiled if candl was configured with CLAN support. */
# ifdef CANDL_SUPPORTS_SCOPLIB
candl_program_p candl_program_convert_scop(scoplib_scop_p, int**);
# endif

# if defined(__cplusplus)
  }
# endif
#endif /* define CANDL_PROGRAM_H */

