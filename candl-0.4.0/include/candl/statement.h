
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                   statement.h                           **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 8th 2003               **
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


#ifndef CANDL_STATEMENT_H
# define CANDL_STATEMENT_H

# include <stdio.h>
# include <candl/matrix.h>

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * CandlStatement structure:
 * this structure contains all the informations about a program statement.
 */
struct candlstatement
{ int label;                    /**< Statement number (it must be the array
                                  *   index of this statement in the "statement"
                                  *   field of the CandlProgram structure).
                                  */
  int type;                     /**< Statement type. */
  int depth;                    /**< Nesting level (nb of surrounding loops).*/
  int * index;                  /**< Iteration domain's iterator labels. */
  CandlMatrix * domain;         /**< Iteration domain. */
  CandlMatrix * written;        /**< Array of written data. */
  CandlMatrix * read;           /**< Array of read data. */
  void* ref;			/**< Reference to another structure
				   describing the same statement. */
};
typedef struct candlstatement CandlStatement;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void candl_statement_print_structure(FILE *, CandlStatement *, int);
void candl_statement_print(FILE *, CandlStatement *);


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/
void candl_statement_free(CandlStatement *);


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/
CandlStatement * candl_statement_read(FILE *, int, int);


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
CandlStatement * candl_statement_malloc();
int              candl_statement_commute(CandlStatement *, CandlStatement *);


# if defined(__cplusplus)
  }
# endif
#endif /* define CANDL_STATEMENT_H */

