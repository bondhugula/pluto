
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                    matrix.h                             **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: december 9th 2005                **
    **--- |"-.-"| -------------------------------------------------------**
          |     |
          |     |
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2005-2008 Cedric Bastoul                                     *
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


#ifndef CANDL_MATRIX_H
# define CANDL_MATRIX_H

# include <stdio.h>
# include <piplib/piplib.h>

# ifdef LINEAR_VALUE_IS_LONG
#  define CLAN_INT_T_IS_LONG
# endif
# ifdef LINEAR_VALUE_IS_LONGLONG
#  define CLAN_INT_T_IS_LONGLONG
# endif
# ifdef LINEAR_VALUE_IS_MP
#  define CLAN_INT_T_IS_MP
# endif

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * The matrix structure comes directly from PipLib (defined in piplib/piplib.h)
 * which is directly the PolyLib Matrix (defined in polylib/types.h)
 * here is how it looks like (at least in PipLib 1.3.5 version):
 *
 * struct pipmatrix
 * { unsigned NbRows;    // The number of rows (= NbConstraints in Polyhedron).
 *   unsigned NbColumns; // The number of columns (= Dimension+2 in Polyhedron).
 *   Value **p;          // An array of pointers to the beginning of each row.
 *   Value *p_Init;      // The matrix is stored here, contiguously in memory.
 *   int p_Init_size;    // Needed to free the memory allocated by mpz_init.
 * } ;
 * typedef struct pipmatrix PipMatrix ;
 */

typedef PipMatrix CandlMatrix;


/**
 * CandlMatrixList structure:
 * this structure reprensents a node of a linked list of CandlMatrix structures.
 */
struct candlmatrixlist
{ CandlMatrix * matrix;         /**< An element of the list. */
  struct candlmatrixlist * next;/**< Pointer to the next element of the list.*/
};
typedef struct candlmatrixlist CandlMatrixList;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void candl_matrix_print_structure(FILE *, CandlMatrix *, int);
void candl_matrix_print(FILE *, CandlMatrix *);
void candl_matrix_print_data(FILE *, CandlMatrix *);
void candl_matrix_list_print_structure(FILE *, CandlMatrixList *, int);
void candl_matrix_list_print(FILE *, CandlMatrixList *);

/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/
void candl_matrix_free(CandlMatrix *);
void candl_matrix_list_free(CandlMatrixList *);


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/
CandlMatrix * candl_matrix_read(FILE *);
CandlMatrixList * candl_matrix_list_read(FILE *);


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
CandlMatrix     * candl_matrix_malloc(int, int);
CandlMatrixList * candl_matrix_list_malloc();
CandlMatrix     * candl_matrix_violation(CandlMatrix  *, CandlMatrix  *,
                                         CandlMatrix  *, int, int);
int		candl_matrix_check_point (CandlMatrix* , CandlMatrix* );

# if defined(__cplusplus)
  }
# endif
#endif /* define CANDL_DEPENDENCE_H */

