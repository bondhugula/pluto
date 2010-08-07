
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                  violation.h                            **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: december 12th 2005               **
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


#ifndef CANDL_VIOLATION_H
# define CANDL_VIOLATION_H

# include <stdio.h>
# include <candl/dependence.h>
# include <candl/matrix.h>

# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * CandlViolation structure:
 * this structure contains all informations about a data dependence violation.
 */
struct candlviolation
{ CandlDependence * dependence; /**< Pointer to violated dependence. */
  int dimension;                /**< Violation dimension. */
  CandlMatrix * domain;         /**< Violation polyhedron. */
  struct candlviolation * next; /**< Pointer to next violation. */
};
typedef struct candlviolation CandlViolation;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void candl_violation_print_structure(FILE *, CandlViolation *, int);
void candl_violation_print(FILE *, CandlViolation *);
void candl_violation_pprint(FILE *, CandlViolation *);
void candl_violation_view(CandlViolation *);


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/
void candl_violation_free(CandlViolation *);


/******************************************************************************
 *                             Processing functions                           *
 ******************************************************************************/
CandlViolation * candl_violation_malloc();
void             candl_violation_add(CandlViolation **, CandlViolation **,
                                     CandlViolation *);
CandlViolation * candl_violation(CandlProgram *, CandlDependence *,
                                 CandlOptions *);


# if defined(__cplusplus)
  }
# endif
#endif /* define CANDL_VIOLATION_H */

