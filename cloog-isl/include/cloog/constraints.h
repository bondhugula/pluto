
   /**-------------------------------------------------------------------**
    **                               CLooG                               **
    **-------------------------------------------------------------------**
    **                           constraints.h                           **
    **-------------------------------------------------------------------**
    **                    First version: april 17th 2005                 **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2005 Cedric Bastoul                                          *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License as published by the Free Software Foundation; either               *
 * version 2.1 of the License, or (at your option) any later version.         *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public           *
 * License along with this library; if not, write to the Free Software        *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,                         *
 * Boston, MA  02110-1301  USA                                                *
 *                                                                            *
 * CLooG, the Chunky Loop Generator                                           *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


#ifndef CLOOG_CONSTRAINTS_H
#define CLOOG_CONSTRAINTS_H

#if defined(CLOOG_POLYLIB)
#include <cloog/polylib/matrix.h>
#elif defined(CLOOG_ISL)
#include <cloog/isl/matrix.h>
#else
struct cloogconstraintset;
typedef struct cloogconstraintset CloogConstraintSet;
struct cloogequalities;
typedef struct cloogequalities CloogEqualities;
#endif

#if defined(__cplusplus)
extern "C" 
  {
#endif 

/******************************************************************************
 *                        Equalities spreading functions                      *
 ******************************************************************************/
CloogEqualities *cloog_equal_alloc(int n, int nb_levels,
			int nb_parameters);
void		 cloog_equal_free(CloogEqualities *equal);
int              cloog_equal_count(CloogEqualities *equal);
int              cloog_equal_type(CloogEqualities *equal, int level);
void             cloog_equal_del(CloogEqualities *equal, int level);
int              cloog_equal_total_dimension(CloogEqualities *equal);

/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
void          cloog_constraint_set_normalize(CloogConstraintSet *, int);
void          cloog_constraint_set_free(CloogConstraintSet *);
int           cloog_constraint_set_contains_level(CloogConstraintSet *constraints,
			int level, int nb_parameters);
int           cloog_constraint_set_total_dimension(CloogConstraintSet *constraints);
int           cloog_constraint_set_n_iterators(CloogConstraintSet *constraints,
			int nb_parameters);
CloogConstraintSet *cloog_constraint_set_copy(CloogConstraintSet *);
CloogConstraintSet *cloog_constraint_set_simplify(CloogConstraintSet *, CloogEqualities *, int, int);

#if defined(CLOOG_POLYLIB) || defined(CLOOG_ISL)

CloogConstraintSet *cloog_constraint_set_for_reduction(CloogConstraint upper,
	       CloogConstraint lower);
CloogConstraintSet *cloog_constraint_set_reduce(CloogConstraintSet *constraints,
	int level, CloogEqualities *equal, int nb_par, cloog_int_t *bound);
CloogConstraint cloog_constraint_first(CloogConstraintSet *constraints);
int             cloog_constraint_is_valid(CloogConstraint constraint);
CloogConstraint cloog_constraint_next(CloogConstraint constraint);
CloogConstraint cloog_constraint_copy(CloogConstraint constraint);
void            cloog_constraint_release(CloogConstraint constraint);
CloogConstraint cloog_constraint_invalid(void);
int             cloog_constraint_total_dimension(CloogConstraint constraint);

CloogConstraint cloog_equal_constraint(CloogEqualities *equal, int j);
void            cloog_equal_add(CloogEqualities *equal,
				  CloogConstraintSet *constraints,
				  int level, CloogConstraint line, int nb_par);

CloogConstraint cloog_constraint_set_defining_equality(
			CloogConstraintSet *constraints, int level);
CloogConstraint cloog_constraint_set_defining_inequalities(
			CloogConstraintSet *constraints,
			int level, CloogConstraint *lower, int nb_parameters);
int           cloog_constraint_involves(CloogConstraint constraint, int v);
int           cloog_constraint_is_lower_bound(CloogConstraint constraint, int v);
int           cloog_constraint_is_upper_bound(CloogConstraint constraint, int v);
int           cloog_constraint_is_equality(CloogConstraint constraint);
void          cloog_constraint_constant_get(CloogConstraint constraint,
			cloog_int_t *val);
void          cloog_constraint_coefficient_get(CloogConstraint constraint,
			int var, cloog_int_t *val);
void          cloog_constraint_coefficient_set(CloogConstraint constraint,
			int var, cloog_int_t val);
void          cloog_constraint_clear(CloogConstraint constraint);
void          cloog_constraint_copy_coefficients(CloogConstraint constraint,
			cloog_int_t *dst);

struct clast_expr *cloog_constraint_variable_expr(CloogConstraint constraint,
			int level, CloogNames *names);

#endif

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
