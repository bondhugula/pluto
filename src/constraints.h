/*
 * PLuTo: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#ifndef _CONSTRAINTS_H
#define _CONSTRAINTS_H

#include "isl/set.h"

/* A system of linear inequalities and equalities; all inequalities in
 * the >= 0 form. The constant term is on the LHS as well, i.e.,
 *  c_1*x_1 + c_2*x_2 + ... + c_n*x_n + c_0 >= / = 0 */
typedef struct {
    /* Can be accessed as a double-subscripted array */
    int **val;

    /* Internal contigous buffer, val is set up to point into it */
    int *buf;

    /* Number of inequalities/equalities */
    int nrows;
    /* Number of columns (number of vars + 1) */
    int ncols;

    /* Is row i an inequality? 1 yes, 0 no */
    int *is_eq;

    /* Number of rows allocated a-priori */
    int alloc_nrows;
    int alloc_ncols;
} PlutoConstraints;


PlutoConstraints *pluto_constraints_alloc(int nrows, int ncols);
void pluto_constraints_free(PlutoConstraints *);
void pluto_constraints_resize(PlutoConstraints *, int, int);
PlutoConstraints *pluto_constraints_copy(PlutoConstraints *dest, const PlutoConstraints *src);

int best_elim_candidate(const PlutoConstraints *, int);
void fourier_motzkin_eliminate(PlutoConstraints *, int n);

PlutoMatrix *pluto2pip(const PlutoConstraints *, PlutoMatrix *pipmat);

PlutoConstraints *pluto_constraints_add(PlutoConstraints *, const PlutoConstraints *);
void pluto_constraints_simplify(PlutoConstraints *const cst);

int *pluto_constraints_solve(const PlutoConstraints *, int use_isl);
void pluto_constraints_add_inequality(PlutoConstraints *cst, int pos);
void pluto_constraints_add_equality(PlutoConstraints *cst, int pos);
void pluto_constraints_add_dim(PlutoConstraints *cst, int pos);
void pluto_constraints_remove_row(PlutoConstraints *, int);
void pluto_constraints_remove_dim(PlutoConstraints *, int);

void pluto_constraints_add_lb(PlutoConstraints *cst, int varnum, int lb);
void pluto_constraints_add_ub(PlutoConstraints *cst, int varnum, int ub);
void pluto_constraints_set_var(PlutoConstraints *cst, int varnum, int val);

void pluto_constraints_zero_row(PlutoConstraints *, int);
void pluto_constraints_normalize_row(PlutoConstraints *cst, int pos);

void pluto_constraints_print(FILE *fp, const PlutoConstraints *);
void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst);
void pluto_constraints_print_polylib(FILE *fp, const PlutoConstraints *cst);

/*
 * Construct a non-parametric basic set from the constraints in cst.
 */
__isl_give isl_basic_set *isl_basic_set_from_pluto_constraints(isl_ctx *ctx,
  const PlutoConstraints *cst);
#endif
