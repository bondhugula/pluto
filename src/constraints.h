/*
 * PLuTo: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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

#ifdef GLPK
#include <glpk.h>
#endif

#include "isl/set.h"
#include "math_support.h"

#define int64 long long int

/* A system of linear inequalities and equalities; all inequalities in
 * the >= 0 form. The constant term is on the LHS as well, i.e.,
 *  c_1*x_1 + c_2*x_2 + ... + c_n*x_n + c_0 >= / = 0 */
struct pluto_constraints {
    /* Can be accessed as a double-subscripted array */
    int64 **val;

    /* Internal contiguous buffer, val is set up to point into it */
    int64 *buf;

    /* Number of inequalities/equalities */
    int nrows;
    /* Number of columns (number of vars + 1) */
    int ncols;

    /* Is row i an equality? 1 yes, 0 no */
    int *is_eq;

    /* Number of rows allocated a-priori */
    int alloc_nrows;
    int alloc_ncols;

    /* Names of the dimensions (optional) */
    char **names;

    struct pluto_constraints *next;
}; 

typedef struct pluto_constraints PlutoConstraints;

/* A constraint with one equality */
typedef PlutoConstraints PlutoEquality;
typedef PlutoConstraints Hyperplane;


PlutoConstraints *pluto_constraints_alloc(int nrows, int ncols);
void pluto_constraints_free(PlutoConstraints *);
PlutoConstraints *pluto_constraints_from_equalities(const PlutoMatrix *mat);
void pluto_constraints_resize(PlutoConstraints *, int, int);
void pluto_constraints_resize_single(PlutoConstraints *cst, int nrows, int ncols);
PlutoConstraints *pluto_constraints_copy(PlutoConstraints *dest, const PlutoConstraints *src);
PlutoConstraints *pluto_constraints_copy_single(PlutoConstraints *dest, const PlutoConstraints *src);
PlutoConstraints *pluto_constraints_dup(const PlutoConstraints *src);
PlutoConstraints *pluto_constraints_dup_single(const PlutoConstraints *src);

void fourier_motzkin_eliminate(PlutoConstraints *, int n);
void fourier_motzkin_eliminate_smart(PlutoConstraints *cst, int pos);

PlutoMatrix *pluto_constraints_to_pip_matrix(const PlutoConstraints *cst, PlutoMatrix *pmat);
PlutoConstraints *pluto_constraints_to_pure_inequalities_single(const PlutoConstraints *cst);
PlutoConstraints *pluto_constraints_from_inequalities(const PlutoMatrix *mat);

PlutoConstraints *pluto_constraints_add(PlutoConstraints *, const PlutoConstraints *);
PlutoConstraints *pluto_constraints_add_to_each(PlutoConstraints *cst1, const PlutoConstraints *cst2);

void pluto_constraints_simplify(PlutoConstraints *const cst);

int64 *pluto_constraints_lexmin(const PlutoConstraints *,int);
int64 *pluto_constraints_lexmin_isl(const PlutoConstraints *cst, int negvar);
int64 *pluto_constraints_lexmin_pip(const PlutoConstraints *cst, int negvar);
void pluto_constraints_add_inequality(PlutoConstraints *cst);
void pluto_constraints_add_equality(PlutoConstraints *cst);
void pluto_constraints_add_constraint(PlutoConstraints *cst, int is_eq);
void pluto_constraints_add_dim(PlutoConstraints *cst, int pos, const char *name);
void pluto_constraints_remove_row(PlutoConstraints *, int);
void pluto_constraints_remove_dim(PlutoConstraints *, int);

void pluto_constraints_add_lb(PlutoConstraints *cst, int varnum, int64 lb);
void pluto_constraints_add_ub(PlutoConstraints *cst, int varnum, int64 ub);
void pluto_constraints_set_var(PlutoConstraints *cst, int varnum, int64 val);

void pluto_constraints_zero_row(PlutoConstraints *, int);
void pluto_constraints_normalize_row(PlutoConstraints *cst, int pos);
PlutoConstraints *pluto_constraints_select_row(const PlutoConstraints *cst, int pos);
void pluto_constraints_negate_row(PlutoConstraints *cst, int pos);
void pluto_constraints_negate_constraint(PlutoConstraints *cst, int pos);
void pluto_constraints_interchange_cols(PlutoConstraints *cst, int col1, int col2);

PlutoConstraints *pluto_constraints_read(FILE *fp);

void pluto_constraints_print(FILE *fp, const PlutoConstraints *);
void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst);
void pluto_constraints_compact_print(FILE *fp, const PlutoConstraints *cst);
void pluto_constraints_print_polylib(FILE *fp, const PlutoConstraints *cst);
PlutoMatrix *pluto_constraints_to_matrix(const PlutoConstraints *cst);
PlutoConstraints *pluto_constraints_from_matrix(const PlutoMatrix *mat);
PlutoConstraints *pluto_constraints_from_mixed_matrix(const PlutoMatrix *mat, int *is_eq);
PlutoConstraints *pluto_constraints_image(const PlutoConstraints *cst, const PlutoMatrix *func);
void pluto_constraints_project_out(PlutoConstraints *cst, int start, int end);

void pluto_constraints_project_out_isl_single(
        PlutoConstraints **cst, 
        int start, 
        int num);


int pluto_constraints_num_in_list(const PlutoConstraints *const cst);

PlutoConstraints *pluto_constraints_intersection(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_intersect(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_difference(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_subtract(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_union(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_unionize(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_unionize_simple(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);

PlutoConstraints *pluto_constraints_intersect_isl(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);

int pluto_constraints_get_const_ub(const PlutoConstraints *cnst, int depth, int64 *ub);
int pluto_constraints_get_const_lb(const PlutoConstraints *cnst, int depth, int64 *lb);

int pluto_constraints_is_empty(const PlutoConstraints *cst);
int pluto_constraints_are_equal(const PlutoConstraints *cst1, const PlutoConstraints *cst2);

PlutoConstraints *pluto_constraints_empty(int ncols);
PlutoConstraints *pluto_constraints_universe(int ncols);
void pluto_constraints_set_names(PlutoConstraints *cst, char **names);
void pluto_constraints_set_names_range(PlutoConstraints *cst, char **names,
        int dest_offset, int src_offset, int num);

void pluto_constraints_set_names(PlutoConstraints *cst, char **names);
void pluto_constraints_set_names_range(PlutoConstraints *cst, char **names,
        int dest_offset, int src_offset, int num);

void print_polylib_visual_sets(char* name, PlutoConstraints *cst);
void print_polylib_visual_sets_new(char* name, PlutoConstraints *cst);

__isl_give isl_set *isl_set_from_pluto_constraints(const PlutoConstraints *cst, 
        isl_ctx *ctx);
PlutoConstraints *isl_set_to_pluto_constraints(__isl_keep isl_set *set);
__isl_give isl_basic_set *isl_basic_set_from_pluto_constraints(isl_ctx *ctx,
        const PlutoConstraints *cst);
PlutoConstraints *isl_basic_set_to_pluto_constraints(
        __isl_keep isl_basic_set *bset);
PlutoConstraints *isl_basic_map_to_pluto_constraints(
        __isl_keep isl_basic_map *bmap);
int isl_basic_map_to_pluto_constraints_func_arg(__isl_keep isl_basic_map *bmap, 
        void *user);
__isl_give isl_basic_map *isl_basic_map_from_pluto_constraints(
       isl_ctx *ctx, const PlutoConstraints *cst, int n_par, int n_in, int n_out);
void pluto_constraints_remove_names_single(PlutoConstraints *cst);

PlutoConstraints *pluto_constraints_unionize_isl(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);
PlutoConstraints *pluto_constraints_union_isl(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2);

PlutoMatrix *pluto_constraints_extract_equalities(const PlutoConstraints *cst);
int pluto_constraints_best_elim_candidate(const PlutoConstraints *cst, int max_elim);
int isl_map_count(__isl_take isl_map *map, void *user);
PlutoMatrix *isl_map_to_pluto_func(isl_map *map, int stmt_dim, int npar);
PlutoConstraints *pluto_hyperplane_get_negative_half_space(Hyperplane *h);
PlutoConstraints *pluto_hyperplane_get_non_negative_half_space(Hyperplane *h);
void pluto_constraints_shift_dim(PlutoConstraints *cst, int pos, PlutoMatrix *func);
void pluto_constraints_remove_names_single(PlutoConstraints *cst);

void pluto_constraints_remove_names_single(PlutoConstraints *cst);

void pluto_constraints_cplex_print(FILE *fp, const PlutoConstraints *cst);
PlutoConstraints *farkas_lemma_affine(const PlutoConstraints *dom, const PlutoMatrix *phi);
void pluto_constraints_gaussian_eliminate(PlutoConstraints *cst, int pos);

int pluto_constraints_get_num_non_zero_coeffs(const PlutoConstraints* cst);
#ifdef GLPK
int64* pluto_prog_constraints_lexmin_glpk(const PlutoConstraints *cst, PlutoMatrix* obj);
#endif
#endif
