/**
 * ISL-based operations for Pluto constraints, etc.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"

#include "isl/map.h"
#include "isl/set.h"
#include "isl/val_gmp.h"

/* start: 0-indexed */
void pluto_constraints_project_out_isl(
    PlutoConstraints *cst,
    int start,
    int num) {
    int end;
    isl_set *set;

    assert(num>= 0);
    if (num == 0)   return;

    end = start + num - 1;
    assert(start >= 0 && end <= cst->ncols-2);

    set = isl_set_from_pluto_constraints(cst, NULL);
    set = isl_set_project_out(set, isl_dim_set, start, num);
    pluto_constraints_free(cst);
    cst = isl_set_to_pluto_constraints(set);
    isl_set_free(set);
}

/* start: 0-indexed */
void pluto_constraints_project_out_isl_single(
    PlutoConstraints **cst,
    int start,
    int num) {
    int end;
    isl_basic_set *bset;

    assert(num >= 0);
    if (num == 0)   return;

    end = start + num - 1;
    assert(start >= 0 && end <= (*cst)->ncols-2);

    isl_ctx *ctx = isl_ctx_alloc();
    bset = isl_basic_set_from_pluto_constraints(ctx, *cst);
    bset = isl_basic_set_project_out(bset, isl_dim_set, start, num);
    pluto_constraints_free(*cst);
    *cst = isl_basic_set_to_pluto_constraints(bset);
    isl_basic_set_free(bset);
    isl_ctx_free(ctx);
}

/*
 * Construct a non-parametric basic set from the constraints in cst;
 * uses the first element in cst
 */
__isl_give isl_basic_set *isl_basic_set_from_pluto_constraints(
    isl_ctx *ctx, const PlutoConstraints *cst) {
    int i, j;
    int n_eq = 0, n_ineq = 0;
    isl_dim *dim;
    isl_mat *eq, *ineq;
    isl_basic_set *bset;

    for (i = 0; i < cst->nrows; ++i)
        if (cst->is_eq[i])
            n_eq++;
        else
            n_ineq++;

    eq = isl_mat_alloc(ctx, n_eq, cst->ncols);
    ineq = isl_mat_alloc(ctx, n_ineq, cst->ncols);

    dim = isl_dim_set_alloc(ctx, 0, cst->ncols - 1);

    n_eq = n_ineq = 0;
    for (i = 0; i < cst->nrows; ++i) {
        isl_mat **m;
        int row;

        if (cst->is_eq[i]) {
            m = &eq;
            row = n_eq++;
        } else {
            m = &ineq;
            row = n_ineq++;
        }

        for (j = 0; j < cst->ncols; ++j) {
            mpz_t tmp, one;
            mpz_init(tmp);
            mpz_init(one);
            mpz_set_ui(one, 1);
            mpz_set_sll(tmp, cst->val[i][j]);
            isl_val *v = isl_val_from_gmp(ctx, tmp, one);
            *m = isl_mat_set_element_val(*m, row, j, v);
            mpz_clear(tmp);
            mpz_clear(one);
        }
    }

    bset = isl_basic_set_from_constraint_matrices(dim, eq, ineq,
            isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
    return bset;
}

/*
 * Construct a non-parametric set from the constraints in cst
 */
__isl_give isl_set *isl_set_from_pluto_constraints(const PlutoConstraints *cst,
        isl_ctx *ctx) {
    isl_set *set;

    isl_space *dim = isl_dim_set_alloc(ctx, 0, cst->ncols - 1);
    set = isl_set_empty(dim);

    while (cst != NULL) {
        isl_basic_set *bset = isl_basic_set_from_pluto_constraints(ctx, cst);
        set = isl_set_union(set, isl_set_from_basic_set(bset));
        cst = cst->next;
    }

    return set;
}

static int extract_basic_set_constraints(__isl_take isl_basic_set *bset, void *usr) {
    PlutoConstraints **cst = (PlutoConstraints **) usr;

    PlutoConstraints *bcst = isl_basic_set_to_pluto_constraints(bset);
    isl_basic_set_free(bset);

    if (*cst == NULL) *cst = bcst;
    else {
        /* Careful - not to use pluto_constraints_union_isl (will lead to a
         * cyclic dependence */
        pluto_constraints_unionize_simple(*cst, bcst);
        pluto_constraints_free(bcst);
    }

    return 0;
}

/* Convert an isl_set to PlutoConstraints */
PlutoConstraints *isl_set_to_pluto_constraints(__isl_keep isl_set *set) {
    PlutoConstraints *cst = NULL;
    assert(set != NULL);
    isl_set_foreach_basic_set(set, &extract_basic_set_constraints, &cst);
    if (cst == NULL) cst = pluto_constraints_empty(isl_set_dim(set,
                               isl_dim_set)+ isl_set_dim(set, isl_dim_param)+ 1);
    return cst;
}


/*
 * Construct a non-parametric basic map from the constraints in cst
 */
__isl_give isl_basic_map *isl_basic_map_from_pluto_constraints(
    isl_ctx *ctx, const PlutoConstraints *cst, int n_in, int n_out,
    int n_par) {
    int i, j;
    int n_eq = 0, n_ineq = 0;
    isl_mat *eq, *ineq;
    isl_basic_map *bmap;
    isl_space *space;

    assert(cst->ncols == n_in + n_out + n_par + 1);

    space = isl_space_alloc(ctx, n_par, n_in, n_out);

    for (i = 0; i < cst->nrows; ++i) {
        if (cst->is_eq[i]) n_eq++;
        else n_ineq++;
    }

    eq = isl_mat_alloc(ctx, n_eq, cst->ncols);
    ineq = isl_mat_alloc(ctx, n_ineq, cst->ncols);

    n_eq = n_ineq = 0;
    for (i = 0; i < cst->nrows; ++i) {
        isl_mat **m;
        int row;

        if (cst->is_eq[i]) {
            m = &eq;
            row = n_eq++;
        } else {
            m = &ineq;
            row = n_ineq++;
        }

        for (j = 0; j < cst->ncols; ++j) {
            *m = isl_mat_set_element_si(*m, row, j, cst->val[i][j]);
        }
    }

    bmap = isl_basic_map_from_constraint_matrices(space, eq, ineq,
            isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_param, isl_dim_cst);
    return bmap;
}


/*
 * Convert an isl_basic_set to PlutoConstraints; drop any
 * existentially quantified variables (leads to an overapproximation in some
 * cases */
PlutoConstraints *isl_basic_set_to_pluto_constraints(
    __isl_keep isl_basic_set *bset) {
    int i, j;
    int eq_row;
    int ineq_row;
    int n_col;
    isl_mat *eq, *ineq;
    PlutoConstraints *cons;
    isl_basic_set *bset_c;

    bset_c = isl_basic_set_remove_divs(isl_basic_set_copy(bset));

    eq = isl_basic_set_equalities_matrix(bset_c,
                                         isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
    ineq = isl_basic_set_inequalities_matrix(bset_c,
            isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
    isl_basic_set_free(bset_c);

    eq_row = isl_mat_rows(eq);
    ineq_row = isl_mat_rows(ineq);
    n_col = isl_mat_cols(eq);

    cons = pluto_constraints_alloc(eq_row + ineq_row, n_col);
    cons->nrows = eq_row + ineq_row;

    for (i = 0; i < eq_row; ++i) {
        cons->is_eq[i] = 1;
        for (j = 0; j < n_col; ++j) {
            isl_val *v = isl_mat_get_element_val(eq, i, j);
            cons->val[i][j] = isl_val_get_num_ll(v);
            isl_val_free(v);
        }
    }

    for (i = 0; i < ineq_row; ++i) {
        cons->is_eq[eq_row+i] = 0;
        for (j = 0; j < n_col; ++j) {
            isl_val *v = isl_mat_get_element_val(ineq, i, j);
            cons->val[eq_row + i][j] = isl_val_get_num_ll(v);
            isl_val_free(v);
        }
    }

    isl_mat_free(eq);
    isl_mat_free(ineq);

    return cons;
}

/* Convert an isl_basic_map to a PlutoConstraints object */
PlutoConstraints *isl_basic_map_to_pluto_constraints(
    __isl_keep isl_basic_map *bmap) {
    PlutoConstraints *cst;

    isl_basic_map_to_pluto_constraints_func_arg(
        isl_basic_map_copy(bmap), &cst);

    return cst;
}

/* Convert an isl_basic_map to a PlutoConstraints object */
int isl_basic_map_to_pluto_constraints_func_arg(
    __isl_take isl_basic_map *bmap, void *user) {
    int i, j;
    int eq_row;
    int ineq_row;
    int n_col;
    isl_mat *eq, *ineq;
    PlutoConstraints *cons;

    eq = isl_basic_map_equalities_matrix(bmap,
                                         isl_dim_in, isl_dim_out, isl_dim_div, isl_dim_param, isl_dim_cst);
    ineq = isl_basic_map_inequalities_matrix(bmap,
            isl_dim_in, isl_dim_out, isl_dim_div, isl_dim_param, isl_dim_cst);

    eq_row = isl_mat_rows(eq);
    ineq_row = isl_mat_rows(ineq);
    n_col = isl_mat_cols(eq);

    cons = pluto_constraints_alloc(eq_row + ineq_row, n_col);
    cons->nrows = eq_row + ineq_row;

    for (i = 0; i < eq_row; ++i) {
        cons->is_eq[i] = 1;
        for (j = 0; j < n_col; ++j) {
            isl_val *v = isl_mat_get_element_val(eq, i, j);
            cons->val[i][j] = isl_val_get_num_ll(v);
            isl_val_free(v);
        }
    }

    for (i = 0; i < ineq_row; ++i) {
        cons->is_eq[eq_row+i] = 0;
        for (j = 0; j < n_col; ++j) {
            isl_val *v = isl_mat_get_element_val(ineq, i, j);
            cons->val[eq_row + i][j] = isl_val_get_num_ll(v);
            isl_val_free(v);
        }
    }

    isl_mat_free(eq);
    isl_mat_free(ineq);
    isl_basic_map_free(bmap);

    *(PlutoConstraints **)user = cons;
    return 0;
}


/* Use isl to solve these constraints (solves just for the first element if
 * it's a list of constraints */
int64 *pluto_constraints_lexmin_isl(const PlutoConstraints *cst, int negvar) {
    int i;
    int64 *sol;
    isl_ctx *ctx;
    isl_basic_set *bset, *all_positive;
    isl_set *domain, *all_positive_set, *lexmin;

    IF_DEBUG2(printf("[pluto] pluto_constraints_lexmin_isl (%d variables, %d constraints)\n",
                     cst->ncols-1, cst->nrows););

    ctx = isl_ctx_alloc();
    bset = isl_basic_set_from_pluto_constraints(ctx, cst);
    domain = isl_set_from_basic_set(bset);

    // Allow only positive values.
    if(negvar == 0) {
        all_positive = isl_basic_set_positive_orthant(isl_set_get_dim(domain));
        all_positive_set = isl_set_from_basic_set(all_positive);
        domain = isl_set_intersect(domain, all_positive_set);
    }
    // isl_set_print(domain, stdout, 0, ISL_FORMAT_ISL);
    // isl_set_dump(domain);
    lexmin = isl_set_lexmin(domain);

    if (isl_set_is_empty(lexmin)) {
        isl_set_free(lexmin);
        isl_ctx_free(ctx);
        return NULL;
    }

    int num_dimensions = isl_set_n_dim(lexmin);
    sol = (int64 *) malloc((num_dimensions)*sizeof(int64));

    // As the set is non parametric, there is only a single point in the set.
    // This point is the lexicographic minimum of the set.
    isl_point *p = isl_set_sample_point(lexmin);

    for (i = 0; i < num_dimensions; i++) {
        isl_val *v;
        v = isl_point_get_coordinate_val(p, isl_dim_set, i);
        sol[i] = isl_val_get_num_ll(v);
        isl_val_free(v);
    }

    isl_point_free(p);
    isl_ctx_free(ctx);

    return sol;
}

PlutoConstraints *pluto_constraints_union_isl(const PlutoConstraints *cst1,
        const PlutoConstraints *cst2) {
    isl_set *set1 = isl_set_from_pluto_constraints(cst1, NULL);
    isl_set *set2 = isl_set_from_pluto_constraints(cst2, NULL);
    isl_set *set3 = isl_set_union(set1, set2);

    PlutoConstraints *ucst = isl_set_to_pluto_constraints(set3);

    isl_set_free(set3);

    return ucst;
}


static int basic_map_count(__isl_take isl_basic_map *bmap, void *user) {
    int *count = user;

    *count += 1;
    isl_basic_map_free(bmap);
    return 0;
}


int isl_map_count(__isl_take isl_map *map, void *user) {
    int r;

    r = isl_map_foreach_basic_map(map, &basic_map_count, user);
    isl_map_free(map);
    return r;
}


PlutoConstraints *pluto_constraints_intersection_isl(const PlutoConstraints *cst1,
        const PlutoConstraints *cst2) {
    isl_set *iset1, *iset2, *iset3;
    PlutoConstraints *icst;

    isl_ctx *ctx = isl_ctx_alloc();

    iset1 = isl_set_from_pluto_constraints(cst1, ctx);
    iset2 = isl_set_from_pluto_constraints(cst2, ctx);

    iset3 = isl_set_intersect(iset1, iset2);

    icst = isl_set_to_pluto_constraints(iset3);
    isl_set_free(iset3);

    isl_ctx_free(ctx);

    return icst;
}


/* In-place intersection: first argument is modified */
PlutoConstraints *pluto_constraints_intersect_isl(PlutoConstraints *cst1,
        const PlutoConstraints *cst2) {
    PlutoConstraints *icst = pluto_constraints_intersection_isl(cst1, cst2);
    pluto_constraints_copy(cst1, icst);
    pluto_constraints_free(icst);

    return cst1;
}
