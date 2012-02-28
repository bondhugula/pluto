/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
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
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
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

#ifdef PIP_WIDTH_MP
#include "piplib/piplibMP.h"
#else
#include "piplib/piplib64.h"
#endif

#define UB 0
#define LB 1
#define NB 2

PipMatrix *pip_matrix_populate(int **cst, int nrows, int ncols);

/*
 * Allocate with a max size of max_rows and max_cols;
 * initialized all to zero and is_eq to 0
 *
 * nrows set to 0, and ncols to max_cols;
 * As rows are added, increase nrows
 */
PlutoConstraints *pluto_constraints_alloc(int max_rows, int max_cols)
{
    PlutoConstraints *cst;

    cst = (PlutoConstraints *)malloc(sizeof(PlutoConstraints));

    int size = PLMAX(1,max_rows)*PLMAX(1,max_cols)*sizeof(int);

    cst->buf = (int *) malloc(size);

    if (cst->buf == NULL) {
        fprintf(stderr, "Not enough memory for allocating constraints\n");
        exit(1);
    }

    bzero(cst->buf, size);

    cst->is_eq = malloc(max_rows*sizeof(int));
    bzero(cst->is_eq, max_rows*sizeof(int));

    cst->val = (int **)malloc(max_rows*sizeof(int *));

    int i;
    for(i=0; i<max_rows; i++)   {
        cst->val[i] = &cst->buf[i*max_cols];
    }

    cst->alloc_nrows = max_rows;
    cst->alloc_ncols = max_cols;

    cst->nrows = 0;
    cst->ncols = max_cols;

    return cst;
}


void pluto_constraints_free(PlutoConstraints *cst)
{
    free(cst->buf);
    free(cst->val);
    free(cst->is_eq);
    free(cst);
}

/* Non-destructive resize */
void pluto_constraints_resize(PlutoConstraints *cst, int nrows, int ncols)
{
    int i, j;
    PlutoConstraints *newCst;

    assert(nrows >= 0 && ncols >= 0);

    newCst = pluto_constraints_alloc(PLMAX(nrows,cst->alloc_nrows), PLMAX(ncols,cst->alloc_ncols));

    newCst->nrows = nrows;
    newCst->ncols = ncols;

    for (i=0; i<PLMIN(newCst->nrows, cst->nrows); i++) {
        for (j=0; j<PLMIN(newCst->ncols, cst->ncols); j++) {
            newCst->val[i][j]= cst->val[i][j];
        }
        newCst->is_eq[i] = cst->is_eq[i];
    }

    free(cst->val);
    free(cst->buf);
    free(cst->is_eq);

    cst->nrows = nrows;
    cst->ncols = ncols;
    cst->alloc_nrows = newCst->alloc_nrows;
    cst->alloc_ncols = newCst->alloc_ncols;
    
    cst->val = newCst->val;
    cst->buf = newCst->buf;
    cst->is_eq = newCst->is_eq;

    free(newCst);
}


/* Adds cs1 and cs2 and puts them in cs1 -> returns cs1 itself */
/* Just adds to the first element if it's a list */
PlutoConstraints *pluto_constraints_add(PlutoConstraints *cst1, const PlutoConstraints *cst2)
{
    assert(cst1->ncols == cst2->ncols);

    if (cst1->nrows+cst2->nrows > cst1->alloc_nrows) {
        pluto_constraints_resize(cst1, cst1->nrows+cst2->nrows, cst1->ncols);
    }else{
        cst1->nrows = cst1->nrows + cst2->nrows;
    }

    int i;
    for (i=0; i < cst2->nrows; i++) {
        memcpy(cst1->val[cst1->nrows-cst2->nrows+i], 
                cst2->val[i], cst1->ncols*sizeof(int));
    }

    memcpy(&cst1->is_eq[cst1->nrows-cst2->nrows], cst2->is_eq, cst2->nrows*sizeof(int));

    cst1->ncols = cst1->ncols;

    return cst1;
}


/* Temporary data structure to compare two rows */
struct row_info {
    int *row;
    int ncols;
};

static int row_compar(const void *e1, const void *e2)
{
    int i, ncols, *row1, *row2;
    struct row_info *u1, *u2;

    u1 = *(struct row_info **)e1;
    u2 = *(struct row_info **)e2;
    row1 = u1->row;
    row2 = u2->row;
    ncols = u1->ncols;
    assert(ncols == u2->ncols);

    for (i=0; i<ncols; i++)  {
        if (row1[i] != row2[i]) break;
    }

    if (i==ncols) return 0;
    else if (row1[i] < row2[i])    {
        return -1;
    }else{
        return 1;
    }

    /* memcmp is inefficient compared to what's above when compiled with gcc */
    /* return memcmp(row1, row2, ncols*sizeof(int)); */
}


/* 
 * Eliminates duplicate constraints; the simplified constraints
 * are still at the same memory location but the number of constraints 
 * in it will decrease
 */
void pluto_constraints_simplify(PlutoConstraints *const cst)
{
    int i, j, p, _gcd;

    if (cst->nrows == 0)    {
        return;
    }

    PlutoConstraints *tmpcst = pluto_constraints_alloc(cst->nrows, cst->ncols);

    int *is_redun = (int *) malloc(sizeof(int)*cst->nrows);
    bzero(is_redun, cst->nrows*sizeof(int));

    /* Normalize cst - will help find redundancy */
    for (i=0; i<cst->nrows; i++)   {
        assert(cst->is_eq[i] == 0);
        _gcd = -1;
        for (j=0; j<cst->ncols; j++)  {
            if (cst->val[i][j] != 0)   {
                if (_gcd == -1) _gcd = PLABS(cst->val[i][j]);
                else _gcd = gcd(PLABS(cst->val[i][j]), _gcd);
            }
        }

        if (_gcd != -1 && _gcd != 1) {
            /* Normalize by gcd */
            for (j=0; j< cst->ncols; j++)   {
                cst->val[i][j] /= _gcd;
            }
        }
    }

    struct row_info **rows;
    rows = (struct row_info **) malloc(cst->nrows*sizeof(struct row_info *));
    for (i=0; i<cst->nrows; i++) {
        rows[i] = (struct row_info *) malloc(sizeof(struct row_info));
        rows[i]->row = cst->val[i];
        rows[i]->ncols = cst->ncols;
    }
    qsort(rows, cst->nrows, sizeof(struct row_info *), row_compar);

    for (i=0; i<cst->nrows; i++) {
        cst->val[i] = rows[i]->row;
        free(rows[i]);
    }
    free(rows);

    is_redun[0] = 0;
    for (i=1; i<cst->nrows; i++)    {
        is_redun[i] = 0;

        for (j=0; j<cst->ncols; j++)    {
            if (cst->val[i-1][j] != cst->val[i][j]) break;
        }

        if (j==cst->ncols) is_redun[i] = 1;

        for (j=0; j<cst->ncols; j++)    {
            if (cst->val[i][j] != 0) break;
        }

        if (j==cst->ncols)  {
            /* All zeros */
            is_redun[i] = 1;
        }
    }

    p = 0;
    for (i=0; i<cst->nrows; i++)    {
        if (!is_redun[i])  {
            /*
            for (j=0; j<cst->ncols; j++)    {
                tmpcst->val[p][j] = cst->val[i][j];
            }
            */
            memcpy(tmpcst->val[p], cst->val[i], cst->ncols*sizeof(int));
            p++;
        }
    }
    tmpcst->nrows = p;
    tmpcst->ncols = cst->ncols;

    pluto_constraints_copy(cst, tmpcst);
    pluto_constraints_free(tmpcst);

    free(is_redun);
}


/* 
 * Eliminates the pos^th variable, where pos has to be between 0 and cst->ncols-2;
 * Remember that the last column is for the constant. The implementation does not 
 * have a redundancy check; it just  eliminates duplicates after gcd normalization
 * cst will be resized if necessary
 */
void fourier_motzkin_eliminate(PlutoConstraints *cst, int pos)
{
    int i, j, k, l, p, q;
    int lb, ub, nb;
    int **csm;
    int *bound;

    for (i=0; i<cst->nrows; i++) {
        /* Does not support equalities: please extend if necessary */
        assert(cst->is_eq[i] != 1);
    }

    PlutoConstraints *newcst;

    // newcst = pluto_constraints_alloc(cst->nrows*cst->nrows/4, cst->ncols);

    assert(pos >= 0);
    assert(pos <= cst->ncols-2);
	// At least one variable
	assert(cst->ncols >= 2);

    csm = cst->val;

    for (i=0; i<cst->nrows; i++)    {
        if (csm[i][pos] != 0) break;
    }

    if (i==cst->nrows) {

        newcst = pluto_constraints_alloc(cst->nrows, cst->ncols);
        int **newcsm = newcst->val;

        for (j=0; j<cst->nrows; j++)    {
            q=0;
            for (k=0; k<cst->ncols; k++){
                if (k!=pos)   {
                    newcsm[j][q] = csm[j][k];
                    q++;
                }
            }
        }
        newcst->nrows = cst->nrows;
        newcst->ncols = cst->ncols-1;
    }else{
        bound = (int *) malloc(cst->nrows*sizeof(int));
        // printf("bound size: %lf\pos", ceil((cst->nrows*cst->nrows)/4));

        lb=0;
        ub=0;
        nb=0;
        /* Variable does appear */
        for (j=0; j<cst->nrows; j++)    {
            if (csm[j][pos] == 0) {
                bound[j] = NB;
                nb++;
            }else if (csm[j][pos] > 0) {
                // printf("accessing: %d\pos", j);
                bound[j] = LB;
                lb++;
            }
            else{
                bound[j] = UB;
                ub++;
            }
        }
        newcst = pluto_constraints_alloc(lb*ub+nb, cst->ncols);
        int **newcsm = newcst->val;

        p=0;
        for (j=0; j<cst->nrows; j++)    {
            if (bound[j] == UB) {
                for (k=0; k<cst->nrows; k++)    {
                    if (bound[k] == LB) {
                        q = 0;
                        for(l=0; l < cst->ncols; l++)  {
                            if (l!=pos)   {
                                newcsm[p][q] = 
                                    csm[j][l]*(lcm(csm[k][pos], 
                                                -csm[j][pos])/(-csm[j][pos])) 
                                    + csm[k][l]*(lcm(-csm[j][pos], 
                                                csm[k][pos])/csm[k][pos]); 
                                q++;
                            }
                        }
                        p++;
                    }
                }
            }else if (bound[j] == NB)   {
                q = 0;
                for (l=0; l<cst->ncols; l++)    {
                    if (l!=pos)   {
                        newcsm[p][q] = csm[j][l];
                        q++;
                    }
                }
                p++;
            }
        }
        newcst->nrows = p;
        newcst->ncols = cst->ncols-1;
        free(bound);
    }

    pluto_constraints_simplify(newcst);
    pluto_constraints_copy(cst, newcst);
    pluto_constraints_free(newcst);
}


/* Copy constraints from src into dest; if dest does not have enough space,
 * resize it */
PlutoConstraints *pluto_constraints_copy(PlutoConstraints *dest, const PlutoConstraints *src)
{
    int i;

    if (src->nrows > dest->alloc_nrows || src->ncols > dest->alloc_ncols) {
        pluto_constraints_resize(dest, PLMAX(src->nrows,dest->alloc_nrows), 
                PLMAX(src->ncols,dest->alloc_ncols));
    }   

    dest->nrows = src->nrows;
    dest->ncols = src->ncols;

    for (i=0; i<dest->nrows; i++) {
        memcpy(dest->val[i], src->val[i], src->ncols*sizeof(int));
    }

    memcpy(dest->is_eq, src->is_eq, dest->nrows*sizeof(int));

    return dest;
}



/* Duplicate constraints; returned constraints should be freed with
 * pluto_constraints_free */
PlutoConstraints *pluto_constraints_dup(const PlutoConstraints *src)
{
    assert(src != NULL);
    PlutoConstraints *dup = pluto_constraints_alloc(src->alloc_nrows, src->alloc_ncols);

    pluto_constraints_copy(dup, src);

    return dup;
}


void pluto_constraints_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j;

    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols);

    for (i=0; i<cst->nrows; i++) {
        for (j=0; j<cst->ncols; j++) {
            fprintf(fp, "%s%d ", cst->val[i][j]>=0? " ":"", cst->val[i][j]);
        }
        fprintf(fp, "\t %s 0\n", cst->is_eq[i]? "==": ">=");
    }
}


/* Print in polylib format */
void pluto_constraints_print_polylib(FILE *fp, const PlutoConstraints *const cst)
{
    int i, j;

    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols+1);

    for (i=0; i<cst->nrows; i++)    {
        fprintf(fp, "%s ", cst->is_eq[i]? "0": "1");
        for (j=0; j<cst->ncols; j++)    {
            fprintf(fp, "%d ", cst->val[i][j]);
        }
        fprintf(fp, "\n");
    }
}


/* Converts to polylib style matrix */
PlutoMatrix *pluto_constraints_to_matrix(const PlutoConstraints *cst)
{
    int i, j;
    PlutoMatrix *mat;

    mat = pluto_matrix_alloc(cst->nrows, cst->ncols+1);

    for (i=0; i<cst->nrows; i++)    {
        mat->val[i][0] = (cst->is_eq[i])? 0: 1;
        for (j=0; j<cst->ncols; j++)    {
            mat->val[i][j+1] = cst->val[i][j];
        }
    }

    return mat;
}



void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j, var;

    int nrows = cst->nrows;
    int ncols = cst->ncols;

    for (i=0; i<nrows; i++) {
        var = 'i';
        for (j=0; j<ncols; j++)    {
            if (j==ncols-1) {
                /* constant */
                fprintf(fp, "%s%d ", (cst->val[i][j]>=0)?"+":"", cst->val[i][j]);
            }else{
                if (cst->val[i][j] > 0)
                    fprintf(fp, "+%dc_%c ", cst->val[i][j], var);
                else if (cst->val[i][j] < 0)
                    fprintf(fp, "%dc_%c ", cst->val[i][j], var);
                else fprintf(fp, "      ");
            }
            var++;
        }
        fprintf(fp, "%s 0\n", cst->is_eq[i]? "==": ">=");
    }
    fprintf(fp, "\n");
}


/* Convert Pluto constraints into PIP format (first column is
 * 0/1 based on equality/inequality */
PlutoMatrix *pluto_constraints_to_pip_matrix(const PlutoConstraints *cst, PlutoMatrix *pmat)
{
    int i, j;

    assert(pmat != NULL);

    for (i=0; i<cst->nrows; i++)    {
        /* Equality or Inequality */
        pmat->val[i][0] = cst->is_eq[i]? 0: 1;
        for (j=1; j< cst->ncols+1; j++)    {
            pmat->val[i][j] = cst->val[i][j-1];
        }
    }
    pmat->nrows = cst->nrows;
    pmat->ncols = cst->ncols+1;

    return pmat;
}

/*
 * Construct a non-parametric basic set from the constraints in cst.
 */
__isl_give isl_basic_set *isl_basic_set_from_pluto_constraints(
       isl_ctx *ctx, const PlutoConstraints *cst)
{
    int i, j;
    int n_eq = 0, n_ineq = 0;
    isl_int v;
    isl_dim *dim;
    isl_mat *eq, *ineq;
    isl_basic_set *bset;

    isl_int_init(v);

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
            isl_int_set_si(v, cst->val[i][j]);
            *m = isl_mat_set_element(*m, row, j, v);
        }
    }

    isl_int_clear(v);

    bset = isl_basic_set_from_constraint_matrices(dim, eq, ineq,
                isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
    return bset;
}

/* Use isl to solve these constraints */
int *pluto_constraints_solve_isl(const PlutoConstraints *cst) {
  int i, *sol;
  isl_ctx *ctx;
  isl_basic_set *bset, *all_positive;
  isl_set *domain, *all_positive_set, *lexmin;

  ctx = isl_ctx_alloc();
  bset = isl_basic_set_from_pluto_constraints(ctx, cst);
  domain = isl_set_from_basic_set(bset);

  // Allow only positive values.
  all_positive = isl_basic_set_positive_orthant(isl_set_get_dim(domain));
  all_positive_set = isl_set_from_basic_set(all_positive);
  domain = isl_set_intersect(domain, all_positive_set);

  // isl_set_print(domain, stdout, 0, ISL_FORMAT_ISL);
  lexmin = isl_set_lexmin(domain);

  if (isl_set_is_empty(lexmin))
    return NULL;

  int num_dimensions = isl_set_n_dim(lexmin);
  sol = (int *) malloc((num_dimensions)*sizeof(int));

  // As the set is non parametric, there is only a single point in the set.
  // This point is the lexicographic minimum of the set.
  isl_point *p = isl_set_sample_point(lexmin);

  for (i = 0; i < num_dimensions; i++) {
    isl_int v;
    isl_int_init(v);
    isl_point_get_coordinate(p, isl_dim_set, i, &v);
    sol[i] = isl_int_get_si(v);
    isl_int_clear(v);
  }

  isl_point_free(p);
  isl_ctx_free(ctx);

  return sol;
}

/* Use PIP to solve these constraints */
int *pluto_constraints_solve_pip(const PlutoConstraints *cst)
{
    int bignum, i;
    PipMatrix  *domain, *context;
    PipQuast   *solution;
    PipOptions *pipOptions;
    PipList *listPtr;
    int *sol;
    PlutoMatrix *pipmat;

        pipmat = pluto_matrix_alloc(cst->nrows, cst->ncols+1);

    /* Convert constraints to PIP format */
    pluto_constraints_to_pip_matrix(cst, pipmat);

    /* First column says whether it's an inequality and the last is for the constant;
     * so we have ncols-2 variables */

    context = NULL;
    bignum = -1;

    domain = pip_matrix_populate(pipmat->val, pipmat->nrows, 
            pipmat->ncols);
    pipOptions = pip_options_init();

    /* IF_DEBUG2(fprintf(stdout, "Calling PIP on a %dx%d formulation\n", 
                pipmat->nrows, pipmat->ncols)); */

    solution = pip_solve(domain,context,bignum,pipOptions) ;

    // IF_DEBUG2(pip_quast_print(stdout, solution, 0));

    assert(solution->condition==NULL);

    listPtr = solution->list;

    sol = NULL;
    if (listPtr != NULL)    {
        sol = (int *) malloc((pipmat->ncols-2)*sizeof(int));

        for (i=0; i<pipmat->ncols-2; i++)    {
            /* This is just a lexmin and not a parametrix lexmin and so each
             * vector in the list is actually a constant */
#ifdef PIP_WIDTH_MP
            sol[i] = mpz_get_si(*listPtr->vector->the_vector);
#else
            sol[i] = (int) *listPtr->vector->the_vector;
#endif
            listPtr = listPtr->next;
        }
    }

    pip_options_free(pipOptions);
    pip_matrix_free(domain);
    pip_quast_free(solution);

    pluto_matrix_free(pipmat);

    return sol;
}

/* Solve these constraints */
int *pluto_constraints_solve(const PlutoConstraints *cst) {
    if (options->islsolve) {
        return pluto_constraints_solve_isl(cst);
    }else{
        return pluto_constraints_solve_pip(cst);
    }
}


/* All equalities */
PlutoConstraints *pluto_constraints_from_equalities(const PlutoMatrix *mat)
{
    int i, j;

    PlutoConstraints *cst; 

    cst = pluto_constraints_alloc(mat->nrows, mat->ncols);

    for (i=0; i<mat->nrows; i++)    {
        cst->is_eq[i] = 1;
        for (j=0; j<mat->ncols; j++)    {
            cst->val[i][j] = mat->val[i][j];
        }
    }
    cst->nrows = mat->nrows;

    return cst;
}

/* Add an inequality (>= 0); initialize it to all zero */
void pluto_constraints_add_inequality(PlutoConstraints *cst, int pos)
{
    int i, j;

    assert(pos >= 0 && pos <= cst->nrows);

    if (cst->nrows == cst->alloc_nrows)   {
        pluto_constraints_resize(cst, cst->nrows+1, cst->ncols);
    }else{
        cst->nrows++;
    }

    for (i=cst->nrows-2; i>=pos; i--) {
        for (j=0; j<cst->ncols; j++) {
            cst->val[i+1][j] = cst->val[i][j];
        }
        cst->is_eq[i+1] = cst->is_eq[i];
    }

    for (j=0; j<cst->ncols; j++) {
        cst->val[pos][j] = 0;
    }
    cst->is_eq[pos] = 0;

}


/* Add an equality; initialize it to all zero */
void pluto_constraints_add_equality(PlutoConstraints *cst, int pos)
{
    int i, j;

    assert(pos >= 0 && pos <= cst->nrows);
    assert(cst->nrows <= cst->alloc_nrows);

    if (cst->nrows == cst->alloc_nrows)   {
        pluto_constraints_resize(cst, cst->nrows+1, cst->ncols);
    }else{
        cst->nrows++;
    }

    for (i=cst->nrows-2; i>=pos; i--) {
        for (j=0; j<cst->ncols; j++) {
            cst->val[i+1][j] = cst->val[i][j];
        }
        cst->is_eq[i+1] = cst->is_eq[i];
    }

    for (j=0; j<cst->ncols; j++) {
        cst->val[pos][j] = 0;
    }
    cst->is_eq[pos] = 1;
}



/* Remove a row; pos is 0-indexed */
void pluto_constraints_remove_row(PlutoConstraints *cst, int pos) 
{
    int i, j;

    assert(pos >= 0 && pos <= cst->nrows-1);

    for (i=pos; i<cst->nrows-1; i++) {
        for (j=0; j<cst->ncols; j++) {
            cst->val[i][j] = cst->val[i+1][j];
        }
        cst->is_eq[i] = cst->is_eq[i+1];
    }
    cst->nrows--;
}


/* Remove a variable */
void pluto_constraints_remove_dim(PlutoConstraints *cst, int pos)
{
    int i, j;

    assert(pos >= 0 && pos <= cst->ncols-2);

    for (j=pos; j<cst->ncols-1; j++) {
        for (i=0; i<cst->nrows; i++) {
            cst->val[i][j] = cst->val[i][j+1];
        }
    }
    cst->ncols--;
}


/* Blank 'pos' th row */
void pluto_constraints_zero_row(PlutoConstraints *cst, int pos)
{
    int j;

    assert(pos >= 0 && pos <= cst->nrows-1);

    for (j=0; j<cst->ncols; j++)    {
        cst->val[pos][j] = 0;
    }
}


/* Normalize row by its gcd */
void pluto_constraints_normalize_row(PlutoConstraints *cst, int pos)
{
    int i, j, k;

    assert(pos >= 0 && pos <= cst->nrows-1);

    /* Normalize cst first */
    for (i=0; i<cst->nrows; i++)    {
        if (cst->val[i][0] == 0) continue;
        int rowgcd = abs(cst->val[i][0]);
        for (j=1; j<cst->ncols; j++)    {
            if (cst->val[i][j] == 0)  break;
            rowgcd = gcd(rowgcd,abs(cst->val[i][j]));
        }
        if (i == cst->nrows)   {
            if (rowgcd > 1)    {
                for (k=0; k<cst->ncols; k++)    {
                    cst->val[i][k] /= rowgcd;
                }
            }
        }
    }
}


/* Add a variable; resize if necessary; initialize to zero */
void pluto_constraints_add_dim(PlutoConstraints *cst, int pos)
{
    int i, j;

    /* Has to be a new variable */
    assert(pos >= 0 && pos <= cst->ncols-1);

    if (cst->ncols == cst->alloc_ncols)  {
        pluto_constraints_resize(cst, cst->nrows, cst->ncols+1);
    }else{
        cst->ncols++;
    }

    for (j=cst->ncols-2; j>=pos; j--) {
        for (i=0; i<cst->nrows; i++) {
            cst->val[i][j+1] = cst->val[i][j];
        }
    }

    /* Initialize to zero */
    for (i=0; i<cst->nrows; i++) {
        cst->val[i][pos] = 0;
    }
}

/* Add a lower bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_lb(PlutoConstraints *cst, int varnum, int lb)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_inequality(cst, cst->nrows);

    cst->val[cst->nrows-1][varnum] = 1;
    cst->val[cst->nrows-1][cst->ncols-1] = -lb;
}

/* Add an upper bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_ub(PlutoConstraints *cst, int varnum, int ub)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_inequality(cst, cst->nrows);

    cst->val[cst->nrows-1][varnum] = -1;
    cst->val[cst->nrows-1][cst->ncols-1] = ub;
}



/* Set a value for a variable: varnum: 0-indexed */
void pluto_constraints_set_var(PlutoConstraints *cst, int varnum, int val)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_equality(cst, cst->nrows);

    cst->val[cst->nrows-1][varnum] = 1;
    cst->val[cst->nrows-1][cst->ncols-1] = -val;
}



/*
 * Returns the best candidate to eliminate (exact index in cst)
 * max_elim: maximum number of variables to eliminate (from the right)
 */
int best_elim_candidate(const PlutoConstraints *cst, int max_elim)
{
    int **csm, i, j, ub, lb, cost;

    int min_cost = cst->nrows*cst->nrows/4;
    int best_candidate = cst->ncols-2;

    csm = cst->val;

    for (j=cst->ncols-2; j > cst->ncols-2-max_elim; j--)    {
        ub=0;
        lb=0;
        for (i=0; i < cst->nrows; i++)    {
            if (csm[i][j] > 0) ub++;
            else if (csm[i][j] < 0) lb++;
        }
        /* cost = MIN(lb, ub); */
        cost = lb*ub;
        if (cost < min_cost)    {
            min_cost = cost;
            best_candidate = j;
        }
    }

    return best_candidate;
}


/* Populate a PIP matrix */
PipMatrix *pip_matrix_populate(int **cst, int nrows, int ncols)
{
    int i, j;
    PipMatrix *matrix ;
    Entier *p ;

    matrix = pip_matrix_alloc(nrows, ncols) ;

    p = matrix->p_Init ;
    for (i=0;i<matrix->NbRows;i++)  {
        for (j=0;j<matrix->NbColumns;j++)   {
#ifdef PIP_WIDTH_MP
            mpz_set_si(*(p++), cst[i][j]) ;
#else
            *(p++) = cst[i][j];
#endif
        }
    }
    return matrix ;
}

PlutoConstraints *pluto_constraints_select_row(const PlutoConstraints *cst, int pos)
{
    int j;

    PlutoConstraints *row = pluto_constraints_alloc(1, cst->ncols);
    row->is_eq[0] = cst->is_eq[pos];
    for (j=0; j<cst->ncols; j++)  {
        row->val[0][j] = cst->val[pos][j];
    }
    row->nrows = 1;
    return row;
}

/*
 * Negate a single row in cst.
 */
void pluto_constraints_negate_row(PlutoConstraints *cst, int pos)
{
    int j;

    for (j=0; j<cst->ncols; j++)    {
        cst->val[pos][j] = -cst->val[pos][j];
    }
}


/*
 * Negation of a single constraint in cst.
 */
void pluto_constraints_negate_constraint(PlutoConstraints *cst, int pos)
{
    int j;

    for (j=0; j<cst->ncols-1; j++)    {
        cst->val[pos][j] = -cst->val[pos][j];
    }
    cst->val[pos][cst->ncols-1]--;
}


/* Convert everything to >= 0 form */
PlutoConstraints *pluto_constraints_to_pure_inequalities(const PlutoConstraints *cst)
{
    int i, j;

    PlutoConstraints *ineq = pluto_constraints_dup(cst);

    for (i=0; i<ineq->nrows; i++) {
        ineq->is_eq[i] = 0;
    }

    /* Add a constraint to make sum of all equalities <= 0 */
    PlutoConstraints *neg_eq = pluto_constraints_alloc(1, cst->ncols);

    int num_eq = 0;
    for (i=0; i<cst->nrows; i++) {
        if (cst->is_eq[i])  {
            num_eq++;
            for (j=0; j<cst->ncols; j++)    {
                neg_eq->val[0][j] -= cst->val[i][j];
            }
        }
    }

    if (num_eq >= 1)    {
        neg_eq->nrows = 1;
        pluto_constraints_add(ineq, neg_eq);
    }
    pluto_constraints_free(neg_eq);

    return ineq;
}

void check_redundancy(PlutoConstraints *cst)
{
    int i;
    PlutoConstraints *check = pluto_constraints_alloc(cst->nrows, cst->ncols);
    int count;

    count = 0;

    for (i=0; i<cst->nrows; i++)    {
        pluto_constraints_copy(check, cst);
        PlutoConstraints *row = pluto_constraints_select_row(cst, i); 
        pluto_constraints_remove_row(check, i); 
        pluto_constraints_negate_constraint(row, 0);
        pluto_constraints_add(check, row);
        if (!pluto_constraints_solve(check))  {
            // printf("%dth constraint is redundant\n", i);
            count++;
        }else{
            // printf("%dth constraint is not redundant\n", i);
        }
        pluto_constraints_free(row);
    }
    IF_DEBUG(printf("%d constraints redundant\n", count););
}
