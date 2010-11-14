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
 * initialized all to zero and inequalities (>= 0)
 */
PlutoConstraints *pluto_constraints_alloc(int max_rows, int max_cols)
{
    PlutoConstraints *cst;

    cst = (PlutoConstraints *)malloc(sizeof(PlutoConstraints));

    int size = PLMAX(1,max_rows)*PLMAX(1,max_cols)*sizeof(int);

    cst->buf = (int *) malloc(size);
    bzero(cst->buf, size);

    if (cst->buf == NULL) {
        fprintf(stderr, "Not enough memory\n");
        exit(1);
    }

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
/* TODO: automatic non-destructive resizing */
PlutoConstraints *pluto_constraints_add(PlutoConstraints *cst1, const PlutoConstraints *cst2)
{
    assert(cst1->ncols == cst2->ncols);
    assert(cst1->nrows+cst2->nrows <= cst1->alloc_nrows);

    int i;
    for (i=0; i < cst2->nrows; i++) {
        memcpy(cst1->val[cst1->nrows+i], 
                cst2->val[i], cst1->ncols*sizeof(int));
    }

    memcpy(&cst1->is_eq[cst1->nrows], cst2->is_eq, cst2->nrows*sizeof(int));

    cst1->nrows = cst1->nrows + cst2->nrows;
    cst1->ncols = cst1->ncols;

    return cst1;
}


static int cols_compar=-1;

static int compar (const void *e1, const void *e2)
{
    int *row1 = *(int **) e1;
    int *row2 = *(int **) e2;

    /* cols_compar is a global variable that shd be set prior to calling qsort which
     * uses this function; can't pass it to this function since the argument
     * list is what qsort needs; the problem comes since the compar function
     * we have here works on arrays as opposed to scalar data and so size can only 
     * be passed externally */
    assert(cols_compar != -1);

    int i;

    for (i=0; i<cols_compar; i++)  {
        if (row1[i] != row2[i]) break;
    }

    if (i==cols_compar) return 0;
    else if (row1[i] < row2[i])    {
        return -1;
    }else{
        return 1;
    }

    /* memcmp is inefficient compared to what's above when compiled with gcc */
    /* return memcmp(row1, row2, cols_compar*sizeof(int)); */
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

    /* A race condition here when compiling with -fopenmp */
    cols_compar = cst->ncols;
    qsort(cst->val, cst->nrows, sizeof(int *), compar);

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
 * Eliminates the nth variable, where n has to be between 0 and cst->ncols-2;
 * Remember that the last column is for the constant. Caller has to make sure
 * that cst is large enough to hold the resulting constraints (it may or may
 * not grow). The implementation does not have a redundancy check; it just
 * eliminates duplicates after gcd normalization
 */
void fourier_motzkin_eliminate(PlutoConstraints *cst, int n)
{
    int i, j, k, l, p, q;
    int lb, ub, nb;
    int **csm;
    int *bound;

    for (i=0; i<cst->nrows; i++) {
        /* Does not support equalities: please extend if necessary */
        assert(cst->is_eq[i] != 1);
    }

    static PlutoConstraints *newcst=NULL;
    int **newcsm=NULL;

    // newcst = pluto_constraints_alloc(cst->nrows*cst->nrows/4, cst->ncols);

    assert(n >= 0);
    assert(n <= cst->ncols-2);
	// At least one variable
	assert(cst->ncols >= 2);

    csm = cst->val;

    for (i=0; i<cst->nrows; i++)    {
        if (csm[i][n] != 0) break;
    }

    if (i==cst->nrows)  {

        if (!newcst || newcst->alloc_nrows < cst->nrows || newcst->alloc_ncols < cst->ncols)    {
            /* Reallocate */
            if (newcst) pluto_constraints_free(newcst);
            newcst = pluto_constraints_alloc(cst->nrows, cst->ncols);
        }
        newcsm = newcst->val;

        for (j=0; j<cst->nrows; j++)    {
            q=0;
            for (k=0; k<cst->ncols; k++){
                if (k!=n)   {
                    newcsm[j][q] = csm[j][k];
                    q++;
                }
            }
        }
        newcst->nrows = cst->nrows;
        newcst->ncols = cst->ncols-1;
    }else{
        bound = (int *) malloc(cst->nrows*sizeof(int));
        // printf("bound size: %lf\n", ceil((cst->nrows*cst->nrows)/4));

        lb=0;
        ub=0;
        nb=0;
        /* Variable does appear */
        for (j=0; j<cst->nrows; j++)    {
            if (csm[j][n] == 0) {
                bound[j] = NB;
                nb++;
            }else if (csm[j][n] > 0) {
                // printf("accessing: %d\n", j);
                bound[j] = LB;
                lb++;
            }
            else{
                bound[j] = UB;
                ub++;
            }
        }
        if (!newcst || newcst->alloc_nrows < lb*ub+nb || newcst->alloc_ncols < cst->ncols){
            if (newcst) pluto_constraints_free(newcst);
            newcst = pluto_constraints_alloc(lb*ub+nb, cst->ncols);
        }
        newcsm = newcst->val;

        p=0;
        for (j=0; j<cst->nrows; j++)    {
            if (bound[j] == UB) {
                for (k=0; k<cst->nrows; k++)    {
                    if (bound[k] == LB) {
                        q = 0;
                        for(l=0; l < cst->ncols; l++)  {
                            if (l!=n)   {
                                newcsm[p][q] = 
                                    csm[j][l]*(lcm(csm[k][n], 
                                                -csm[j][n])/(-csm[j][n])) 
                                    + csm[k][l]*(lcm(-csm[j][n], 
                                                csm[k][n])/csm[k][n]); 
                                q++;
                            }
                        }
                        p++;
                    }
                }
            }else if (bound[j] == NB)   {
                q = 0;
                for (l=0; l<cst->ncols; l++)    {
                    if (l!=n)   {
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

    if (cst->alloc_nrows <= newcst->nrows)   {
        fprintf(stderr, "[FM eliminated] Not sufficient space allocated by caller\n");
        fprintf(stderr, "Aborting\n");
        exit(2);
    }

    pluto_constraints_copy(cst, newcst);
}


/* Copy constraints from src into dest; dest should have enough space */
PlutoConstraints *pluto_constraints_copy(PlutoConstraints *dest, const PlutoConstraints *src)
{
    int i;

    assert(src->nrows <= dest->alloc_nrows);
    assert(src->ncols <= dest->alloc_ncols);

    dest->nrows = src->nrows;
    dest->ncols = src->ncols;

    for (i=0; i<dest->nrows; i++) {
        memcpy(dest->val[i], src->val[i], src->ncols*sizeof(int));
    }

    memcpy(dest->is_eq, src->is_eq, dest->nrows*sizeof(int));

    return dest;
}


void pluto_constraints_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j;

    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols);

    for (i=0; i<cst->nrows; i++) {
        for (j=0; j<cst->ncols; j++) {
            fprintf(fp, "%s%d ", cst->val[i][j]>=0? " ":"", cst->val[i][j]);
        }
        fprintf(fp, "%s0\n", cst->is_eq[i]? " == ": ">= ");
    }
    fprintf(fp, "\n");
}


void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j;

    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols);

    for (i=0; i<cst->nrows; i++) {
        for (j=0; j<cst->ncols; j++) {
            fprintf(fp, "%s%d ", cst->val[i][j]>=0? " ":"", cst->val[i][j]);
        }
        fprintf(fp, "\t %s 0\n", cst->is_eq[i]? ">=": "==");
    }
    fprintf(fp, "\n");
}


/* Little more than */
void pluto_constraints_pretty_print2(FILE *fp, PlutoConstraints *cst)
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
        fprintf(fp, "%s 0\n", cst->is_eq[i]? ">=": "==");
    }
    fprintf(fp, "\n");
}


/* Convert Pluto constraints into PIP format (first column is
 * 0/1 based on equality/inequality */
PlutoMatrix *pluto2pip(const PlutoConstraints *cst, PlutoMatrix *pmat)
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


/* Use PIP to solve these constraints */
int *pluto_constraints_solve(const PlutoConstraints *cst)
{
    int bignum, i;
    PipMatrix  *domain, *context;
    PipQuast   *solution;
    PipOptions *pipOptions;
    PipList *listPtr;
    int *sol;
    static PlutoMatrix *pipmat = NULL;

    if (!pipmat || pipmat->alloc_nrows < cst->nrows || pipmat->alloc_ncols < cst->ncols+1)    {
        if (pipmat) pluto_matrix_free(pipmat);
        pipmat = pluto_matrix_alloc(cst->nrows, cst->ncols+1);
    }

    /* Convert constraints to PIP format */
    pluto2pip(cst, pipmat);

    /* First column says whether it's an inequality and the last is for the constant;
     * so we have ncols-2 variables */

    context = NULL;
    bignum = -1;

    domain = pip_matrix_populate(pipmat->val, pipmat->nrows, 
            pipmat->ncols);
    pipOptions = pip_options_init() ;

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

    return sol;
}


/* Dump in polylib format */
void pluto_constraints_dump_polylib(PlutoConstraints *cst)
{
    int i, j;

    printf("%d %d\n", cst->nrows, cst->ncols+1);

    for (i=0; i<cst->nrows; i++)    {
        printf("%s ", cst->is_eq[i]? "0": "1");
        for (j=0; j<cst->ncols; j++)    {
            printf("%d ", cst->val[i][j]);
        }
        printf("\n");
    }
}


/* Add an inequality (>= 0); initialize it to all zero */
void pluto_constraints_add_inequality(PlutoConstraints *cst, int pos)
{
    int i, j;

    assert (pos >= 0 && pos <= cst->nrows);
    assert (cst->nrows <= cst->alloc_nrows);

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

    assert (pos >= 0 && pos <= cst->nrows);
    assert (cst->nrows <= cst->alloc_nrows);

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



/* Remove a row */
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


/* Add a variable; resize if necessary */
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
    cst->val[cst->nrows-1][cst->ncols] = -lb;
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
