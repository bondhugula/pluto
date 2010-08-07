/*
 * PLUTO: An automatic parallelizer and locality optimizer
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

#define UB 0
#define LB 1
#define NB 2

#ifdef PIP_WIDTH_MP
#include "piplib/piplibMP.h"
#else
#include "piplib/piplib64.h"
#endif


PipMatrix *pip_matrix_populate(int **cst, int nrows, int ncols);

/*
 * Allocate with a max size of max_rows and max_cols
 * initialize them to zero
 */
PlutoInequalities *constraints_alloc(int max_rows, int max_cols)
{
    PlutoInequalities *cst;

    cst = (PlutoInequalities *)malloc(sizeof(PlutoInequalities));
    cst->val = (int **)malloc(max_rows*sizeof(int *));

    if (!cst)   {
        fprintf(stderr, "Not enough memory\n");
        exit(1);
    }

    int i;
    for(i=0; i<max_rows; i++)   {
        cst->val[i] = (int *)malloc(max_cols*sizeof(int));
        if (cst->val[i] == NULL)   {
            fprintf(stderr, "Not enough memory\n");
            exit(1);
        }
        bzero(cst->val[i], max_cols*sizeof(int));
    }

    cst->alloc_nrows = max_rows;
    cst->alloc_ncols = max_cols;

    cst->nrows = 0;
    cst->ncols = max_cols;

    return cst;
}


void constraints_free(PlutoInequalities *cst)
{
    pluto_matrix_free(cst);
}


/* Adds cs1 and cs2 and puts them in cs1 -> returns cs1 itself */
/* TODO: automatic non-destructive resizing */
PlutoInequalities *constraints_add(PlutoInequalities *cst1, PlutoInequalities *cst2)
{
    assert(cst1->ncols == cst2->ncols);
    assert(cst1->nrows+cst2->nrows <= cst1->alloc_nrows);

    int i;
    for (i=0; i < cst2->nrows; i++) {
        memcpy(cst1->val[cst1->nrows+i], 
                cst2->val[i], cst1->ncols*sizeof(int));
    }

    cst1->nrows = cst1->nrows + cst2->nrows;
    cst1->ncols = cst1->ncols;

    return cst1;
}


/*
 * Returns the best candidate to eliminate (exact index in cst)
 *
 * max_elim: maximum number of variables to eliminate (from the right)
 */
int best_elim_candidate(PlutoInequalities *cst, int max_elim)
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
void constraints_simplify(PlutoInequalities *const cst)
{
    if (cst->nrows == 0)    {
        return;
    }

    int i, j, p, _gcd;

    PlutoInequalities *tmpcst = constraints_alloc(cst->nrows, cst->ncols);

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

    constraints_copy(tmpcst, cst);
    constraints_free(tmpcst);

    free(is_redun);
}



/* 
 * Eliminates the nth variable, where n has to be between 0 and cst->ncols-2;
 * Remember that the last column is for the constant. Caller has to make sure
 * that cst is large enough to hold the resulting constraints (it may or may
 * not grow). The implementation does not have a redundancy check; it just
 * eliminates duplicates after gcd normalization
 */
void fourier_motzkin_eliminate(PlutoInequalities *cst, int n)
{
    int i, j, k, l, p, q;
    int lb, ub, nb;
    int **csm;
    int *bound;

    static PlutoInequalities *newcst=NULL;
    int **newcsm=NULL;

    // newcst = constraints_alloc(cst->nrows*cst->nrows/4, cst->ncols);

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
            if (newcst) constraints_free(newcst);
            newcst = constraints_alloc(cst->nrows, cst->ncols);
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
            if (newcst) constraints_free(newcst);
            newcst = constraints_alloc(lb*ub+nb, cst->ncols);
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

    constraints_simplify(newcst);

    if (cst->alloc_nrows <= newcst->nrows)   {
        fprintf(stderr, "[FM eliminated] Not sufficient space allocated by caller\n");
        fprintf(stderr, "Aborting\n");
        exit(2);
    }

    constraints_copy(newcst, cst);
}



PlutoInequalities *constraints_copy(PlutoInequalities *src, PlutoInequalities *dest)
{
    int i;

    assert(src->nrows <= dest->alloc_nrows);
    assert(src->ncols <= dest->alloc_ncols);

    dest->nrows = src->nrows;
    dest->ncols = src->ncols;

    for (i=0; i<dest->nrows; i++)
        memcpy(dest->val[i], src->val[i], src->ncols*sizeof(int));

    return dest;
}

void constraints_print(FILE *fp, PlutoInequalities *cst)
{
    pluto_matrix_print(fp, cst);
}

void constraints_add_lower_bound(int varnum)
{

}

void constraints_pretty_print(FILE *fp, PlutoInequalities *cst)
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
        fprintf(fp, ">= 0\n");
    }
    fprintf(fp, "\n");
}


/* Convert constraints from Pluto format to PIP format */
PlutoInequalities *pluto2pip(PlutoInequalities *cst, PlutoInequalities *pipcst)
{
    int i, j;

    assert(pipcst != NULL);

    for (i=0; i<cst->nrows; i++)    {
        /* Inequality */
        pipcst->val[i][0] = 1;
        for (j=1; j< cst->ncols+1; j++)    {
            pipcst->val[i][j] = cst->val[i][j-1];
        }
    }
    pipcst->nrows = cst->nrows;
    pipcst->ncols = cst->ncols+1;

    return pipcst;
}



/*
 * Use PIP to solve these constraints
 */
int *constraints_solve (PlutoInequalities *cst)
{
    int bignum, i;
    PipMatrix  *domain, *context;
    PipQuast   *solution;
    PipOptions *pipOptions;
    PipList *listPtr;
    int *sol;
    static PlutoInequalities *pipcst = NULL;

    if (!pipcst || pipcst->alloc_nrows < cst->nrows || pipcst->alloc_ncols < cst->ncols+1)    {
        if (pipcst) constraints_free(pipcst);
        pipcst = constraints_alloc(cst->nrows, cst->ncols+1);
    }

    /* Convert constraints to PIP format */
    pluto2pip(cst, pipcst);

    /* First column says whether it's an inequality and the last is for the constant;
     * so we have ncols-2 variables */

    context = NULL;
    bignum = -1;

    domain = pip_matrix_populate(pipcst->val, pipcst->nrows, 
            pipcst->ncols);
    pipOptions = pip_options_init() ;

    /* IF_DEBUG2(fprintf(stdout, "Calling PIP on a %dx%d formulation\n", 
                pipcst->nrows, pipcst->ncols)); */

    solution = pip_solve(domain,context,bignum,pipOptions) ;

    // IF_DEBUG2(pip_quast_print(stdout, solution, 0));

    assert(solution->condition==NULL);

    listPtr = solution->list;

    sol = NULL;
    if (listPtr != NULL)    {
        sol = (int *) malloc((pipcst->ncols-2)*sizeof(int));

        for (i=0; i<pipcst->ncols-2; i++)    {
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

void dump_polylib(PlutoInequalities *cst)
{
    int i, j;

    printf("%d %d\n", cst->nrows, cst->ncols+1);

    for (i=0; i<cst->nrows; i++)    {
        printf("1 ");
        for (j=0; j<cst->ncols; j++)    {
            printf("%d ", cst->val[i][j]);
        }
        printf("\n");
    }

}
