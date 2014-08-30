/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
#include <limits.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"

#include "piplib/piplib64.h"

#define UB 0
#define LB 1
#define NB 2

PipMatrix *pip_matrix_populate(int64 **cst, int nrows, int ncols);

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

    int size = PLMAX(1,max_rows)*PLMAX(1,max_cols)*sizeof(int64);

    cst->buf = (int64 *) malloc(size);

    if (cst->buf == NULL) {
        fprintf(stderr, "Not enough memory for allocating constraints\n");
        exit(1);
    }

    bzero(cst->buf, size);

    cst->is_eq = malloc(max_rows*sizeof(int));
    bzero(cst->is_eq, max_rows*sizeof(int));

    cst->val = (int64 **)malloc(max_rows*sizeof(int64 *));

    int i;
    for(i=0; i<max_rows; i++)   {
        cst->val[i] = &cst->buf[i*max_cols];
    }

    cst->alloc_nrows = max_rows;
    cst->alloc_ncols = max_cols;
    cst->names = NULL;
    cst->next = NULL;

    cst->nrows = 0;
    cst->ncols = max_cols;

    return cst;
}


void pluto_constraints_free(PlutoConstraints *cst)
{
    int i;

    if (cst == NULL) return;
    free(cst->buf);
    free(cst->val);
    free(cst->is_eq);
    if (cst->names) {
        for (i=0; i<cst->ncols-1; i++) {
            free(cst->names[i]);
        }
        free(cst->names);
    }

    if (cst->next != NULL)  pluto_constraints_free(cst->next);

    free(cst);
}


/* Non-destructive resize */
void pluto_constraints_resize(PlutoConstraints *cst, int nrows, int ncols)
{
    pluto_constraints_resize_single(cst, nrows, ncols);

    if (cst->next != NULL)  {
        pluto_constraints_resize(cst->next, nrows, ncols);
    }
}


/*
 * Non-destructive resize (single element of constraint list)
 * Not just change the allocated size, but also increase/reduce the number
 * of constraints/dimensions to match nrows and ncols
 */
void pluto_constraints_resize_single(PlutoConstraints *cst, int nrows, int ncols)
{
    int i, j;
    PlutoConstraints *newCst;

    assert(nrows >= 0 && ncols >= 0);

    newCst = pluto_constraints_alloc(PLMAX(nrows,cst->alloc_nrows), 
            PLMAX(ncols,cst->alloc_ncols));

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

    if (cst->names) {
        for (i=ncols-1; i<cst->ncols-1; i++) {
            free(cst->names[i]);
        }
    }

    if (cst->names) {
        cst->names = realloc(cst->names, (ncols-1)*sizeof(char *));
        for (i=cst->ncols-1; i<ncols-1; i++) {
            cst->names[i] = NULL;
        }
    }

    cst->nrows = nrows;
    cst->ncols = ncols;
    cst->alloc_nrows = newCst->alloc_nrows;
    cst->alloc_ncols = newCst->alloc_ncols;

    cst->val = newCst->val;
    cst->buf = newCst->buf;
    cst->is_eq = newCst->is_eq;

    free(newCst);
}


/* Adds cs1 and cs2 and puts them in cs1 -> returns cs1 itself; only for
 * single element list. Multiple elements doesn't make sense; you may want to
 * use intersect
 */
PlutoConstraints *pluto_constraints_add(PlutoConstraints *cst1, const PlutoConstraints *cst2)
{
    assert(cst2 != NULL);
    assert(cst1->ncols == cst2->ncols);
    assert(cst1->next == NULL);
    assert(cst2->next == NULL);

    if (cst1->nrows+cst2->nrows > cst1->alloc_nrows) {
        pluto_constraints_resize(cst1, cst1->nrows+cst2->nrows, cst1->ncols);
    }else{
        cst1->nrows = cst1->nrows + cst2->nrows;
    }

    int i;
    for (i=0; i < cst2->nrows; i++) {
        memcpy(cst1->val[cst1->nrows-cst2->nrows+i], 
                cst2->val[i], cst1->ncols*sizeof(int64));
    }

    memcpy(&cst1->is_eq[cst1->nrows-cst2->nrows], cst2->is_eq, cst2->nrows*sizeof(int));

    return cst1;
}


/* Adds cs2 to each element in cst1's list; cst2 should have a single element */
PlutoConstraints *pluto_constraints_add_to_each(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    int i;

    assert(cst2->next == NULL);
    assert(cst1->ncols == cst2->ncols);

    if (cst1->nrows+cst2->nrows > cst1->alloc_nrows) {
        pluto_constraints_resize(cst1, cst1->nrows+cst2->nrows, cst1->ncols);
    }else{
        cst1->nrows = cst1->nrows + cst2->nrows;
    }

    for (i=0; i < cst2->nrows; i++) {
        memcpy(cst1->val[cst1->nrows-cst2->nrows+i], 
                cst2->val[i], cst1->ncols*sizeof(int64));
    }

    memcpy(&cst1->is_eq[cst1->nrows-cst2->nrows], cst2->is_eq, cst2->nrows*sizeof(int));

    if (cst1->next != NULL) {
        pluto_constraints_add_to_each(cst1->next, cst2);
    }
    return cst1;
}

/* Temporary data structure to compare two rows */
struct row_info {
    int64 *row;
    int ncols;
};

static int row_compar(const void *e1, const void *e2)
{
    int i, ncols;
    int64 *row1, *row2;
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
    /* return memcmp(row1, row2, ncols*sizeof(int64)); */
}


/* 
 * Eliminates duplicate constraints; the simplified constraints
 * are still at the same memory location but the number of constraints 
 * in it will decrease
 */
void pluto_constraints_simplify(PlutoConstraints *const cst)
{
    int i, j, p;
    int64 _gcd;

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

        assert(_gcd != 0);
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
            memcpy(tmpcst->val[p], cst->val[i], cst->ncols*sizeof(int64));
            p++;
        }
    }
    tmpcst->nrows = p;
    tmpcst->ncols = cst->ncols;

    pluto_constraints_copy_single(cst, tmpcst);
    pluto_constraints_free(tmpcst);

    free(is_redun);

    if (cst->next != NULL) pluto_constraints_simplify(cst->next);
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
    int64 lb, ub, nb;
    int *bound;

    // At least one variable
    assert(cst->ncols >= 2);
    assert(pos >= 0);
    assert(pos <= cst->ncols-2);

    for (i=0; i<cst->nrows; i++) {
        if (cst->is_eq[i]) {
            PlutoConstraints *tmpcst = pluto_constraints_to_pure_inequalities_single(cst);
            pluto_constraints_copy_single(cst, tmpcst);
            pluto_constraints_free(tmpcst);
            break;
        }
    }

    PlutoConstraints *newcst;

    // newcst = pluto_constraints_alloc(cst->nrows*cst->nrows/4, cst->ncols);

    for (i=0; i<cst->nrows; i++)    {
        if (cst->val[i][pos] != 0) break;
    }

    if (i==cst->nrows) {
        newcst = pluto_constraints_dup(cst);
        pluto_constraints_remove_dim(newcst, pos);
    }else{
        bound = (int *) malloc(cst->nrows*sizeof(int));

        lb=0;
        ub=0;
        nb=0;
        /* Variable does appear */
        for (j=0; j<cst->nrows; j++)    {
            if (cst->val[j][pos] == 0) {
                bound[j] = NB;
                nb++;
            }else if (cst->val[j][pos] > 0) {
                bound[j] = LB;
                lb++;
            }else{
                bound[j] = UB;
                ub++;
            }
        }
        newcst = pluto_constraints_alloc(lb*ub+nb, cst->ncols-1);

        p=0;
        for (j=0; j<cst->nrows; j++)    {
            if (bound[j] == UB) {
                for (k=0; k<cst->nrows; k++)    {
                    if (bound[k] == LB) {
                        q = 0;
                        for(l=0; l < cst->ncols; l++)  {
                            if (l!=pos)   {
                                newcst->val[p][q] = 
                                    cst->val[j][l]*(lcm(cst->val[k][pos], 
                                                -cst->val[j][pos])/(-cst->val[j][pos])) 
                                    + cst->val[k][l]*(lcm(-cst->val[j][pos], 
                                                cst->val[k][pos])/cst->val[k][pos]); 
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
                        newcst->val[p][q] = cst->val[j][l];
                        q++;
                    }
                }
                p++;
            }
        }
        assert(p <= lb*ub + nb);
        newcst->nrows = p;
        newcst->ncols = cst->ncols-1;
        free(bound);
    }

    pluto_constraints_simplify(newcst);
    pluto_constraints_copy_single(cst, newcst);
    pluto_constraints_free(newcst);

    if (cst->next != NULL) fourier_motzkin_eliminate(cst->next,pos);
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
        memcpy(dest->val[i], src->val[i], src->ncols*sizeof(int64));
    }

    memcpy(dest->is_eq, src->is_eq, dest->nrows*sizeof(int));
    if (src->names) pluto_constraints_set_names(dest, src->names);


    if (src->next != NULL && dest->next != NULL)  {
        pluto_constraints_copy(dest->next, src->next);
    }
    if (src->next == NULL && dest->next != NULL)  {
        pluto_constraints_free(dest->next);
        dest->next = NULL;
    }
    if (src->next != NULL && dest->next == NULL)  {
        dest->next = pluto_constraints_dup(src->next);
    }

    return dest;
}

/* Copy constraints from the first element of src into the first element
 * of dest; if dest does not have enough space, resize it. The _single
 * signifies only this src in the list is copied */
PlutoConstraints *pluto_constraints_copy_single(PlutoConstraints *dest, const PlutoConstraints *src)
{
    int i;

    if (src->nrows > dest->alloc_nrows || src->ncols > dest->alloc_ncols) {
        pluto_constraints_resize_single(dest, PLMAX(src->nrows,dest->alloc_nrows),
                PLMAX(src->ncols,dest->alloc_ncols));
    }   

    dest->nrows = src->nrows;
    dest->ncols = src->ncols;

    for (i=0; i<dest->nrows; i++) {
        memcpy(dest->val[i], src->val[i], src->ncols*sizeof(int64));
    }

    memcpy(dest->is_eq, src->is_eq, dest->nrows*sizeof(int));
    if (src->names) {
        pluto_constraints_set_names(dest, src->names);
    }else{
        pluto_constraints_remove_names_single(dest);
    }

    return dest;
}

/* Duplicate constraints; returned constraints should be freed with
 * pluto_constraints_free */
PlutoConstraints *pluto_constraints_dup_single(const PlutoConstraints *src)
{
    assert(src != NULL);

    PlutoConstraints *dup = pluto_constraints_alloc(src->alloc_nrows, src->alloc_ncols);

    pluto_constraints_copy_single(dup, src);

    return dup;
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


static void pluto_constraints_print_single(FILE *fp, const PlutoConstraints *cst, int set_num)
{
    int i, j;

    if(cst == NULL)
        return;

    fprintf(fp, "Set #%d\n", set_num+1);
    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols);

    for (i=0; i<cst->nrows; i++) {
        for (j=0; j<cst->ncols; j++) {
            fprintf(fp, "% 3lld ", cst->val[i][j]);
        }
        fprintf(fp, "\t %s 0\n", cst->is_eq[i]? "==": ">=");
    }

    if (cst->next != NULL)  {
        pluto_constraints_print_single(fp, cst->next, set_num+1);
    }
}

void pluto_constraints_print(FILE *fp, const PlutoConstraints *cst)
{
    pluto_constraints_print_single(fp, cst, 0);
}


/* Print in polylib format */
static void pluto_constraints_print_polylib_without_num(FILE *fp,
        const PlutoConstraints *const cst)
{
    int i, j;

    fprintf(fp, "%d %d\n", cst->nrows, cst->ncols+1);

    for (i=0; i<cst->nrows; i++)    {
        fprintf(fp, "%s ", cst->is_eq[i]? "0": "1");
        for (j=0; j<cst->ncols; j++)    {
            fprintf(fp, "%lld ", cst->val[i][j]);
        }
        fprintf(fp, "\n");
    }

    if (cst->next != NULL) {
        pluto_constraints_print_polylib_without_num(fp, cst->next);
    }
}

void pluto_constraints_print_polylib(FILE *fp, const PlutoConstraints *const cst)
{
    int num;

    num =  pluto_constraints_num_in_list(cst);
    if (num >= 2) fprintf(fp, "%d\n", num);

    pluto_constraints_print_polylib_without_num(fp, cst);
}




/* Converts to polylib style matrix (first element of cst if it's a list) */
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


/* Create pluto_constraints from polylib-style matrix  */
PlutoConstraints *pluto_constraints_from_matrix(const PlutoMatrix *mat)
{
    int i, j;
    PlutoConstraints *cst;

    cst = pluto_constraints_alloc(mat->nrows, mat->ncols-1);

    cst->nrows = mat->nrows;

    for (i=0; i<cst->nrows; i++)    {
        cst->is_eq[i] = (mat->val[i][0] == 0);
        for (j=0; j<cst->ncols; j++)    {
            cst->val[i][j] = mat->val[i][j+1];
        }
    }

    return cst;
}

/* Create pluto_constraints from matrix  */
PlutoConstraints *pluto_constraints_from_inequalities(const PlutoMatrix *mat)
{
    int i, j;
    PlutoConstraints *cst;

    cst = pluto_constraints_alloc(mat->nrows, mat->ncols);
    cst->nrows = mat->nrows;

    for (i=0; i<cst->nrows; i++)    {
        cst->is_eq[i] = 0;
        for (j=0; j<cst->ncols; j++)    {
            cst->val[i][j] = mat->val[i][j];
        }
    }
    return cst;
}





/* Read constraints in polylib format from a file */
PlutoConstraints *pluto_constraints_read(FILE *fp)
{
    int i, j, ineq, nrows, ncols;

    fscanf(fp, "%d", &nrows);
    fscanf(fp, "%d", &ncols);
    ncols--;

    PlutoConstraints *cst = pluto_constraints_alloc(nrows, ncols);
    cst->nrows = nrows;

    for (i=0; i<cst->nrows; i++)    {
        fscanf(fp, "%d ", &ineq);
        cst->is_eq[i] = ineq^1;
        for (j=0; j<cst->ncols; j++)    {
            fscanf(fp, "%lld", &cst->val[i][j]);
        }
    }

    int num;
    if (fscanf(fp, "%d", &num) != EOF) {
        if (num >= 1) cst->names = malloc(num*sizeof(char *));
        for (i=0; i<num; i++) {
            cst->names[i] = malloc(6);
            fscanf(fp, "%s", cst->names[i]);
        }
    }

    return cst;
}


void pluto_constraints_compact_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j;

    int nrows = cst->nrows;
    int ncols = cst->ncols;

    assert(cst->next == NULL);

    if (nrows == 0) {
        printf("No constraints!\n");
    }

    for (i=0; i<nrows; i++) {
        int first = 1;
        for (j=0; j<ncols; j++)    {
            if (j==ncols-1) {
                /* constant */
                if (cst->val[i][j] == 0) fprintf(fp, " ");
                else fprintf(fp, "%s%lld ", 
                        (cst->val[i][j]>=0)?"+":"", cst->val[i][j]);
            }else{
                char var[6];
                var[5] = '\0';
                if (cst->names) strncpy(var, cst->names[j], 5);
                else sprintf(var, "c_%c", 'i'+j);

                if (cst->val[i][j] == 1) {
                    fprintf(fp, "%s%s", first?"":"+", var);
                    first = 0;
                }else if (cst->val[i][j] == -1) {
                    fprintf(fp, "-%s", var);
                    first = 0;
                }else if (cst->val[i][j] >= 2) {
                    fprintf(fp, "%s%lld%s", first?"":"+", cst->val[i][j], var);
                    first = 0;
                }else if (cst->val[i][j] <= -2) {
                    fprintf(fp, "%lld%s", cst->val[i][j], var);
                    first = 0;
                }
            }
        }
        fprintf(fp, "%s 0\n", cst->is_eq[i]? "==": ">=");
    }
    fprintf(fp, "\n");
}



void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst)
{
    int i, j;

    int nrows = cst->nrows;
    int ncols = cst->ncols;

    assert(cst->next == NULL);

    if (nrows == 0) {
        printf("No constraints!\n");
    }

    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++)    {
            if (j==ncols-1) {
                /* constant */
                if (cst->val[i][j] == 0) fprintf(fp, "     ");
                else fprintf(fp, "%s%lld ", (cst->val[i][j]>=0)?"+":"", cst->val[i][j]);
            }else{
                char var[6];
                var[5] = '\0';
                if (cst->names) strncpy(var, cst->names[j], 5);
                else sprintf(var, "c_%c", 'i'+j);

                if (cst->val[i][j] == 1) {
                    fprintf(fp, "+%s ", var);
                }else if (cst->val[i][j] == -1) {
                    fprintf(fp, "-%s ", var);
                }else if (cst->val[i][j] >= 2) {
                    fprintf(fp, "+%lld%s ", cst->val[i][j], var);
                }else if (cst->val[i][j] <= -2) {
                    fprintf(fp, "%lld%s ", cst->val[i][j], var);
                }else fprintf(fp, "     ");
            }
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

    //assert(cst->next == NULL);
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

/* Use PIP to solve these constraints (solves for the first element
 * if it's a list of constraints) */
int64 *pluto_constraints_solve_pip(const PlutoConstraints *cst, int negvar)
{
    int bignum, i;
    PipMatrix  *domain, *context;
    PipQuast   *solution;
    PipOptions *pipOptions;
    PipList *listPtr;
    int64 *sol;
    PlutoMatrix *pipmat;

    pipmat = pluto_matrix_alloc(cst->nrows, cst->ncols+1);

    /* Convert constraints to PIP format */
    pluto_constraints_to_pip_matrix(cst, pipmat);

    /* First column says whether it's an inequality and  
     * the last is for the constant; so we have ncols-2 variables */
    context = NULL;
    bignum = -1;

    domain = pip_matrix_populate(pipmat->val, pipmat->nrows, 
            pipmat->ncols);

    pipOptions = pip_options_init();

    if(negvar == 1) {
        pipOptions->Urs_parms=1;
        pipOptions->Urs_unknowns=1;
    }

    /* IF_DEBUG2(fprintf(stdout, "Calling PIP on a %dx%d formulation\n", 
       pipmat->nrows, pipmat->ncols)); */

    solution = pip_solve(domain, context, bignum, pipOptions) ;

    // IF_DEBUG2(pip_quast_print(stdout, solution, 0));

    assert(solution->condition==NULL);

    listPtr = solution->list;

    sol = NULL;
    if (listPtr != NULL)    {
        sol = (int64 *) malloc((pipmat->ncols-2)*sizeof(int64));

        for (i=0; i<pipmat->ncols-2; i++)    {
            /* This is just a lexmin and not a parametrix lexmin and so each
             * vector in the list is actually a constant */
#ifdef PIP_WIDTH_MP
            sol[i] = mpz_get_si(*listPtr->vector->the_vector);
#else
            sol[i] = (int64) *listPtr->vector->the_vector;
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
int64 *pluto_constraints_solve(const PlutoConstraints *cst, int negvar)
{
    if (options->islsolve) {
        return pluto_constraints_solve_isl(cst, negvar);
    }else{
        return pluto_constraints_solve_pip(cst, negvar);
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

/* Add an inequality (>= 0); initialize it to all zero; will be added
 * as the last row */
void pluto_constraints_add_inequality(PlutoConstraints *cst)
{
    int j;

    if (cst->nrows == cst->alloc_nrows)   {
        pluto_constraints_resize_single(cst, cst->nrows+1, cst->ncols);
    }else{
        cst->nrows++;
    }

    for (j=0; j<cst->ncols; j++) {
        cst->val[cst->nrows-1][j] = 0;
    }
    cst->is_eq[cst->nrows-1] = 0;

    if (cst->next != NULL) {
        pluto_constraints_add_inequality(cst->next);
    }

}


/* Add an equality (== 0); initialize it to all zero; will be added
 * as the last row */
void pluto_constraints_add_equality(PlutoConstraints *cst)
{
    int j;

    assert(cst->nrows <= cst->alloc_nrows);

    if (cst->nrows == cst->alloc_nrows)   {
        pluto_constraints_resize_single(cst, cst->nrows+1, cst->ncols);
    }else{
        cst->nrows++;
    }

    for (j=0; j<cst->ncols; j++) {
        cst->val[cst->nrows-1][j] = 0;
    }
    cst->is_eq[cst->nrows-1] = 1;

    if (cst->next != NULL) {
        pluto_constraints_add_equality(cst->next);
    }
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

    if (cst->names) free(cst->names[pos]);
    for (j=pos; j<cst->ncols-1; j++) {
        for (i=0; i<cst->nrows; i++) {
            cst->val[i][j] = cst->val[i][j+1];
        }
        if (j < cst->ncols-2 && cst->names[j+1]) {
            cst->names[j] = cst->names[j+1];
        }
    }
    if (cst->names) {
        cst->names = realloc(cst->names, (cst->ncols-2)*sizeof(char *));
    }

    cst->ncols--;

    if (cst->next != NULL) pluto_constraints_remove_dim(cst->next, pos);
}


/* Blank 'pos' th row */
void pluto_constraints_zero_row(PlutoConstraints *cst, int pos)
{
    int j;

    assert(cst->next == NULL);

    assert(pos >= 0 && pos <= cst->nrows-1);

    for (j=0; j<cst->ncols; j++)    {
        cst->val[pos][j] = 0;
    }
}


/* Normalize row by its gcd */
void pluto_constraints_normalize_row(PlutoConstraints *cst, int pos)
{
    int i, j, k;

    assert(cst->next == NULL);
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
void pluto_constraints_add_dim(PlutoConstraints *cst, int pos, const char *name)
{
    int i, j;

    /* Has to be a new variable */
    assert(pos >= 0 && pos <= cst->ncols-1);

    if (cst->ncols == cst->alloc_ncols)  {
        pluto_constraints_resize_single(cst, cst->nrows, cst->ncols+1);
    }else{
        cst->ncols++;
    }
    if (cst->names) cst->names = realloc(cst->names, (cst->ncols-1)*sizeof(char *));

    for (j=cst->ncols-2; j>=pos; j--) {
        for (i=0; i<cst->nrows; i++) {
            cst->val[i][j+1] = cst->val[i][j];
        }
        if (cst->names && j <= cst->ncols-3) {
            cst->names[j+1] = cst->names[j];
        }
    }

    /* Initialize to zero */
    for (i=0; i<cst->nrows; i++) {
        cst->val[i][pos] = 0;
    }
    if (cst->names) cst->names[pos] = strdup(name? name:"_u");

    if (cst->next != NULL)  {
        pluto_constraints_add_dim(cst->next, pos, name);
    }
}

/* Add a lower bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_lb(PlutoConstraints *cst, int varnum, int64 lb)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_inequality(cst);

    while (cst != NULL) {
        cst->val[cst->nrows-1][varnum] = 1;
        cst->val[cst->nrows-1][cst->ncols-1] = -lb;
        cst = cst->next;
    }
}

/* Add an upper bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_ub(PlutoConstraints *cst, int varnum, int64 ub)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_inequality(cst);

    while (cst != NULL) {
        cst->val[cst->nrows-1][varnum] = -1;
        cst->val[cst->nrows-1][cst->ncols-1] = ub;
        cst = cst->next;
    }
}



/* Set a value for a variable: varnum: 0-indexed */
void pluto_constraints_set_var(PlutoConstraints *cst, int varnum, int64 val)
{
    assert(varnum >=0 && varnum <= cst->ncols-2);

    pluto_constraints_add_equality(cst);

    while (cst != NULL)  {
        cst->val[cst->nrows-1][varnum] = 1;
        cst->val[cst->nrows-1][cst->ncols-1] = -val;
        cst = cst->next;
    }
}




/* Populate a PIP matrix */
PipMatrix *pip_matrix_populate(int64 **cst, int nrows, int ncols)
{
    int i, j;
    PipMatrix *matrix ;
    Entier *p ;

    matrix = pip_matrix_alloc(nrows, ncols) ;

    p = matrix->p_Init ;
    for (i=0;i<matrix->NbRows;i++)  {
        for (j=0;j<matrix->NbColumns;j++)   {
            *(p++) = cst[i][j];
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
    assert(cst->next == NULL);

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

    assert(cst->next == NULL);

    for (j=0; j<cst->ncols-1; j++)    {
        cst->val[pos][j] = -cst->val[pos][j];
    }
    cst->val[pos][cst->ncols-1]--;
}


/* Convert everything to >= 0 form; for all equalities, add a single extra 
 * constraint that puts the sum of those to be <= 0; it is only performed on a
 * single PlutoConstraints in a list */
PlutoConstraints *pluto_constraints_to_pure_inequalities_single(const PlutoConstraints *cst)
{
    int i, j;

    PlutoConstraints *ineq = pluto_constraints_dup_single(cst);

    for (i=0; i<ineq->nrows; i++) {
        ineq->is_eq[i] = 0;
    }

    int has_eq = 0;
    for (i=0; i<cst->nrows; i++) {
        if (cst->is_eq[i])  {
            has_eq = 1;
            break;
        }
    }

    if (has_eq) {
        /* Add a constraint to make sum of all equalities <= 0 */
        PlutoConstraints *neg_eq = pluto_constraints_alloc(1, cst->ncols);

        for (i=0; i<cst->nrows; i++) {
            if (cst->is_eq[i])  {
                for (j=0; j<cst->ncols; j++)    {
                    neg_eq->val[0][j] -= cst->val[i][j];
                }
            }
        }
        neg_eq->nrows = 1;
        pluto_constraints_add(ineq, neg_eq);
        pluto_constraints_free(neg_eq);
    }

    return ineq;
}


/* start: 0-indexed */
void pluto_constraints_project_out_single(
        PlutoConstraints *cst, 
        int start, 
        int num)
{
    int i, end;

    assert(num>= 0);

    if (num == 0)   return;

    end = start + num - 1;

    assert(start >= 0 && end <= cst->ncols-2);

    PlutoMatrix *func = pluto_matrix_alloc(cst->ncols-num, cst->ncols);
    pluto_matrix_initialize(func, 0);
    for (i=0; i<start; i++) {
        func->val[i][i] = 1;
    }
    for (i=end+1; i<cst->ncols; i++) {
        func->val[i-num][i] = 1;
    }

    PlutoConstraints *img =  pluto_constraints_image(cst, func);
    pluto_constraints_copy_single(cst, img);

    pluto_constraints_free(img);
    pluto_matrix_free(func);
}

/* start: 0-indexed */
void pluto_constraints_project_out(
        PlutoConstraints *cst, 
        int start, 
        int num)
{
    pluto_constraints_project_out_single(cst, start, num);

    if (cst->next != NULL)  {
        pluto_constraints_project_out(cst->next, start, num);
    }
}


void pluto_constraints_interchange_cols(PlutoConstraints *cst, int col1,
        int col2)
{
    int r, tmp;

    //assert(cst->next == NULL);

    for (r=0; r<cst->nrows; r++)    {
        tmp = cst->val[r][col1];
        cst->val[r][col1] = cst->val[r][col2];
        cst->val[r][col2] = tmp;
    }

    if(cst->next != NULL) {
        pluto_constraints_interchange_cols(cst->next, col1, col2);
    }

}


void check_redundancy(PlutoConstraints *cst)
{
    int i;
    PlutoConstraints *check = pluto_constraints_alloc(cst->nrows, cst->ncols);
    int count;

    assert(cst->next == NULL);

    count = 0;

    for (i=0; i<cst->nrows; i++)    {
        pluto_constraints_copy(check, cst);
        PlutoConstraints *row = pluto_constraints_select_row(cst, i); 
        pluto_constraints_remove_row(check, i); 
        pluto_constraints_negate_constraint(row, 0);
        pluto_constraints_add(check, row);
        if (!pluto_constraints_solve(check,DO_NOT_ALLOW_NEGATIVE_COEFF))  {
            // printf("%dth constraint is redundant\n", i);
            count++;
        }else{
            // printf("%dth constraint is not redundant\n", i);
        }
        pluto_constraints_free(row);
    }
    IF_DEBUG(printf("%d constraints redundant\n", count););
}

int pluto_constraints_num_in_list(const PlutoConstraints *const cst)
{
    if (cst == NULL)  return 0;

    return 1 + pluto_constraints_num_in_list(cst->next);
}

/* In-place union: first argument is modified */
PlutoConstraints *pluto_constraints_unionize(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    PlutoConstraints *ucst = pluto_constraints_union(cst1, cst2);
    pluto_constraints_copy(cst1, ucst);
    pluto_constraints_free(ucst);
    return cst1;
}

/* In-place union: first argument is modified */
PlutoConstraints *pluto_constraints_unionize_isl(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    PlutoConstraints *ucst = pluto_constraints_union_isl(cst1, cst2);
    pluto_constraints_copy(cst1, ucst);
    pluto_constraints_free(ucst);
    return cst1;
}



/* In-place intersection: first argument is modified */
PlutoConstraints *pluto_constraints_intersect(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    PlutoConstraints *icst = pluto_constraints_intersection(cst1, cst2);
    pluto_constraints_copy(cst1, icst);
    pluto_constraints_free(icst);

    return cst1;
}


PlutoConstraints *pluto_constraints_universe(int ncols)
{
    PlutoConstraints *universe = pluto_constraints_alloc(1, ncols);
    return universe;
}


/* Return an empty polyhedron */
PlutoConstraints *pluto_constraints_empty(int ncols)
{
    PlutoConstraints *empty = pluto_constraints_alloc(1, ncols);
    pluto_constraints_add_equality(empty);
    empty->val[0][ncols-1] = 1;
    return empty;
}


int pluto_constraints_is_empty(const PlutoConstraints *cst)
{
    int64 *sol;
    bool is_empty;

    sol = pluto_constraints_solve(cst, DO_NOT_ALLOW_NEGATIVE_COEFF);

    is_empty = (sol == NULL);
    if (sol != NULL) free(sol);

    if (!is_empty) return 0;

    if (cst->next != NULL)  return pluto_constraints_is_empty(cst->next);

    /* This one is empty and there are no more */
    return 1;
}

/* Append the constraints together */
PlutoConstraints *pluto_constraints_unionize_simple(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    while (cst1->next != NULL) {
        cst1 = cst1->next;
    }
    cst1->next = pluto_constraints_dup(cst2);
    return cst1;
}


/* Get lower bound for pos^th variable if it's a (single) constant: return 1 
 * if constant, 0 otherwise */
int pluto_constraints_get_const_lb(const PlutoConstraints *cnst, int pos, 
        int64 *lb)
{
    int i, retval;

    PlutoConstraints *cst = pluto_constraints_dup(cnst);

    pluto_constraints_project_out_single(cst, 0, pos);
    pluto_constraints_project_out_single(cst, 1, cst->ncols-2);

    retval = 0;
    *lb = INT_MIN;
    for (i=0; i<cst->nrows; i++) {
        if (cst->is_eq[i]) {
            *lb = cst->val[i][cst->ncols-1];
            retval = 1;
            break;
        }
        if (!cst->is_eq[i] && cst->val[i][0] >= 1) {
            retval = 1;
            *lb = PLMAX(*lb,-cst->val[i][cst->ncols-1]);
        }
    }
    pluto_constraints_free(cst);

    if (retval && cnst->next != NULL) {
        int64 next_lb;
        retval = pluto_constraints_get_const_lb(cnst->next, pos, &next_lb);
        if (*lb != next_lb)  retval = 0;
    }

    return retval;
}


/* Get upper bound for the pos^th variable if it's a (single) constant: return 
 * 0 if not constant,  1 otherwise */
int pluto_constraints_get_const_ub(const PlutoConstraints *cnst, int pos, 
        int64 *ub)
{
    int i, retval;

    PlutoConstraints *cst = pluto_constraints_dup(cnst);

    pluto_constraints_project_out_single(cst, 0, pos);
    pluto_constraints_project_out_single(cst, 1, cst->ncols-2);
    //pluto_constraints_project_out_single_isl(cst, 0, pos);
    //pluto_constraints_project_out_single_isl(cst, 1, cst->ncols-2);

    retval = 0;
    *ub = LONG_LONG_INT_MAX;
    for (i=0; i<cst->nrows; i++) {
        if (cst->is_eq[i]) {
            *ub = cst->val[i][cst->ncols-1];
            retval = 1;
            break;
        }
        if (!cst->is_eq[i] && cst->val[i][0] <= -1) {
            *ub = PLMIN(*ub,cst->val[i][cst->ncols-1]);
            retval = 1;
        }
    }
    pluto_constraints_free(cst);

    if (retval && cnst->next != NULL) {
        int64 next_ub;
        retval = pluto_constraints_get_const_ub(cnst->next, pos, &next_ub);
        if (*ub != next_ub)  retval = 0;
    }

    return retval;
}

void print_polylib_visual_sets_internal(char* str, int k,  PlutoConstraints *cst){
    int i, j, first = 0;
    char name[100];

	if(!options->moredebug)
		return;

	if(cst == NULL)
		return;

	sprintf(name, "%si%d", str, k);

    printf("%s := { ", name );
    for( i=0; i<cst->ncols-1; i++){
        if(i == cst->ncols-2)
            printf("t%d | ", i);
        else
            printf("t%d,", i);
    }

    for(i=0; i<cst->nrows; i++ ){
        first = 0;
        for(j=0; j<cst->ncols-1; j++){
            if(cst->val[i][j] == 0)
                continue;

            if(first == 0){
                first = 1;
                if(cst->val[i][j] == 1)
                    printf("t%d", j);
                else if(cst->val[i][j] == -1)
                    printf("-t%d", j);
                else if(cst->val[i][j] != 0)
                    printf("%lld*t%d", cst->val[i][j], j);
            }
            else {
                first = 1;
                if(cst->val[i][j] == 1)
                    printf("+t%d", j);
                else if(cst->val[i][j] == -1)
                    printf("-t%d", j);
                else if(cst->val[i][j] != 0)
                    printf("%+lld*t%d", cst->val[i][j], j);
            }
        }

        if(cst->val[i][j] != 0){
            if(first != 0)
                printf("%+lld", cst->val[i][j]);
            else
                printf("%lld", cst->val[i][j]);
        }

        printf("%s0 ", cst->is_eq[i]?"=":">=" );
        if(i != cst->nrows -1)
            printf(",");
    }

    printf("};\n\n");

    if(cst->next != NULL)
        print_polylib_visual_sets_internal(str,k+1,  cst->next);
}



void print_polylib_visual_sets_internal_new(char* str, int k,  PlutoConstraints *cst) {

    int i, j, first = 0;
	char name[100];

	sprintf(name, "%s%d", str, k);

    printf("%s := { ", name );
    for( i=0; i<cst->ncols-1; i++){
        if(i == cst->ncols-2)
            printf("t%d | ", i);
        else
            printf("t%d,", i);
    }

    for(i=0; i<cst->nrows; i++ ){
        first = 0;
        for(j=0; j<cst->ncols-1; j++){
            if(cst->val[i][j] == 0)
                continue;

            if(first == 0){
                first = 1;
                if(cst->val[i][j] == 1)
                    printf("t%d", j);
                else if(cst->val[i][j] == -1)
                    printf("-t%d", j);
                else if(cst->val[i][j] != 0)
                    printf("%lld*t%d", cst->val[i][j], j);
            }
            else {
                first = 1;
                if(cst->val[i][j] == 1)
                    printf("+t%d", j);
                else if(cst->val[i][j] == -1)
                    printf("-t%d", j);
                else if(cst->val[i][j] != 0)
                    printf("%+lld*t%d", cst->val[i][j], j);
            }
        }

        if(cst->val[i][j] != 0){
            if(first != 0)
                printf("%+lld", cst->val[i][j]);
            else
                printf("%lld", cst->val[i][j]);
        }

        printf("%s0 ", cst->is_eq[i]?"=":">=" );
        if(i != cst->nrows -1)
            printf(",");
    }

    printf("};\n\n");

    if(cst->next != NULL)
        print_polylib_visual_sets_internal(str,k+1,  cst->next);
}

void print_polylib_visual_sets_new(char* name, PlutoConstraints *cst){
	print_polylib_visual_sets_internal_new(name, 0, cst);
	return;
}

void print_polylib_visual_sets(char* name, PlutoConstraints *cst)
{
	print_polylib_visual_sets_internal(name, 0, cst);
	return;
}

PlutoDepList* pluto_dep_list_alloc(Dep *dep){
    PlutoDepList *list = (PlutoDepList *)malloc(sizeof(PlutoDepList));
    assert(dep!= NULL);
    //list->dep = pluto_dependence_dup(dep, dep->dpolytope);
    list->dep = dep;
    list->next = NULL;
    return list;
}

void pluto_deps_list_free(PlutoDepList *deplist) {
    if (deplist == NULL)
        return;
    pluto_deps_list_free(deplist->next);
}

PlutoDepList* pluto_deps_list_dup(PlutoDepList *src){

    assert(src != NULL);

    PlutoDepList *new = pluto_dep_list_alloc(src->dep);

    if(src->next != NULL) {
        new->next = pluto_deps_list_dup(src->next);
    }

    return new;
}


void pluto_deps_list_append(PlutoDepList *list, Dep *dep) {

    assert(list != NULL);
    assert(dep != NULL);

    PlutoDepList *new = pluto_dep_list_alloc(dep);
    new->next = list->next;
    list->next = new;

    return;
}

PlutoConstraintsList* pluto_constraints_list_alloc(PlutoConstraints *cst) {

    PlutoConstraintsList *list = (PlutoConstraintsList *)malloc(sizeof(PlutoConstraintsList));
    list->constraints = cst;
    list->deps = NULL;
    list->next = NULL;
    return list;
}

void pluto_constraints_list_free(PlutoConstraintsList *cstlist) {
    if (cstlist == NULL)
        return;
    pluto_constraints_list_free(cstlist->next);
    pluto_constraints_free(cstlist->constraints);
    pluto_deps_list_free(cstlist->deps);
}

/*adds a new element to constraints list
*/
void pluto_constraints_list_add(PlutoConstraintsList *list,const PlutoConstraints *cst, Dep *dep, int copyDep) {

    assert(list != NULL);
    assert(cst != NULL);

    PlutoConstraintsList *new = pluto_constraints_list_alloc(pluto_constraints_dup(cst));
    new->next = list->next;
    list->next = new;

    if(copyDep) {
        new->deps = pluto_deps_list_dup(list->deps);
        pluto_deps_list_append(new->deps,dep);
    }
    else {
        new->deps = pluto_dep_list_alloc(dep);
    }
    return;
}

void pluto_constraints_list_replace(PlutoConstraintsList *list, PlutoConstraints *cst) {

    assert(list != NULL);
    assert(cst != NULL);

    pluto_constraints_free(list->constraints);
    list->constraints = pluto_constraints_dup(cst);

    return;
}

<<<<<<< HEAD
int pluto_constraints_get_num_equalities(const PlutoConstraints *cst)
{
    int r, count;

    count = 0;
    for (r=0; r<cst->nrows; r++) {
        if (cst->is_eq[r]) count++;
    }
    return count;
}


PlutoMatrix *pluto_constraints_extract_equalities(const PlutoConstraints *cst)
{
    int r, c, count;

    PlutoMatrix *mat = pluto_matrix_alloc(pluto_constraints_get_num_equalities(cst),
            cst->ncols);

    count=0;
    for (r=0; r<cst->nrows; r++) {
        if (cst->is_eq[r]) {
            for (c=0; c<cst->ncols; c++) {
                mat->val[count][c] = cst->val[r][c]; 
            }
            count++;
        }
    }

    return mat;
}


=======

void pluto_constraints_remove_names_single(PlutoConstraints *cst)
{
    int i;

    if (cst->names) {
        for (i=0; i<cst->ncols-1; i++) {
            free(cst->names[i]);
        }
    }
    free(cst->names);
    cst->names = NULL;

    // if (cst->next != NULL) pluto_constraints_remove_names(cst->next);
}


/*
 * Existing names if any are free
 */
>>>>>>> origin/master
void pluto_constraints_set_names(PlutoConstraints *cst, char **names)
{
    int i;

    if (cst->names) {
        for (i=0; i<cst->ncols-1; i++) {
            free(cst->names[i]);
        }
    }

    if (!cst->names) {
        cst->names = malloc((cst->ncols-1)*sizeof(char *));
    }

    for (i=0; i<cst->ncols-1; i++) {
        cst->names[i] = strdup(names[i]);
    }
}

void pluto_constraints_set_names_range(PlutoConstraints *cst, char **names,
        int dest_offset, int src_offset, int num)
{
    int i;

    if (num == 0) return;

    assert(names != NULL);
    assert(dest_offset+num <= cst->ncols-1);

    if (!cst->names) {
        cst->names = malloc((cst->ncols-1)*sizeof(char *));
        for (i=0; i<cst->ncols-1; i++) {
            cst->names[i] = NULL;
        }
    }else{
        for (i=0; i<num; i++) {
            free(cst->names[dest_offset+i]);
        }
    }

    for (i=0; i<num; i++) {
        cst->names[dest_offset+i] = strdup(names[src_offset+i]);
    }
}


/*
 * Returns the best candidate to eliminate (exact index in cst)
 * max_elim: maximum number of variables to eliminate (from the right)
 *
 * FIXME: update it for constraints with equalities
 */
int pluto_constraints_best_elim_candidate(const PlutoConstraints *cst, int max_elim)
{
    int64 **csm, i, j, ub, lb, cost;

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
