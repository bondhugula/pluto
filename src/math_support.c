/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2008 Uday Kumar Bondhugula
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

#include "math_support.h"
#include "constraints.h"

/*
 * Allocated; not initialized
 */
PlutoMatrix *pluto_matrix_alloc(int alloc_nrows, int alloc_ncols)
{
    int i;
    PlutoMatrix *mat;

    mat = (PlutoMatrix *) malloc(sizeof(PlutoMatrix));

    mat->val = (int **) malloc(alloc_nrows*sizeof(int *));

    for (i=0; i<alloc_nrows; i++) {
        mat->val[i] = (int *) malloc(alloc_ncols*sizeof(int));
    }

    mat->alloc_nrows = alloc_nrows;
    mat->alloc_ncols = alloc_ncols;

    mat->nrows = alloc_nrows;
    mat->ncols = alloc_ncols;

    return mat;
}

void pluto_matrix_free(PlutoMatrix *mat)
{
    int i;

    for (i=0; i<mat->alloc_nrows; i++) {
        free(mat->val[i]);
    }

    free(mat->val);
    free(mat);
}


/* Remove column <pos> from the matrix: pos starts from 0 */
void pluto_matrix_remove_col(PlutoMatrix *mat, int pos)
{
    int i, j;

    assert(pos <= mat->ncols-1);

    for (j=pos; j<mat->ncols-1; j++) {
        for (i=0; i<mat->nrows; i++) {
            mat->val[i][j] = mat->val[i][j+1];
        }
    }
    mat->ncols--;
}

/* Remove row <pos> from the matrix: pos starts from 0 */
void pluto_matrix_remove_row(PlutoMatrix *mat, int pos)
{
    int i, j;

    assert(pos <= mat->nrows-1);

    for (i=pos; i<mat->nrows-1; i++) {
        for (j=0; j<mat->ncols; j++) {
            mat->val[i][j] = mat->val[i+1][j];
        }
    }
    mat->nrows--;
}


PlutoMatrix *pluto_matrix_resize(PlutoMatrix *mat, int nrows, int ncols)
{
    int i, j;
    PlutoMatrix *newMat;

    newMat = pluto_matrix_alloc(PLMAX(nrows,mat->alloc_nrows), PLMAX(ncols,mat->alloc_ncols));

    newMat->nrows = nrows;
    newMat->ncols = ncols;

    for (i=0; i<PLMIN(newMat->nrows, mat->nrows); i++) {
        for (j=0; j<PLMIN(newMat->ncols, mat->ncols); j++) {
            newMat->val[i][j]= mat->val[i][j];
        }
    }

    pluto_matrix_free(mat);

    return newMat;
}


/* Add column to the matrix at <pos>: pos starts from 0;
 * New column is initialized to zero */
void pluto_matrix_add_col(PlutoMatrix **mat, int pos)
{
    int i, j;

    assert(pos <= (*mat)->ncols);

    if ((*mat)->ncols == (*mat)->alloc_ncols)  {
        *mat = pluto_matrix_resize(*mat, (*mat)->nrows, (*mat)->ncols+1);
    }else{
        (*mat)->ncols++;
    }

    for (j=(*mat)->ncols-2; j>=pos; j--) {
        for (i=0; i<(*mat)->nrows; i++) {
            (*mat)->val[i][j+1] = (*mat)->val[i][j];
        }
    }

    /* Initialize to zero */
    for (i=0; i<(*mat)->nrows; i++) {
        (*mat)->val[i][pos] = 0;
    }
}


/* Add row to the matrix at <pos>: pos starts from 0; row is
 * initialized to zero */
void pluto_matrix_add_row(PlutoMatrix **mat, int pos)
{
    int i, j;

    assert (pos <= (*mat)->nrows);

    if ((*mat)->nrows == (*mat)->alloc_nrows)   {
        *mat = pluto_matrix_resize(*mat, (*mat)->nrows+1, (*mat)->ncols);
    }else{
        (*mat)->nrows++;
    }

    for (i=(*mat)->nrows-2; i>=pos; i--) {
        for (j=0; j<(*mat)->ncols; j++) {
            (*mat)->val[i+1][j] = (*mat)->val[i][j];
        }
    }

    for (j=0; j<(*mat)->ncols; j++) {
        (*mat)->val[pos][j] = 0;
    }

}



/* Return a copy of src */
PlutoMatrix *pluto_matrix_copy(const PlutoMatrix *src)
{
    int i, j;

    PlutoMatrix *dest = pluto_matrix_alloc(src->alloc_nrows, src->alloc_ncols);

    for (i=0; i<src->nrows; i++)    {
        for (j=0; j<src->ncols; j++)    {
            dest->val[i][j] = src->val[i][j];
        }
    }

    dest->nrows = src->nrows;
    dest->ncols = src->ncols;

    return dest;
}

void pluto_matrix_zero_row(PlutoMatrix *mat, int pos)
{
    int j;

    assert(pos >= 0 && pos <= mat->nrows-1);

    for (j=0; j<mat->ncols; j++)    {
        mat->val[pos][j] = 0;
    }
}

void pluto_matrix_zero_col (PlutoMatrix *mat, int pos)
{
    int i;

    assert(pos >= 0 && pos <= mat->ncols-1);

    for (i=0; i<mat->nrows; i++)    {
        mat->val[i][pos] = 0;
    }
}


void pluto_matrix_read(FILE *fp, const PlutoMatrix *mat)
{
    int i, j;

    for (i=0; i<mat->nrows; i++)
        for (j=0; j<mat->ncols; j++)
            fscanf(fp, "%d", &mat->val[i][j]);
}


void pluto_matrix_print(FILE *fp, const PlutoMatrix *mat)
{
    int i, j;

    fprintf(fp, "%d %d\n", mat->nrows, mat->ncols);

    for (i=0; i<mat->nrows; i++) {
        for (j=0; j<mat->ncols; j++) {
            fprintf(fp, "%s%d ", mat->val[i][j]>=0? " ":"", mat->val[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}


/* Normalize row by its gcd */
void pluto_matrix_normalize_row(PlutoMatrix *mat, int pos)
{
    int i, j, k;

    /* Normalize mat first */
    for (i=0; i<mat->nrows; i++)    {
        if (mat->val[i][0] == 0) continue;
        int rowgcd = abs(mat->val[i][0]);
        for (j=1; j<mat->ncols; j++)    {
            if (mat->val[i][j] == 0)  break;
            rowgcd = gcd(rowgcd,abs(mat->val[i][j]));
        }
        if (i == mat->nrows) {
            if (rowgcd > 1) {
                for (k=0; k<mat->ncols; k++)    {
                    mat->val[i][k] /= rowgcd;
                }
            }
        }
    }
}


inline int lcm(int a, int b)
{
    return (a*b)/gcd(a,b);
}


inline int gcd(int a, int b)
{
    a = abs(a);
    b = abs(b);

    if (a==0)   return b;
    if (b==0)   return a;

    if (a == b) return a;

    return ((a > b)? gcd(a%b,b): gcd(a,b%a));
}


int *min_lexical(int *a, int *b, int num)   
{
    int i;

    for (i=0; i<num; i++)   {
        if (a[i] > b[i]) return b;
        else if (a[i] < b[i]) return a;
    }

    /* both are equal */
    return a;
}
