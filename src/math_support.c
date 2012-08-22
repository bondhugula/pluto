/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
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
 *
 * nrows and ncols initialized to allocated number of rows and cols
 */
PlutoMatrix *pluto_matrix_alloc(int alloc_nrows, int alloc_ncols)
{
    int i;
    PlutoMatrix *mat;

    assert(alloc_nrows >= 0);
    assert(alloc_ncols >= 0);

    mat = (PlutoMatrix *) malloc(sizeof(PlutoMatrix));

    mat->val = (int **) malloc(PLMAX(alloc_nrows,1)*sizeof(int *));

    mat->alloc_nrows = PLMAX(alloc_nrows,1);
    mat->alloc_ncols = PLMAX(alloc_ncols,1);

    for (i=0; i<mat->alloc_nrows; i++) {
        mat->val[i] = (int *) malloc(mat->alloc_ncols*sizeof(int));
    }

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


/* Non-destructive resize */
void pluto_matrix_resize(PlutoMatrix *mat, int nrows, int ncols)
{
    int i;

    int alloc_nrows = PLMAX(nrows,mat->alloc_nrows);
    int alloc_ncols = PLMAX(ncols,mat->alloc_ncols);

    mat->val = (int **) realloc(mat->val, alloc_nrows*sizeof(int *));

    for (i=mat->alloc_nrows; i<alloc_nrows; i++)    {
        mat->val[i] = NULL;
    }

    for (i=0; i<alloc_nrows; i++) {
        mat->val[i] = (int *) realloc(mat->val[i], alloc_ncols*sizeof(int));
    }

    mat->alloc_nrows = alloc_nrows;
    mat->alloc_ncols = alloc_ncols;

    mat->nrows = nrows;
    mat->ncols = ncols;
}


/* Add column to the matrix at <pos>: pos starts from 0;
 * New column is initialized to zero */
void pluto_matrix_add_col(PlutoMatrix *mat, int pos)
{
    int i, j;

    assert(pos >= 0 && pos <= mat->ncols);

    if (mat->ncols == mat->alloc_ncols)  {
        pluto_matrix_resize(mat, mat->nrows, mat->ncols+1);
    }else{
        mat->ncols++;
    }

    for (j=mat->ncols-2; j>=pos; j--) {
        for (i=0; i<mat->nrows; i++) {
            mat->val[i][j+1] = mat->val[i][j];
        }
    }

    /* Initialize to zero */
    for (i=0; i<mat->nrows; i++) {
        mat->val[i][pos] = 0;
    }
}

/* Negate entire row; pos is 0-indexed */
void pluto_matrix_negate_row(PlutoMatrix *mat, int pos)
{
    int j;

    for (j=0; j<mat->ncols; j++)    {
        mat->val[pos][j] = -mat->val[pos][j];
    }
}

/* Add rows of mat2 to mat1 */
void pluto_matrix_add(PlutoMatrix *mat1, const PlutoMatrix *mat2)
{
    int i, j;

    assert(mat1->ncols == mat2->ncols);

    pluto_matrix_resize(mat1, mat1->nrows+mat2->nrows, mat1->ncols);

    for (i=mat1->nrows-mat2->nrows; i<mat1->nrows; i++) {
        for (j=0; j<mat1->ncols; j++)   {
            mat1->val[i][j] = mat2->val[i-(mat1->nrows-mat2->nrows)][j];
        }
    }
}


/* Add row to the matrix at <pos>: pos starts from 0; row is
 * initialized to zero */
void pluto_matrix_add_row(PlutoMatrix *mat, int pos)
{
    int i, j;

    assert(mat != NULL);
    assert(pos <= mat->nrows);

    if (mat->nrows == mat->alloc_nrows)   {
        pluto_matrix_resize(mat, mat->nrows+1, mat->ncols);
    }else{
        mat->nrows++;
    }

    for (i=mat->nrows-2; i>=pos; i--) {
        for (j=0; j<mat->ncols; j++) {
            mat->val[i+1][j] = mat->val[i][j];
        }
    }

    for (j=0; j<mat->ncols; j++) {
        mat->val[pos][j] = 0;
    }

}

void pluto_matrix_interchange_rows(PlutoMatrix *mat, int r1, int r2)
{
    int tmp, j;

    for (j=0; j<mat->ncols; j++) {
        tmp = mat->val[r1][j];
        mat->val[r1][j] = mat->val[r2][j];
        mat->val[r2][j] = tmp;
    }
}


void pluto_matrix_interchange_cols(PlutoMatrix *mat, int c1, int c2)
{
    int tmp, i;

    for (i=0; i<mat->nrows; i++) {
        tmp = mat->val[i][c1];
        mat->val[i][c1] = mat->val[i][c2];
        mat->val[i][c2] = tmp;
    }
}

/* Move column from position c1 to c2 */
void pluto_matrix_move_col(PlutoMatrix *mat, int c1, int c2)
{
    int j;

    if (c1 < c2) {
        for (j=c1; j<c2; j++) {
            pluto_matrix_interchange_cols(mat, j, j+1);
        }
    }else{
        for (j=c1; j>c2; j--) {
            pluto_matrix_interchange_cols(mat, j, j-1);
        }
    }
}



/* Return a duplicate of src */
PlutoMatrix *pluto_matrix_dup(const PlutoMatrix *src)
{
    int i, j;

    assert(src != NULL);

    PlutoMatrix *dup = pluto_matrix_alloc(src->alloc_nrows, src->alloc_ncols);

    for (i=0; i<src->nrows; i++)    {
        for (j=0; j<src->ncols; j++)    {
            dup->val[i][j] = src->val[i][j];
        }
    }

    dup->nrows = src->nrows;
    dup->ncols = src->ncols;

    return dup;
}

/* Initialize matrix with val */
void pluto_matrix_initialize(PlutoMatrix *mat, int val)
{
    int i, j;

    for (i=0; i<mat->nrows; i++)    {
        for (j=0; j<mat->ncols; j++)    {
            mat->val[i][j] = val;
        }
    }
}

/* Return an identity matrix of size: size x size */
PlutoMatrix *pluto_matrix_identity(int size)
{
    int i;

    PlutoMatrix *mat = pluto_matrix_alloc(size, size);
    pluto_matrix_initialize(mat, 0);

    for (i=0; i<size; i++)  {
        mat->val[i][i] = 1;
    }

    return mat;
}


/* Zero out a row */
void pluto_matrix_zero_row(PlutoMatrix *mat, int pos)
{
    int j;

    assert(pos >= 0 && pos <= mat->nrows-1);

    for (j=0; j<mat->ncols; j++)    {
        mat->val[pos][j] = 0;
    }
}

/* Zero out a column */
void pluto_matrix_zero_col(PlutoMatrix *mat, int pos)
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

PlutoMatrix *pluto_matrix_input(FILE *fp)
{
    int i, j, nrows, ncols;
    fscanf(fp, "%d %d", &nrows, &ncols);

    PlutoMatrix *mat = pluto_matrix_alloc(nrows, ncols);

    for (i=0; i<mat->nrows; i++)
        for (j=0; j<mat->ncols; j++)
            fscanf(fp, "%d", &mat->val[i][j]);

    return mat;
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
    // fprintf(fp, "\n");
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


/* pos: 0-indexed */
void gaussian_eliminate_var(PlutoMatrix *mat, int pos)
{
    int r, r2, c;
    int factor1, factor2;

    // printf("Before gaussian eliminate\n");
    // pluto_matrix_print(stdout, mat);
    // printf("eliminate: %d\n", pos);

    for (r=0; r<mat->nrows; r++)    {
        if (mat->val[r][pos] != 0)  {
            for (r2=0; r2<mat->nrows; r2++) {
                if (r2 == r) continue;
                if (mat->val[r2][pos] != 0) {
                    factor1 = lcm(abs(mat->val[r][pos]),abs(mat->val[r2][pos]))/mat->val[r2][pos];
                    factor2 = lcm(abs(mat->val[r][pos]),abs(mat->val[r2][pos]))/mat->val[r][pos];
                    for (c=0; c<mat->ncols; c++) {
                        // printf("%d %d\n", mat->val[r2][pos], mat->val[r][pos]);
                        // printf("%d\n", factor1);
                        // printf("%d\n", factor2);
                        mat->val[r2][c] = mat->val[r2][c]*factor1
                            - mat->val[r][c]*factor2;
                    }
                }
            }
            pluto_matrix_remove_row(mat, r);
            break;
        }
    }
    // printf("After gaussian eliminate\n");
    // pluto_matrix_print(stdout, mat);

    pluto_matrix_remove_col(mat, pos);
}

/* Eliminate variables from start to end (inclusive); start is 0-indexed
 */
void gaussian_eliminate(PlutoMatrix *mat, int start, int num_elim)
{
    int i;

    for (i=0; i<num_elim; i++) {
        gaussian_eliminate_var(mat, start);
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

/* Free returned string with free */
char *concat(const char *prefix, const char *suffix)
{
    char *concat = malloc(strlen(prefix) + strlen(suffix) + 1);
    sprintf(concat, "%s%s", prefix, suffix);
    return concat;
}

PlutoMatrix *pluto_matrix_product(const PlutoMatrix *mat1, 
        const PlutoMatrix *mat2)
{
    assert(mat1->ncols == mat2->nrows);

    int i, j, k;

    PlutoMatrix *mat3 = pluto_matrix_alloc(mat1->nrows, mat2->ncols);

    for (i=0; i<mat1->nrows; i++)   {
        for (j=0; j<mat2->ncols; j++)   {
            mat3->val[i][j] = 0;
            for (k=0; k<mat1->ncols; k++)   {
                mat3->val[i][j] += mat1->val[i][k]*mat2->val[k][j];
            }
        }
    }
    return mat3;
}

/* Converts matrix to row-echelon form in-place */
PlutoMatrix *pluto_matrix_to_row_echelon(PlutoMatrix *mat)
{
    int i, j, k, r, _lcm, factor1;

    r=0;
    for (i=0; i< PLMIN(mat->ncols,mat->nrows); i++)  {
        //pluto_matrix_print(stdout, sched);
        if (mat->val[r][i] == 0) {
            for (k=r+1; k<mat->nrows; k++) {
                if (mat->val[k][i] != 0) break;
            }
            if (k<mat->nrows)    {
                pluto_matrix_interchange_rows(mat, r, k);
            }
        }
        if (mat->val[r][i] != 0) {
            for (k=r+1; k<mat->nrows; k++) {
                printf("i=%d, r=%d\n", i, r);
                if (mat->val[k][i] == 0) continue;
                _lcm = lcm(mat->val[k][i], mat->val[r][i]);
                factor1 = _lcm/mat->val[k][i];
                for (j=i; j<mat->ncols; j++) {
                    mat->val[k][j] = mat->val[k][j]*factor1
                        - mat->val[r][j]*(_lcm/mat->val[r][i]);
                }
            }
            r++;
        }
    }

    return mat;
}

int pluto_matrix_get_rank(const PlutoMatrix *mat)
{
    int i, j, null, rank;

    PlutoMatrix *re = pluto_matrix_to_row_echelon(pluto_matrix_dup(mat));

    // pluto_matrix_print(stdout, re);

    null = 0;
    for (i=0; i<re->nrows; i++) {
        int sum = 0;
        for (j=0; j<re->ncols; j++) {
            sum += abs(re->val[i][j]);
        }
        if (sum==0) null++;
    }
    rank = re->nrows - null;
    pluto_matrix_free(re);
    return rank;
}
