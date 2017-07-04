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

#include "math_support.h"
#include "constraints.h"

#include "isl/val.h"
#include "isl/val_gmp.h"
#include "isl/deprecated/int.h"

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

    mat->val = (int64 **) malloc(PLMAX(alloc_nrows,1)*sizeof(int64 *));

    mat->alloc_nrows = PLMAX(alloc_nrows,1);
    mat->alloc_ncols = PLMAX(alloc_ncols,1);

    for (i=0; i<mat->alloc_nrows; i++) {
        mat->val[i] = (int64 *) malloc(mat->alloc_ncols*sizeof(int64));
    }

    mat->nrows = alloc_nrows;
    mat->ncols = alloc_ncols;

    return mat;
}

void pluto_matrix_free(PlutoMatrix *mat)
{
    int i;

    if (mat) {
        for (i=0; i<mat->alloc_nrows; i++) {
            free(mat->val[i]);
        }

        free(mat->val);
        free(mat);
    }
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

    mat->val = (int64 **) realloc(mat->val, alloc_nrows*sizeof(int64 *));

    for (i=mat->alloc_nrows; i<alloc_nrows; i++)    {
        mat->val[i] = NULL;
    }

    for (i=0; i<alloc_nrows; i++) {
        mat->val[i] = (int64 *) realloc(mat->val[i], alloc_ncols*sizeof(int64));
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

void pluto_matrix_negate(PlutoMatrix *mat)
{
    int r;
    for (r=0; r<mat->nrows; r++)    {
        pluto_matrix_negate_row(mat, r);
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
void pluto_matrix_set(PlutoMatrix *mat, int val)
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
    pluto_matrix_set(mat, 0);

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
            fscanf(fp, "%lld", &mat->val[i][j]);
}

PlutoMatrix *pluto_matrix_input(FILE *fp)
{
    int i, j, nrows, ncols;
    fscanf(fp, "%d %d", &nrows, &ncols);

    PlutoMatrix *mat = pluto_matrix_alloc(nrows, ncols);

    for (i=0; i<mat->nrows; i++)
        for (j=0; j<mat->ncols; j++)
            fscanf(fp, "%lld", &mat->val[i][j]);

    return mat;
}




void pluto_matrix_print(FILE *fp, const PlutoMatrix *mat)
{
    int i, j;

    fprintf(fp, "%d %d\n", mat->nrows, mat->ncols);

    for (i=0; i<mat->nrows; i++) {
        for (j=0; j<mat->ncols; j++) {
            fprintf(fp, "%s%lld ", mat->val[i][j]>=0? " ":"", mat->val[i][j]);
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
        int rowgcd = llabs(mat->val[i][0]);
        for (j=1; j<mat->ncols; j++)    {
            if (mat->val[i][j] == 0)  break;
            rowgcd = gcd(rowgcd,llabs(mat->val[i][j]));
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
                    factor1 = lcm(llabs(mat->val[r][pos]),llabs(mat->val[r2][pos]))/mat->val[r2][pos];
                    factor2 = lcm(llabs(mat->val[r][pos]),llabs(mat->val[r2][pos]))/mat->val[r][pos];
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


inline int64 lcm(int64 a, int64 b)
{
    if (a*b == 0) return 0;
    return (a*b)/gcd(a,b);
}


/* Assuming both args are not zero */
inline int64 gcd(int64 a, int64 b)
{
    a = llabs(a);
    b = llabs(b);

    /* If at least one of them is zero */
    if (a*b ==0)   return a+b;

    if (a == b) return a;

    return ((a > b)? gcd(a%b,b): gcd(a,b%a));
}


int64 *min_lexical(int64 *a, int64 *b, int64 num)   
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

/* Rank of the matrix */
int pluto_matrix_get_rank(const PlutoMatrix *mat)
{
    int i, j, null, rank;

    PlutoMatrix *re = pluto_matrix_to_row_echelon(pluto_matrix_dup(mat));

    // pluto_matrix_print(stdout, re);

    null = 0;
    for (i=0; i<re->nrows; i++) {
        int sum = 0;
        for (j=0; j<re->ncols; j++) {
            sum += llabs(re->val[i][j]);
        }
        if (sum==0) null++;
    }
    rank = re->nrows - null;
    pluto_matrix_free(re);
    return rank;
}


void pluto_matrix_swap_rows(PlutoMatrix *mat, int r1, int r2)
{
    int64 tmp;
    int j;

    for (j=0; j<mat->ncols; j++) {
        tmp = mat->val[r2][j];
        mat->val[r2][j] = mat->val[r1][j];
        mat->val[r1][j] = tmp;
    }

}

/* Reverse the order of rows in the matrix */
void pluto_matrix_reverse_rows(PlutoMatrix *mat)
{
    int i;

    for (i=0; i<mat->nrows/2; i++) {
        pluto_matrix_swap_rows(mat, i, mat->nrows-1-i);
    }
}

/*
 * Construct a PlutoMatrix with the same content as the given isl_mat.
 */
PlutoMatrix *pluto_matrix_from_isl_mat(__isl_keep isl_mat *mat)
{
    int i, j;
    int rows, cols;
    PlutoMatrix *pluto;

    rows = isl_mat_rows(mat);
    cols = isl_mat_cols(mat);
    pluto = pluto_matrix_alloc(rows, cols);

    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j) {
            isl_val *v = isl_mat_get_element_val(mat, i, j);
            pluto->val[i][j] = isl_val_get_num_si(v);
            isl_val_free(v);
        }

    return pluto;
}

/*
 * Pretty prints a one-dimensional affine function
 * ndims: number of variables
 * func should have ndims+1 elements (affine function)
 * vars: names of the ndims variables; if NULL, x0, x1, ... are used
 */
void pluto_affine_function_print(FILE *fp, int64 *func, int ndims, char **vars)
{
    char *var[ndims];
    int j;

    for (j=0; j<ndims; j++)  {
        if (vars && vars[j]) {
            var[j] = strdup(vars[j]);
        }else{
            var[j] = malloc(5);
            sprintf(var[j], "x%d", j+1);
        }
    }

    int first = 0;
    for (j=0; j<ndims; j++)   {
        if (func[j] == 1)  {
            if (first) fprintf(fp, "+");
            fprintf(fp, "%s", var[j]);
        }else if (func[j] == -1)  {
            fprintf(fp, "-%s", var[j]);
        }else if (func[j] != 0)  {
            if (func[j] > 0) {
                fprintf(fp, "%s%lld%s", first?"+":"", func[j], var[j]);
            }else{
                fprintf(fp, "%lld%s", func[j], var[j]);
            }
        }
        if (func[j] != 0) first = 1;
    }
    /* Constant part */
    if (func[ndims] >= 1)  {
        if (first) fprintf(fp, "+");
        fprintf(fp, "%lld", func[ndims]);
    }else if (func[ndims] <= -1){
        fprintf(fp, "%lld", func[ndims]);
    }else{
        /* 0 */
        if (!first) fprintf(fp, "0");
    }

    for (j=0; j<ndims; j++)  {
        free(var[j]);
    }
}

/* Returned string should be freed with malloc */
char *pluto_affine_function_sprint(int64 *func, int ndims, char **vars)
{
    char *var[ndims], *out;
    int j, n;

    /* max 5 chars for var, 3 for coefficient + 1 if ndims is 0 + 1 null char */
    n = 9*ndims + 1 + 1;
    out = malloc(n);
    *out = '\0';

    for (j=0; j<ndims; j++)  {
        if (vars && vars[j]) {
            var[j] = strdup(vars[j]);
        }else{
            var[j] = malloc(5);
            sprintf(var[j], "x%d", j+1);
        }
    }

    int first = 0;
    for (j=0; j<ndims; j++)   {
        if (func[j] == 1)  {
            if (first) strcat(out, "+");
            snprintf(out+strlen(out), 5, "%s", var[j]);
        }else if (func[j] == -1)  {
            snprintf(out+strlen(out), 6, "-%s", var[j]);
        }else if (func[j] != 0)  {
            if (func[j] >= 1) {
                snprintf(out+strlen(out), 9, "%s%lld%s", first?"+":"", func[j], var[j]);
            }else{
                snprintf(out+strlen(out), 8, "%lld%s", func[j], var[j]);
            }
        }
        if (func[j] != 0) first = 1;
    }
    /* Constant part */
    if (func[ndims] >= 1)  {
        if (first) strcat(out, "+");
        snprintf(out+strlen(out), 3, "%lld", func[ndims]);
    }else if (func[ndims] <= -1){
        snprintf(out+strlen(out), 3, "%lld", func[ndims]);
    }else{
        /* 0 */
        if (!first) strcat(out, "0");
    }

    for (j=0; j<ndims; j++)  {
        free(var[j]);
    }

    return out;
}



/*
 * Convert an isl affine expression to Pluto function
 */
int isl_aff_to_pluto_func(__isl_take isl_set *set, __isl_take isl_aff *aff,
        void *user)
{
    int i, j, npar;

    npar = isl_aff_dim(aff, isl_dim_param);

    PlutoMatrix **mat_p = (PlutoMatrix **) user;
    if (*mat_p != NULL) pluto_matrix_free(*mat_p);
    *mat_p = pluto_matrix_alloc(1, isl_aff_dim(aff, isl_dim_in) + npar + 1);
    PlutoMatrix *mat = *mat_p;

    if (isl_aff_dim(aff, isl_dim_div) >= 1) {
        isl_aff *div = isl_aff_get_div(aff, 0);
        isl_val *v = isl_aff_get_denominator_val(div);
        isl_aff_free(div);
        if (!isl_val_is_one(v)) {
            pluto_matrix_zero_row(mat, 0);
            isl_val_free(v);
            isl_set_free(set);
            isl_aff_free(aff);
            return 0;
        }
        isl_val_free(v);
    }

    for (i=0; i<isl_aff_dim(aff, isl_dim_in); i++) {
        isl_val *v = isl_aff_get_coefficient_val(aff, isl_dim_in, i);
        mat->val[0][i] = isl_val_get_num_si(v);
        isl_val_free(v);
    }
    for (j=0; j<npar; i++,j++) {
        isl_val *v =  isl_aff_get_coefficient_val(aff, isl_dim_param, j);
        mat->val[0][i] = isl_val_get_num_si(v);
        isl_val_free(v);
    }
    isl_val *v = isl_aff_get_constant_val(aff);
    mat->val[0][i] = isl_val_get_num_si(v);
    isl_val_free(v);
    //pluto_matrix_print(stdout, mat);

    isl_set_free(set);
    isl_aff_free(aff);

    return 0;
}



/*
 * Is row r1 of mat1 parallel to row r2 of mat2
 */
int pluto_vector_is_parallel(PlutoMatrix *mat1, int r1, PlutoMatrix *mat2, int r2)
{
    int num, den, j;

    assert(mat1->ncols == mat2->ncols);

    num = 0;
    den = 0;

    for (j=0; j<mat1->ncols; j++) {
        if (mat1->val[r1][j] == 0 && mat2->val[r2][j] == 0) continue;
        if (mat1->val[r1][j] != 0 && mat2->val[r2][j] == 0) return 0;
        if (mat1->val[r1][j] == 0 && mat2->val[r2][j] != 0) return 0;

        /* num and den are always non-zero */
        if (num == 0) {
            /* first time */
            num = mat1->val[r1][j];
            den = mat2->val[r2][j];
        }else{
            if (num*mat2->val[r2][j] != den*mat1->val[r1][j]) return 0;
        }
    }

    return 1;
}


int pluto_vector_is_normal(PlutoMatrix *mat1, int r1, PlutoMatrix *mat2, int r2)
{
    int j, dot;

    assert(mat1->ncols == mat2->ncols);

    dot = 0;
    for (j=0; j<mat1->ncols; j++) {
        dot +=  mat1->val[r1][j]*mat2->val[r2][j];
    }

    return (dot == 0);
}

/* Convert from mpz to signed long long */
void mpz_set_sll(mpz_t n, long long sll)
{   
    /* n = (int)sll >> 32 */
    mpz_set_si(n, (int)(sll >> 32));  
    /* n <<= 32 */
    mpz_mul_2exp(n, n, 32 );
    /* n += (unsigned int)sll */
    mpz_add_ui(n, n, (unsigned int)sll); 
}

/* Convert from mpz to unsigned long long */
void mpz_set_ull(mpz_t n, unsigned long long ull)
{
    /* n = (unsigned int)(ull >> 32) */
    mpz_set_ui(n, (unsigned int)(ull >> 32)); 
    /* n <<= 32 */
    mpz_mul_2exp(n, n, 32);                   
    /* n += (unsigned int)ull */
    mpz_add_ui(n, n, (unsigned int)ull);      
}


/* Does the mpz value fit in a long long? */
static int mpz_fits_ll(mpz_t z)
{
    int sign;

    mpz_t tmp;
    mpz_init(tmp);
    /* tmp = (upper bits of r) */
    if (mpz_sgn(z) > 0) {
        mpz_tdiv_q_2exp(tmp, z, 63);
    }else{
        mpz_add_ui(tmp, z, 1);
        mpz_tdiv_q_2exp(tmp, tmp, 63);
    }
    sign = mpz_sgn(tmp);
    mpz_clear(tmp);
    return sign == 0;

}

/* Extract the numerator of a rational value "v" as a long long int.
 *
 * If "v" is not a rational value, then the result is undefined.
 */
long long isl_val_get_num_ll(__isl_keep isl_val *v)
{
    unsigned lo, hi;
    mpz_t z, tmp;
    int sign;
    long long result;

    mpz_init(tmp);
    // isl_val_get_num_gmp(v, tmp);
    // gmp_printf("isl_int: %Zd\n", tmp );
    // gmp_printf("isl_int hex: %Zx\n", tmp);
    
    if (!v)
        return 0;
    if (!isl_val_is_rat(v))
        isl_die(isl_val_get_ctx(v), isl_error_invalid,
                "expecting rational value", return 0);

    mpz_init(z);
    isl_val_get_num_gmp(v, z);

    if (!mpz_fits_ll(z)) {
        // gmp_printf("isl_int: %Zd\n", z);
        // gmp_printf("isl_int hex: %Zx\n", z);
        printf("[pluto_math_support] numerator too large; returning largest/smallest signed 64-bit number\n");
        sign = mpz_sgn(z);
        mpz_clear(z);

        /* 2^63-1 is the largest positive signed number for 64-bit */
        /* -2^63 is the largest negative signed number for 64-bit */
        return sign? (long long) ((1ULL<<63) -1): (long long) ((1ULL << 63));
    }

    /* tmp = (lower 64 bits of r) */
    mpz_mod_2exp(tmp, z, 64);
    mpz_clear(z);

    /* lo = tmp & 0xffffffff */
    lo = mpz_get_ui(tmp);

    /* tmp >>= 32 */
    mpz_div_2exp(tmp, tmp, 32);

    /* hi = tmp & 0xffffffff */
    hi = mpz_get_ui(tmp);

    mpz_clear(tmp);
    result = (long long)( (((unsigned long long)hi) << 32) + lo);

    // printf("result: %lld\n", result);

    return result;
}
