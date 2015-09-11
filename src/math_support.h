/*
 * PLuTo: An automatic parallelier and locality optimizer
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
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#ifndef _MATH_SUPPORT_H
#define _MATH_SUPPORT_H

#include <stdio.h>

#define PLMAX(a,b) ((a>=b)?(a):(b))
#define PLMIN(a,b) ((a<=b)?(a):(b))
#define PLABS(a) ((a>=0)?(a):(-a))

/* A matrix */
struct plutoMatrix{
    /* The values */
    int **val;

    int nrows;
    int ncols;

    /* Pre-allocated number of rows */
    int alloc_nrows;
    int alloc_ncols;
};
typedef struct plutoMatrix PlutoMatrix;

void pluto_matrix_print(FILE *, const PlutoMatrix *);
void pluto_matrix_read(FILE *, const PlutoMatrix *);
PlutoMatrix *pluto_matrix_alloc(int nrows, int ncols);
void pluto_matrix_free(PlutoMatrix *mat);
PlutoMatrix *pluto_matrix_dup(const PlutoMatrix *src);
void pluto_matrix_add_row(PlutoMatrix *mat, int pos);
void pluto_matrix_add_col(PlutoMatrix *mat, int pos);
void pluto_matrix_remove_row(PlutoMatrix *mat, int pos);
void pluto_matrix_remove_col(PlutoMatrix *, int);
void pluto_matrix_zero_row(PlutoMatrix *mat, int pos);
void pluto_matrix_zero_col(PlutoMatrix *mat, int pos);
void pluto_matrix_normalize_row(PlutoMatrix *mat, int pos);
void pluto_matrix_negate_row(PlutoMatrix *mat, int pos);
void pluto_matrix_initialize(PlutoMatrix *mat, int val);
void pluto_matrix_add(PlutoMatrix *mat1, const PlutoMatrix *mat2);

int lcm(int a, int b);
int gcd(int a, int b);
int *min_lexical(int *a, int *b, int num);

#endif
