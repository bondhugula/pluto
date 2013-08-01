/*
 * Polylib interface for PlutoConstraints
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
#include "constraints_polylib.h"

#include "polylib/polylib64.h"

Matrix *pluto_matrix_to_polylib(const PlutoMatrix *mat)
{
    int r, c;
    Matrix *polymat;

    polymat = Matrix_Alloc(mat->nrows, mat->ncols);

    for (r=0; r<mat->nrows; r++)    {
        for (c=0; c<mat->ncols; c++)    {
            polymat->p[r][c] = mat->val[r][c];
        }
    }

    return polymat;
}


PlutoMatrix *polylib_matrix_to_pluto(Matrix *pmat)
{
    PlutoMatrix *mat;
    int r, c;

    mat = pluto_matrix_alloc(pmat->NbRows, pmat->NbColumns);

    for (r=0; r<mat->nrows; r++)    {
        for (c=0; c<mat->ncols; c++)    {
            mat->val[r][c] = pmat->p[r][c];
        }

    }

    return mat;
}


/* Converts to a polylib polyhedron */
Polyhedron *pluto_constraints_to_polylib(const PlutoConstraints *cst)
{
    Polyhedron *pol;
    PlutoMatrix *mat;
    Matrix *polymat;

    mat = pluto_constraints_to_matrix(cst);

    polymat = pluto_matrix_to_polylib(mat);

    pol = Constraints2Polyhedron(polymat, 50);

    Matrix_Free(polymat);
    pluto_matrix_free(mat);

    if (cst->next != NULL)  {
        pol->next = pluto_constraints_to_polylib(cst->next);
    }

    return pol;
}

PlutoConstraints *polylib_to_pluto_constraints(Polyhedron *pol)
{
    PlutoConstraints *cst;

    // printf("Poly\n");
    // Polyhedron_Print(stdout, "%4d", pol);
    Matrix *polymat = Polyhedron2Constraints(pol);
     //printf("constraints\n");
     //Matrix_Print(stdout, "%4d", polymat);

    cst = polylib_matrix_to_pluto_constraints(polymat);
    Matrix_Free(polymat);

    if (pol->next != NULL) {
        cst->next = polylib_to_pluto_constraints(pol->next);
    }

    return cst;
}


/*
 * Image: if cst is a list of constraints, just its first element's image
 * is taken.
 */
PlutoConstraints *pluto_constraints_image(const PlutoConstraints *cst, const PlutoMatrix *func)
{
    // IF_DEBUG(printf("pluto const image: domain\n"););
    // IF_DEBUG(pluto_constraints_print(stdout, cst););
    // IF_DEBUG(printf("pluto const image: func\n"););
    // IF_DEBUG(pluto_matrix_print(stdout, func););
    assert(func->ncols == cst->ncols);

    PlutoConstraints *imagecst;

    Polyhedron *pol = pluto_constraints_to_polylib(cst);
    Matrix *polymat = pluto_matrix_to_polylib(func);

    // Polyhedron_Print(stdout, "%4d", pol);
    Polyhedron *image = Polyhedron_Image(pol, polymat, 2*cst->nrows);

    // Polyhedron_Print(stdout, "%4d ",  image);
    imagecst = polylib_to_pluto_constraints(image);

    Matrix_Free(polymat);
    Domain_Free(pol);
    Polyhedron_Free(image);

    // IF_DEBUG(printf("pluto const image is\n"););
    // IF_DEBUG(pluto_constraints_print(stdout, imagecst););

    return imagecst;
}

PlutoConstraints *pluto_constraints_union(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
    Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);
    Polyhedron *pol3 = DomainUnion(pol1, pol2, 50);

    PlutoConstraints *ucst = polylib_to_pluto_constraints(pol3);

    Domain_Free(pol1);
    Domain_Free(pol2);
    Domain_Free(pol3);

    return ucst;
}

PlutoConstraints *pluto_constraints_difference(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    assert(cst1->ncols == cst2->ncols);

    Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
    Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);
    Polyhedron *pol3 = DomainDifference(pol1, pol2, 50);

    PlutoConstraints *diffcst = polylib_to_pluto_constraints(pol3);

    Domain_Free(pol1);
    Domain_Free(pol2);
    Domain_Free(pol3);

    return diffcst;
}


PlutoConstraints *pluto_constraints_intersection(const PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
    Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);

    Polyhedron *pol3 = DomainIntersection(pol1, pol2, 50);

    PlutoConstraints *icst = polylib_to_pluto_constraints(pol3);

    Domain_Free(pol1);
    Domain_Free(pol2);
    Domain_Free(pol3);

    return icst;
}


/* Converts polylib matrix to pluto constraints */
PlutoConstraints *polylib_matrix_to_pluto_constraints(Matrix *polymat)
{
    int i, j;
    PlutoConstraints *cst;

    cst = pluto_constraints_alloc(polymat->NbRows, polymat->NbColumns-1);
    cst->nrows = polymat->NbRows;

    for (i=0; i<cst->nrows; i++)    {
        cst->is_eq[i]= (polymat->p[i][0] == 0)? 1: 0;
        for (j=0; j<cst->ncols; j++)    {
            cst->val[i][j] = polymat->p[i][j+1];
        }
    }

    return cst;
}

PlutoMatrix *pluto_matrix_inverse(PlutoMatrix *mat)
{
    assert(mat->nrows == mat->ncols);

    PlutoMatrix *inv;

    int dim = mat->nrows;

    Matrix *pinv = Matrix_Alloc(dim, dim);
    Matrix *pmat = pluto_matrix_to_polylib(mat);

    Matrix_Inverse(pmat, pinv);
    
    inv = polylib_matrix_to_pluto(pinv);

    Matrix_Free(pmat);
    Matrix_Free(pinv);

    return inv;
}

/* In-place difference: Like difference, but first argument is modified */
PlutoConstraints *pluto_constraints_subtract(PlutoConstraints *cst1, 
        const PlutoConstraints *cst2)
{
    PlutoConstraints *dcst = pluto_constraints_difference(cst1,cst2);
    pluto_constraints_copy(cst1,dcst);
    pluto_constraints_free(dcst);
    return cst1;
}

int pluto_constraints_are_equal(const PlutoConstraints *cst1, const PlutoConstraints *cst2) {
    PlutoConstraints *diff = pluto_constraints_difference(cst1,cst2);
    int are_constraints_equal = pluto_constraints_is_empty(diff);
    pluto_constraints_free(diff);
    if (are_constraints_equal) {
        diff = pluto_constraints_difference(cst2, cst1);
        are_constraints_equal = pluto_constraints_is_empty(diff);
        pluto_constraints_free(diff);
    }
    return are_constraints_equal;
}
