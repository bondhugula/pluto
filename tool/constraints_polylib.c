/**
 * This file is part of Pluto.
 *
 * This software is available under the MIT license. Please see LICENSE in the
 * top-level directory for details.
 *
 * Polylib interface for PlutoConstraints.
 *
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "constraints.h"
#include "constraints_polylib.h"
#include "math_support.h"
#include "pluto.h"
#include "pluto/matrix.h"

#include "polylib/polylib64.h"

Matrix *pluto_matrix_to_polylib(const PlutoMatrix *mat) {
  Matrix *polymat;

  polymat = Matrix_Alloc(mat->nrows, mat->ncols);

  for (unsigned r = 0; r < mat->nrows; r++) {
    for (unsigned c = 0; c < mat->ncols; c++) {
      polymat->p[r][c] = mat->val[r][c];
    }
  }

  return polymat;
}

PlutoMatrix *polylib_matrix_to_pluto(Matrix *pmat, PlutoContext *context) {
  PlutoMatrix *mat = pluto_matrix_alloc(pmat->NbRows, pmat->NbColumns, context);

  for (unsigned r = 0; r < mat->nrows; r++) {
    for (unsigned c = 0; c < mat->ncols; c++) {
      mat->val[r][c] = pmat->p[r][c];
    }
  }

  return mat;
}

/* Converts to a polylib polyhedron */
Polyhedron *pluto_constraints_to_polylib(const PlutoConstraints *cst) {
  Polyhedron *pol;
  PlutoMatrix *mat;
  Matrix *polymat;

  mat = pluto_constraints_to_matrix(cst);

  polymat = pluto_matrix_to_polylib(mat);

  pol = Constraints2Polyhedron(polymat, 50);

  Matrix_Free(polymat);
  pluto_matrix_free(mat);

  if (cst->next != NULL) {
    pol->next = pluto_constraints_to_polylib(cst->next);
  }

  return pol;
}

PlutoConstraints *polylib_to_pluto_constraints(Polyhedron *pol,
                                               PlutoContext *context) {
  Matrix *polymat = Polyhedron2Constraints(pol);
  PlutoConstraints *cst = polylib_matrix_to_pluto_constraints(polymat, context);
  Matrix_Free(polymat);

  if (pol->next != NULL) {
    cst->next = polylib_to_pluto_constraints(pol->next, context);
  }

  return cst;
}

/*
 * Image: if cst is a list of constraints, just its first element's image
 * is taken.
 */
PlutoConstraints *pluto_constraints_image(const PlutoConstraints *cst,
                                          const PlutoMatrix *func) {
  assert(func->ncols == cst->ncols);

  Polyhedron *pol = pluto_constraints_to_polylib(cst);
  Matrix *polymat = pluto_matrix_to_polylib(func);

  // Polyhedron_Print(stdout, "%4d", pol);
  Polyhedron *image = Polyhedron_Image(pol, polymat, 2 * cst->nrows);

  // Polyhedron_Print(stdout, "%4d ",  image);
  PlutoConstraints *imagecst =
      polylib_to_pluto_constraints(image, cst->context);

  Matrix_Free(polymat);
  Domain_Free(pol);
  Polyhedron_Free(image);

  return imagecst;
}

PlutoConstraints *pluto_constraints_union(const PlutoConstraints *cst1,
                                          const PlutoConstraints *cst2) {
  Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
  Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);
  Polyhedron *pol3 = DomainUnion(pol1, pol2, 50);

  PlutoConstraints *ucst = polylib_to_pluto_constraints(pol3, cst1->context);

  Domain_Free(pol1);
  Domain_Free(pol2);
  Domain_Free(pol3);

  return ucst;
}

PlutoConstraints *pluto_constraints_difference(const PlutoConstraints *cst1,
                                               const PlutoConstraints *cst2) {
  assert(cst1->ncols == cst2->ncols);

  Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
  Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);
  Polyhedron *pol3 = DomainDifference(pol1, pol2, 50);

  PlutoConstraints *diffcst = polylib_to_pluto_constraints(pol3, cst1->context);

  Domain_Free(pol1);
  Domain_Free(pol2);
  Domain_Free(pol3);

  return diffcst;
}

PlutoConstraints *pluto_constraints_intersection(const PlutoConstraints *cst1,
                                                 const PlutoConstraints *cst2) {
  Polyhedron *pol1 = pluto_constraints_to_polylib(cst1);
  Polyhedron *pol2 = pluto_constraints_to_polylib(cst2);

  Polyhedron *pol3 = DomainIntersection(pol1, pol2, 50);

  PlutoConstraints *icst = polylib_to_pluto_constraints(pol3, cst1->context);

  Domain_Free(pol1);
  Domain_Free(pol2);
  Domain_Free(pol3);

  return icst;
}

/* In-place intersection: first argument is modified */
PlutoConstraints *pluto_constraints_intersect(PlutoConstraints *cst1,
                                              const PlutoConstraints *cst2) {
  PlutoConstraints *icst = pluto_constraints_intersection(cst1, cst2);
  pluto_constraints_copy(cst1, icst);
  pluto_constraints_free(icst);

  return cst1;
}

/* Converts polylib matrix to pluto constraints */
PlutoConstraints *polylib_matrix_to_pluto_constraints(Matrix *polymat,
                                                      PlutoContext *context) {
  PlutoConstraints *cst =
      pluto_constraints_alloc(polymat->NbRows, polymat->NbColumns - 1, context);
  cst->nrows = polymat->NbRows;

  for (unsigned i = 0; i < cst->nrows; i++) {
    cst->is_eq[i] = (polymat->p[i][0] == 0) ? 1 : 0;
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = polymat->p[i][j + 1];
    }
  }

  return cst;
}

PlutoMatrix *pluto_matrix_inverse(PlutoMatrix *mat) {
  assert(mat->nrows == mat->ncols);

  PlutoMatrix *inv;

  int dim = mat->nrows;

  Matrix *pinv = Matrix_Alloc(dim, dim);
  Matrix *pmat = pluto_matrix_to_polylib(mat);

  Matrix_Inverse(pmat, pinv);

  inv = polylib_matrix_to_pluto(pinv, mat->context);

  Matrix_Free(pmat);
  Matrix_Free(pinv);

  return inv;
}

/* In-place difference: Like difference, but first argument is modified */
PlutoConstraints *pluto_constraints_subtract(PlutoConstraints *cst1,
                                             const PlutoConstraints *cst2) {
  PlutoConstraints *dcst = pluto_constraints_difference(cst1, cst2);
  pluto_constraints_copy(cst1, dcst);
  pluto_constraints_free(dcst);
  return cst1;
}

int pluto_constraints_are_equal(const PlutoConstraints *cst1,
                                const PlutoConstraints *cst2) {
  PlutoConstraints *diff = pluto_constraints_difference(cst1, cst2);
  int are_constraints_equal = pluto_constraints_is_empty(diff);
  pluto_constraints_free(diff);
  if (are_constraints_equal) {
    diff = pluto_constraints_difference(cst2, cst1);
    are_constraints_equal = pluto_constraints_is_empty(diff);
    pluto_constraints_free(diff);
  }
  return are_constraints_equal;
}
