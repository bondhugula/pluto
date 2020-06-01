/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE.MIT
 * in the top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>
#include <unistd.h>

#include "constraints.h"
#include "math_support.h"
#include "pluto.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"

#include "piplib/piplib64.h"

#define UB 0
#define LB 1
#define NB 2

PipMatrix *pip_matrix_populate(int64_t **cst, int nrows, int ncols);

/*
 * Allocate with a max size of max_rows and max_cols;
 * initialized all to zero and is_eq to 0
 *
 * nrows set to 0, and ncols to max_cols, i.e., the initially allocated
 * constraints correspond to the universe (no constraints)
 *
 * As rows are added, increase nrows
 */
PlutoConstraints *pluto_constraints_alloc(int max_rows, int max_cols,
                                          PlutoContext *context) {
  assert(context && "null context");
  PlutoConstraints *cst = (PlutoConstraints *)malloc(sizeof(PlutoConstraints));
  cst->context = context;

  size_t size =
      ((size_t)PLMAX(1, max_rows)) * PLMAX(1, max_cols) * sizeof(int64_t);

  cst->buf = (int64_t *)malloc(size);

  if (cst->buf == NULL) {
    fprintf(stderr,
            "[pluto] ERROR: Not enough memory to allocate constraints\n");
    fprintf(stderr, "[pluto] %zd bytes needed\n", size);
    exit(1);
  }

  bzero(cst->buf, size);

  cst->is_eq = (int *)malloc(max_rows * sizeof(int));
  bzero(cst->is_eq, max_rows * sizeof(int));

  cst->val = (int64_t **)malloc(max_rows * sizeof(int64_t *));

  int i;
  for (i = 0; i < max_rows; i++) {
    cst->val[i] = &cst->buf[i * max_cols];
  }

  cst->alloc_nrows = max_rows;
  cst->alloc_ncols = max_cols;
  cst->names = NULL;
  cst->next = NULL;

  cst->nrows = 0;
  cst->ncols = max_cols;

  return cst;
}

/* Initialize entire *allocated* constraints to zero; everything is also
 * initialized to inequality >= 0 */
void pluto_constraints_zero(PlutoConstraints *cst) {
  bzero(cst->buf, cst->alloc_ncols * cst->alloc_nrows * sizeof(int64_t));
  bzero(cst->is_eq, cst->alloc_nrows * sizeof(int));
}

void pluto_constraints_free(PlutoConstraints *cst) {
  if (cst == NULL)
    return;

  free(cst->buf);
  free(cst->val);
  free(cst->is_eq);
  if (cst->names) {
    for (int i = 0; i < (int)cst->ncols - 1; i++) {
      free(cst->names[i]);
    }
    free(cst->names);
  }

  if (cst->next != NULL)
    pluto_constraints_free(cst->next);

  free(cst);
}

/* Non-destructive resize */
void pluto_constraints_resize(PlutoConstraints *cst, int nrows, int ncols) {
  pluto_constraints_resize_single(cst, nrows, ncols);

  if (cst->next != NULL) {
    pluto_constraints_resize(cst->next, nrows, ncols);
  }
}

/*
 * Non-destructive resize (single element of constraint list)
 * Not just change the allocated size, but also increase/reduce the number
 * of constraints/dimensions to match nrows and ncols
 */
void pluto_constraints_resize_single(PlutoConstraints *cst, unsigned nrows,
                                     unsigned ncols) {
  PlutoConstraints *newCst =
      pluto_constraints_alloc(PLMAX(nrows, cst->alloc_nrows),
                              PLMAX(ncols, cst->alloc_ncols), cst->context);

  newCst->nrows = nrows;
  newCst->ncols = ncols;

  for (unsigned i = 0; i < PLMIN(newCst->nrows, cst->nrows); i++) {
    for (unsigned j = 0; j < PLMIN(newCst->ncols, cst->ncols); j++) {
      newCst->val[i][j] = cst->val[i][j];
    }
    newCst->is_eq[i] = cst->is_eq[i];
  }

  free(cst->val);
  free(cst->buf);
  free(cst->is_eq);

  if (cst->names) {
    for (int i = ncols - 1; i < (int)cst->ncols - 1; i++) {
      free(cst->names[i]);
    }
  }

  if (cst->names) {
    cst->names = (char **)realloc(cst->names, (ncols - 1) * sizeof(char *));
    for (int i = cst->ncols - 1; i < (int)ncols - 1; i++) {
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

/*
 * The affine form of the Farkas lemma
 *
 * cst is the domain on which the affine form described in phi is non-negative
 *
 * Returns: constraints on the variables that correspond to the columns of
 * \phi (each row of phi is an affine function of these variables). The output
 * is thus a  constraint set with phi->ncols-1 variables. Farkas multipliers
 * are eliminated  by Fourier-Motzkin elimination. In effect, this allows one
 * to linearize the constraint describing \phi.
 *
 * The rows of phi correspond to coefficients of variables in 'dom' in
 * that order with the last row of phi representing the translation part of
 * the affine form. The number of rows in phi is thus the same as the number
 * of columns in dom (number of dom dimensions + 1).
 *
 * Each row of phi itself is an affine function of a set of variables
 * (phi->ncols-1 variables); the last column corresponds to the constant.
 *
 * Eg:
 * (c_1 + c_2)*i + (c_2 - c3)*j + (c1 + c2 + c3 + 1) >= 0 over a domain (dom) on
 * (i,j), say {(i,j)| 0 <= i <= N-1 and j = i+1 }
 *
 * Here, phi would be
 *
 * [1 1 0  0] <-- i
 * [0 1 -1 0] <-- j
 * [1 1 1  1] <-- 1
 *
 * Let cst have faces (inequalities representing non-negative half-spaces) f1,
 * f2, ..., fn
 *
 * The affine form of the Farkas lemma states that
 *
 * (c_1 + c_2)*i + (c_2 - c3)*j + (c1 + c2 + c3 + 1) = \lambda_1*f1 +
 * \lambda_2*f2 + ... + \lambda_n*fn + \lambda_0,
 * with all \lambda_i >= 0
 *
 * Eliminate Farkas multipliers by FM and return constraints in c_1, c_2, ...
 *
 * */
PlutoConstraints *farkas_lemma_affine(const PlutoConstraints *dom,
                                      const PlutoMatrix *phi) {
  /* Only for a convex constraint set */
  assert(dom->next == NULL);
  PlutoContext *context = dom->context;

  assert(phi->nrows == dom->ncols);

  IF_MORE_DEBUG(printf("[farkas_lemma_affine]\n"););

  /* Convert everything into inequalities of >= 0 form */
  PlutoConstraints *idom = pluto_constraints_to_pure_inequalities_single(dom);

  /* Add a trivial row (1 >= 0) for the translation Farkas
   * multiplier (\lambda_0) so that the non-negative linear
   * combination of the faces is modeled naturally below */
  pluto_constraints_add_inequality(idom);
  idom->val[idom->nrows - 1][idom->ncols - 1] = 1;

  /*
   * Farkas space
   * idom->ncols equalities (one for each of the idom->ncols-1 variables)
   * and one for the constant part, followed by idom->nrows constraints for
   * the farkas multipliers *
   * idom->nrows is the number of Farkas multipliers
   * format: [phi->ncols-1 vars, idom->nrows farkas multipliers, const]
   * the translation farkas multiplier appears last
   *
   *    Eg: [c_1, c_2, c_3, l_1, l_2, ..., l_n, l_0, 1]
   */
  PlutoConstraints *farkas = pluto_constraints_alloc(
      idom->ncols + idom->nrows, phi->ncols + idom->nrows, context);
  farkas->nrows = idom->ncols + idom->nrows;

  int farkas_offset = phi->ncols - 1;

  /* First idom->ncols equalities */
  for (unsigned i = 0; i < idom->ncols; i++) {
    farkas->is_eq[i] = 1;
    for (unsigned j = 0; j < (int)phi->ncols - 1; j++) {
      farkas->val[i][j] = phi->val[i][j];
    }
    for (unsigned j = 0; j < idom->nrows; j++) {
      farkas->val[i][farkas_offset + j] = -idom->val[j][i];
    }
    farkas->val[i][farkas_offset + idom->nrows] = phi->val[i][phi->ncols - 1];
  }

  /* All farkas multipliers are non-negative */
  for (unsigned j = 0; j < idom->nrows; j++) {
    farkas->is_eq[idom->ncols + j] = 0;
    farkas->val[idom->ncols + j][farkas_offset + j] = 1;
  }

  for (unsigned i = 0; i < idom->nrows; i++) {
    int best_elim =
        pluto_constraints_best_elim_candidate(farkas, idom->nrows - i);
    IF_MORE_DEBUG(printf("[farkas_lemma_affine] eliminating multiplier %d "
                         "(c_%c) from %d constraints\n",
                         i, 'i' + best_elim, farkas->nrows));
    fourier_motzkin_eliminate_smart(farkas, best_elim);
  }
  assert(farkas->ncols == phi->ncols);

  pluto_constraints_free(idom);

  return farkas;
}

/// Adds cs1 and cs2 and puts them in cs1 -> returns cs1 itself for convenience.
/// This function is only for  single element lists. For multiple element lists,
/// use pluto_constraints_intersect.
PlutoConstraints *pluto_constraints_add(PlutoConstraints *cst1,
                                        const PlutoConstraints *cst2) {
  assert(cst2 != NULL);
  assert(cst1->ncols == cst2->ncols);
  assert(cst1->next == NULL);
  assert(cst2->next == NULL);

  if (cst1->nrows + cst2->nrows > cst1->alloc_nrows) {
    pluto_constraints_resize(cst1, cst1->nrows + cst2->nrows, cst1->ncols);
  } else {
    cst1->nrows = cst1->nrows + cst2->nrows;
  }

  for (int i = 0; i < (int)cst2->nrows; i++) {
    memcpy(cst1->val[(int)cst1->nrows - (int)cst2->nrows + i], cst2->val[i],
           cst1->ncols * sizeof(int64_t));
  }

  memcpy(&cst1->is_eq[cst1->nrows - cst2->nrows], cst2->is_eq,
         cst2->nrows * sizeof(int));

  return cst1;
}

/* Adds cs2 to each element in cst1's list; cst2 should have a single element */
PlutoConstraints *pluto_constraints_add_to_each(PlutoConstraints *cst1,
                                                const PlutoConstraints *cst2) {
  assert(cst2->next == NULL);
  assert(cst1->ncols == cst2->ncols);

  if (cst1->nrows + cst2->nrows > cst1->alloc_nrows) {
    pluto_constraints_resize(cst1, cst1->nrows + cst2->nrows, cst1->ncols);
  } else {
    cst1->nrows = cst1->nrows + cst2->nrows;
  }

  for (int i = 0; i < (int)cst2->nrows; i++) {
    memcpy(cst1->val[(int)cst1->nrows - (int)cst2->nrows + i], cst2->val[i],
           cst1->ncols * sizeof(int64_t));
  }

  memcpy(&cst1->is_eq[cst1->nrows - cst2->nrows], cst2->is_eq,
         cst2->nrows * sizeof(int));

  if (cst1->next != NULL) {
    pluto_constraints_add_to_each(cst1->next, cst2);
  }
  return cst1;
}

/* Temporary structure to compare two rows */
struct row_info {
  int64_t *row;
  short is_eq;
  int ncols;
};

static int row_compar(const void *e1, const void *e2) {
  int i, ncols;
  int64_t *row1, *row2;
  struct row_info *u1, *u2;

  u1 = *(struct row_info **)e1;
  u2 = *(struct row_info **)e2;
  row1 = u1->row;
  row2 = u2->row;
  ncols = u1->ncols;
  assert(ncols == u2->ncols);

  for (i = 0; i < ncols; i++) {
    if (row1[i] != row2[i])
      break;
  }

  if (i == ncols) {
    /* Equal if both are inequalities or both are equalities;
     * otherwise, equalities will be sorted ahead of (>=0) inequalities */
    return u2->is_eq - u1->is_eq;
  } else if (row1[i] < row2[i]) {
    return -1;
  } else {
    return 1;
  }

  /* memcmp is inefficient compared to what's above when compiled with gcc */
  /* return memcmp(row1, row2, ncols*sizeof(int64_t)); */
}

/*
 * Eliminates duplicate constraints; the simplified constraints
 * are still at the same memory location but the number of constraints
 * will decrease
 */
void pluto_constraints_simplify(PlutoConstraints *const cst) {
  int64_t _gcd;

  if (cst->nrows == 0) {
    return;
  }

  PlutoConstraints *tmpcst =
      pluto_constraints_alloc(cst->nrows, cst->ncols, cst->context);
  tmpcst->nrows = 0;

  int *is_redun = (int *)malloc(sizeof(int) * cst->nrows);
  bzero(is_redun, cst->nrows * sizeof(int));

  /* Normalize cst - will help find redundancy */
  for (unsigned i = 0; i < cst->nrows; i++) {
    unsigned j;
    for (j = 0; j < cst->ncols; j++) {
      if (cst->val[i][j] != 0)
        break;
    }

    if (j < cst->ncols) {
      _gcd = PLABS(cst->val[i][j]);
      for (; j < cst->ncols; j++) {
        _gcd = gcd(PLABS(cst->val[i][j]), _gcd);
      }

      /* Normalize by gcd */
      for (j = 0; j < cst->ncols; j++) {
        cst->val[i][j] /= _gcd;
      }
    }
  }

  struct row_info **rows;
  rows = (struct row_info **)malloc(cst->nrows * sizeof(struct row_info *));
  for (unsigned i = 0; i < cst->nrows; i++) {
    rows[i] = (struct row_info *)malloc(sizeof(struct row_info));
    rows[i]->row = cst->val[i];
    rows[i]->is_eq = cst->is_eq[i];
    rows[i]->ncols = cst->ncols;
  }
  qsort(rows, cst->nrows, sizeof(struct row_info *), row_compar);

  for (unsigned i = 0; i < cst->nrows; i++) {
    cst->val[i] = rows[i]->row;
    cst->is_eq[i] = rows[i]->is_eq;
    free(rows[i]);
  }
  free(rows);

  is_redun[0] = 0;
  for (unsigned i = 1; i < cst->nrows; i++) {
    unsigned j;
    for (j = 0; j < cst->ncols; j++) {
      if (cst->val[i][j] != 0)
        break;
    }

    if (j == cst->ncols) {
      /* All zeros */
      is_redun[i] = 1;
      continue;
    }

    for (j = 0; j < cst->ncols; j++) {
      if (cst->val[i - 1][j] != cst->val[i][j])
        break;
    }

    if (j == cst->ncols && cst->is_eq[i - 1] == cst->is_eq[i]) {
      /* Same as cst(i-1) */
      is_redun[i] = 1;
    } else
      is_redun[i] = 0;
  }

  unsigned p = 0;
  for (unsigned i = 0; i < cst->nrows; i++) {
    if (!is_redun[i]) {
      memcpy(tmpcst->val[p], cst->val[i], cst->ncols * sizeof(int64_t));
      tmpcst->is_eq[p] = cst->is_eq[i];
      p++;
    }
  }
  tmpcst->nrows = p;
  tmpcst->ncols = cst->ncols;

  pluto_constraints_copy_single(cst, tmpcst);
  pluto_constraints_free(tmpcst);

  free(is_redun);

  if (cst->next != NULL)
    pluto_constraints_simplify(cst->next);
}

/*
 * Eliminates the pos^th variable, where pos has to be between 0 and
 * cst->ncols-2;
 * Remember that the last column is for the constant. The implementation does
 * not
 * have a redundancy check; it just  eliminates duplicates after gcd
 * normalization
 * cst will be resized if necessary
 */
void fourier_motzkin_eliminate(PlutoConstraints *cst, unsigned pos) {
  int64_t lb, ub, nb;
  int *bound;

  // At least one variable
  assert(cst->ncols >= 2);
  assert(pos <= (int)cst->ncols - 2);

  for (unsigned i = 0; i < cst->nrows; i++) {
    if (cst->is_eq[i]) {
      PlutoConstraints *tmpcst =
          pluto_constraints_to_pure_inequalities_single(cst);
      pluto_constraints_copy_single(cst, tmpcst);
      pluto_constraints_free(tmpcst);
      break;
    }
  }

  PlutoConstraints *newcst;

  unsigned i;
  for (i = 0; i < cst->nrows; i++) {
    if (cst->val[i][pos] != 0)
      break;
  }

  if (i == cst->nrows) {
    newcst = pluto_constraints_dup_single(cst);
    pluto_constraints_remove_dim(newcst, pos);
  } else {
    bound = (int *)malloc(cst->nrows * sizeof(int));

    lb = 0;
    ub = 0;
    nb = 0;

    /* Variable does appear */
    for (unsigned j = 0; j < cst->nrows; j++) {
      if (cst->val[j][pos] == 0) {
        bound[j] = NB;
        nb++;
      } else if (cst->val[j][pos] > 0) {
        bound[j] = LB;
        lb++;
      } else {
        bound[j] = UB;
        ub++;
      }
    }
    newcst = pluto_constraints_alloc(lb * ub + nb, cst->ncols, cst->context);
    pluto_constraints_copy_single(newcst, cst);
    pluto_constraints_remove_dim(newcst, pos);
    pluto_constraints_zero(newcst);
    newcst->nrows = 0;

    unsigned p = 0;
    for (unsigned j = 0; j < cst->nrows; j++) {
      if (bound[j] == UB) {
        for (unsigned k = 0; k < cst->nrows; k++) {
          if (bound[k] == LB) {
            unsigned q = 0;
            for (unsigned l = 0; l < cst->ncols; l++) {
              if (l != pos) {
                newcst->val[p][q] =
                    cst->val[j][l] * (lcm(cst->val[k][pos], -cst->val[j][pos]) /
                                      (-cst->val[j][pos])) +
                    cst->val[k][l] * (lcm(-cst->val[j][pos], cst->val[k][pos]) /
                                      cst->val[k][pos]);
                q++;
              }
            }
            p++;
          }
        }
      } else if (bound[j] == NB) {
        unsigned q = 0;
        for (unsigned l = 0; l < cst->ncols; l++) {
          if (l != pos) {
            newcst->val[p][q] = cst->val[j][l];
            q++;
          }
        }
        p++;
      }
    }
    assert(p <= lb * ub + nb);
    newcst->nrows = p;
    free(bound);
  }

  pluto_constraints_simplify(newcst);
  pluto_constraints_copy_single(cst, newcst);
  pluto_constraints_free(newcst);

  if (cst->next != NULL)
    fourier_motzkin_eliminate(cst->next, pos);
}

/* Copy constraints from src into dest; if dest does not have enough space,
 * resize it */
PlutoConstraints *pluto_constraints_copy(PlutoConstraints *dest,
                                         const PlutoConstraints *src) {
  pluto_constraints_copy_single(dest, src);

  if (src->next != NULL && dest->next != NULL) {
    pluto_constraints_copy(dest->next, src->next);
  }
  if (src->next == NULL && dest->next != NULL) {
    pluto_constraints_free(dest->next);
    dest->next = NULL;
  }
  if (src->next != NULL && dest->next == NULL) {
    dest->next = pluto_constraints_dup(src->next);
  }

  return dest;
}

/* Copy constraints from the first element of src into the first element
 * of dest; if dest does not have enough space, resize it. The _single
 * signifies only this src in the list is copied */
PlutoConstraints *pluto_constraints_copy_single(PlutoConstraints *dest,
                                                const PlutoConstraints *src) {
  if (src->nrows > dest->alloc_nrows || src->ncols > dest->alloc_ncols) {
    pluto_constraints_resize_single(dest, PLMAX(src->nrows, dest->alloc_nrows),
                                    PLMAX(src->ncols, dest->alloc_ncols));
  }

  /* Resize above may not be needed; need to still free all names;
   * they will be reassigned */
  if (dest->names) {
    for (int i = 0; i < (int)dest->ncols - 1; i++) {
      free(dest->names[i]);
      dest->names[i] = NULL;
    }
  }

  dest->nrows = src->nrows;
  dest->ncols = src->ncols;

  for (unsigned i = 0; i < dest->nrows; i++) {
    memcpy(dest->val[i], src->val[i], src->ncols * sizeof(int64_t));
  }

  memcpy(dest->is_eq, src->is_eq, dest->nrows * sizeof(int));
  if (src->names) {
    pluto_constraints_set_names(dest, src->names);
  } else {
    pluto_constraints_remove_names_single(dest);
  }

  return dest;
}

/* Duplicate constraints; returned constraints should be freed with
 * pluto_constraints_free */
PlutoConstraints *pluto_constraints_dup_single(const PlutoConstraints *src) {
  assert(src != NULL);

  PlutoConstraints *dup =
      pluto_constraints_alloc(src->alloc_nrows, src->alloc_ncols, src->context);

  pluto_constraints_copy_single(dup, src);

  return dup;
}

/* Duplicate constraints; returned constraints should be freed with
 * pluto_constraints_free */
PlutoConstraints *pluto_constraints_dup(const PlutoConstraints *src) {
  assert(src != NULL);

  PlutoConstraints *dup =
      pluto_constraints_alloc(src->alloc_nrows, src->alloc_ncols, src->context);

  pluto_constraints_copy(dup, src);

  return dup;
}

static void pluto_constraints_print_single(FILE *fp,
                                           const PlutoConstraints *cst,
                                           int set_num) {
  if (cst == NULL)
    return;

  fprintf(fp, "Set #%d\n", set_num + 1);
  fprintf(fp, "%d %d\n", cst->nrows, cst->ncols);

  for (unsigned i = 0; i < cst->nrows; i++) {
    for (unsigned j = 0; j < cst->ncols; j++) {
      fprintf(fp, "% 3ld ", cst->val[i][j]);
    }
    fprintf(fp, "\t %s 0\n", cst->is_eq[i] ? "==" : ">=");
  }

  if (cst->next) {
    pluto_constraints_print_single(fp, cst->next, set_num + 1);
  }
}

void pluto_constraints_print(FILE *fp, const PlutoConstraints *cst) {
  pluto_constraints_print_single(fp, cst, 0);
}

/* Print in polylib format */
static void
pluto_constraints_print_polylib_without_num(FILE *fp,
                                            const PlutoConstraints *const cst) {
  fprintf(fp, "%d %d\n", cst->nrows, cst->ncols + 1);

  for (unsigned i = 0; i < cst->nrows; i++) {
    fprintf(fp, "%s ", cst->is_eq[i] ? "0" : "1");
    for (unsigned j = 0; j < cst->ncols; j++) {
      fprintf(fp, "%ld ", cst->val[i][j]);
    }
    fprintf(fp, "\n");
  }

  if (cst->next) {
    pluto_constraints_print_polylib_without_num(fp, cst->next);
  }
}

void pluto_constraints_print_polylib(FILE *fp,
                                     const PlutoConstraints *const cst) {
  int num;

  num = pluto_constraints_num_in_list(cst);
  if (num >= 2)
    fprintf(fp, "%d\n", num);

  pluto_constraints_print_polylib_without_num(fp, cst);
}

/* Converts to polylib style matrix (first element of cst if it's a list) */
PlutoMatrix *pluto_constraints_to_matrix(const PlutoConstraints *cst) {
  PlutoMatrix *mat =
      pluto_matrix_alloc(cst->nrows, cst->ncols + 1, cst->context);

  for (unsigned i = 0; i < cst->nrows; i++) {
    mat->val[i][0] = (cst->is_eq[i]) ? 0 : 1;
    for (unsigned j = 0; j < cst->ncols; j++) {
      mat->val[i][j + 1] = cst->val[i][j];
    }
  }

  return mat;
}

/* Create pluto_constraints from polylib-style matrix  */
PlutoConstraints *pluto_constraints_from_mixed_matrix(const PlutoMatrix *mat,
                                                      int *is_eq) {
  PlutoConstraints *cst;

  cst = pluto_constraints_alloc(mat->nrows, mat->ncols, mat->context);

  cst->nrows = mat->nrows;

  for (unsigned i = 0; i < cst->nrows; i++) {
    cst->is_eq[i] = is_eq[i];
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = mat->val[i][j];
    }
  }

  return cst;
}

/* Create pluto_constraints from polylib-style matrix  */
PlutoConstraints *pluto_constraints_from_matrix(const PlutoMatrix *mat) {
  PlutoConstraints *cst;

  cst = pluto_constraints_alloc(mat->nrows, mat->ncols - 1, mat->context);

  cst->nrows = mat->nrows;

  for (unsigned i = 0; i < cst->nrows; i++) {
    cst->is_eq[i] = (mat->val[i][0] == 0);
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = mat->val[i][j + 1];
    }
  }

  return cst;
}

/* Create pluto_constraints from matrix  */
PlutoConstraints *pluto_constraints_from_inequalities(const PlutoMatrix *mat) {
  PlutoConstraints *cst;

  cst = pluto_constraints_alloc(mat->nrows, mat->ncols, mat->context);
  cst->nrows = mat->nrows;

  for (unsigned i = 0; i < cst->nrows; i++) {
    cst->is_eq[i] = 0;
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = mat->val[i][j];
    }
  }
  return cst;
}

/* Read constraints in polylib format from a file */
PlutoConstraints *pluto_constraints_read(FILE *fp, PlutoContext *context) {
  int ineq, nrows, ncols, retval, num;

  fscanf(fp, "%d", &nrows);
  retval = fscanf(fp, "%d", &ncols);
  ncols--;

  if (retval == EOF || retval == 0) {
    return NULL;
  }

  PlutoConstraints *cst = pluto_constraints_alloc(nrows, ncols, context);
  cst->nrows = nrows;

  for (unsigned i = 0; i < cst->nrows; i++) {
    fscanf(fp, "%d ", &ineq);
    cst->is_eq[i] = ineq ^ 1;
    for (unsigned j = 0; j < cst->ncols; j++) {
      fscanf(fp, "%ld", &cst->val[i][j]);
    }
  }

  if (fscanf(fp, "%d", &num) != EOF) {
    if (num >= 1)
      cst->names = (char **)malloc(num * sizeof(char *));
    for (int i = 0; i < num; i++) {
      cst->names[i] = (char *)malloc(6);
      fscanf(fp, "%s", cst->names[i]);
    }
  }

  return cst;
}

void pluto_constraints_compact_print_single(FILE *fp,
                                            const PlutoConstraints *cst,
                                            int set_num) {
  int i, j, nrows, ncols;

  if (cst == NULL) {
    return;
  }

  nrows = cst->nrows;
  ncols = cst->ncols;

  printf("Set #%d\n", set_num + 1);

  if (nrows == 0) {
    printf("Universal polyhedron -- No constraints (%d dims)!\n",
           cst->ncols - 1);
    return;
  }

  fprintf(fp, "[%d dims; %d constraints]\n", cst->ncols - 1, cst->nrows);

  for (i = 0; i < nrows; i++) {
    int first = 1;
    for (j = 0; j < ncols; j++) {
      if (j == ncols - 1) {
        /* constant */
        if (cst->val[i][j] == 0 && !first)
          fprintf(fp, " ");
        else
          fprintf(fp, "%s%ld ", (cst->val[i][j] >= 0 && !first) ? "+" : "",
                  cst->val[i][j]);
      } else {
        char var[6];
        var[5] = '\0';
        if (cst->names) {
          assert(cst->names[j]);
          strncpy(var, cst->names[j], 5);
        } else
          sprintf(var, "c_%c", 'i' + j);

        if (cst->val[i][j] == 1) {
          fprintf(fp, "%s%s", first ? "" : "+", var);
          first = 0;
        } else if (cst->val[i][j] == -1) {
          fprintf(fp, "-%s", var);
          first = 0;
        } else if (cst->val[i][j] >= 2) {
          fprintf(fp, "%s%ld%s", first ? "" : "+", cst->val[i][j], var);
          first = 0;
        } else if (cst->val[i][j] <= -2) {
          fprintf(fp, "%ld%s", cst->val[i][j], var);
          first = 0;
        }
      }
    }
    fprintf(fp, "%s 0\n", cst->is_eq[i] ? "=" : ">=");
  }
  fprintf(fp, "\n");

  if (cst->next != NULL) {
    pluto_constraints_compact_print_single(fp, cst->next, set_num + 1);
  }
}

void pluto_constraints_compact_print(FILE *fp, const PlutoConstraints *cst) {
  assert(cst != NULL);
  pluto_constraints_compact_print_single(fp, cst, 0);
}

void pluto_constraints_pretty_print(FILE *fp, const PlutoConstraints *cst) {
  int i, j;

  int nrows = cst->nrows;
  int ncols = cst->ncols;

  assert(cst->next == NULL);

  if (nrows == 0) {
    printf("Universal polyhedron -- no constraints (%d dims)!\n",
           cst->ncols - 1);
    return;
  }

  fprintf(fp, "[%d dims; %d constraints]\n", cst->ncols - 1, cst->nrows);

  for (i = 0; i < nrows; i++) {
    /* Is it first non-zero entry */
    int first = 1;
    for (j = 0; j < ncols; j++) {
      if (j == ncols - 1) {
        /* constant */
        if (cst->val[i][j] == 0 && !first)
          fprintf(fp, "     ");
        else
          fprintf(fp, "%s%ld ", (cst->val[i][j] >= 0) ? "+" : "",
                  cst->val[i][j]);
      } else {
        char var[6];
        var[5] = '\0';
        if (cst->names)
          strncpy(var, cst->names[j], 5);
        else
          sprintf(var, "c_%c", 'i' + j);

        if (cst->val[i][j] == 1) {
          fprintf(fp, "+%s ", var);
          first = 0;
        } else if (cst->val[i][j] == -1) {
          fprintf(fp, "-%s ", var);
          first = 0;
        } else if (cst->val[i][j] >= 2) {
          fprintf(fp, "+%ld%s ", cst->val[i][j], var);
          first = 0;
        } else if (cst->val[i][j] <= -2) {
          fprintf(fp, "%ld%s ", cst->val[i][j], var);
          first = 0;
        } else
          fprintf(fp, "      ");

        if (cst->val[i][j] != 0)
          first = 0;
      }
    }
    fprintf(fp, "%s 0\n", cst->is_eq[i] ? "=" : ">=");
  }
  fprintf(fp, "\n");
}

void pluto_constraints_cplex_print(FILE *fp, const PlutoConstraints *cst) {
  int i;

  int nrows = cst->nrows;
  int ncols = cst->ncols;

  assert(cst->next == NULL);

  if (nrows == 0) {
    printf("No constraints!\n");
  }

  for (i = 0; i < nrows; i++) {
    int first = 1;
    for (unsigned char j = 0; j < ncols; j++) {
      if (j == ncols - 1) {
        /* constant */
        /* Not supported in CPLEX format */
        assert(!first || cst->val[i][j] >= 0);
        if (!first)
          fprintf(fp, "%s %ld\n", cst->is_eq[i] ? "=" : ">=", -cst->val[i][j]);
      } else {
        char var[6];
        if (cst->names)
          strncpy(var, cst->names[j], 5);
        else
          snprintf(var, 6, "c_%d", j);

        if (cst->val[i][j] == 1) {
          fprintf(fp, "+%s ", var);
          first = 0;
        } else if (cst->val[i][j] == -1) {
          fprintf(fp, "-%s ", var);
          first = 0;
        } else if (cst->val[i][j] >= 2) {
          fprintf(fp, "+%ld%s ", cst->val[i][j], var);
          first = 0;
        } else if (cst->val[i][j] <= -2) {
          fprintf(fp, "%ld%s ", cst->val[i][j], var);
          first = 0;
        } else
          fprintf(fp, "    ");
      }
    }
  }
  fprintf(fp, "\n");
}

/* Convert Pluto constraints into PIP format (first column is
 * 0/1 based on equality/inequality */
PlutoMatrix *pluto_constraints_to_pip_matrix(const PlutoConstraints *cst,
                                             PlutoMatrix *pmat) {
  assert(pmat);

  for (unsigned i = 0; i < cst->nrows; i++) {
    /* Equality or Inequality */
    pmat->val[i][0] = cst->is_eq[i] ? 0 : 1;
    for (unsigned j = 1; j < cst->ncols + 1; j++) {
      pmat->val[i][j] = cst->val[i][j - 1];
    }
  }
  pmat->nrows = cst->nrows;
  pmat->ncols = cst->ncols + 1;

  return pmat;
}

/* Use PIP to solve these constraints (solves for the first element
 * if it's a list of constraints) */
int64_t *pluto_constraints_lexmin_pip(const PlutoConstraints *cst, int negvar) {
  int bignum, i;
  PipMatrix *domain, *pip_context;
  PipQuast *solution;
  PipOptions *pipOptions;
  PipList *listPtr;
  int64_t *sol;
  PlutoMatrix *pipmat;
  PlutoContext *context = cst->context;

  IF_DEBUG2(printf("[pluto] pluto_constraints_lexmin_pip (%d variables, %d "
                   "constraints)\n",
                   cst->ncols - 1, cst->nrows););

  pipmat = pluto_matrix_alloc(cst->nrows, cst->ncols + 1, cst->context);

  /* Convert constraints to PIP format */
  pluto_constraints_to_pip_matrix(cst, pipmat);

  /* First column says whether it's an inequality and
   * the last is for the constant; so we have ncols-2 variables */
  pip_context = NULL;
  bignum = -1;

  domain = pip_matrix_populate(pipmat->val, pipmat->nrows, pipmat->ncols);

  pipOptions = pip_options_init();

  if (negvar == 1) {
    pipOptions->Urs_parms = 1;
    pipOptions->Urs_unknowns = 1;
  }

  solution = pip_solve(domain, pip_context, bignum, pipOptions);

  assert(solution->condition == NULL);

  listPtr = solution->list;

  sol = NULL;
  if (listPtr != NULL) {
    sol = (int64_t *)malloc((pipmat->ncols - 2) * sizeof(int64_t));

    for (i = 0; i < pipmat->ncols - 2; i++) {
      /* This is just a lexmin and not a parametrix lexmin and so each
       * vector in the list is actually a constant */
#ifdef PIP_WIDTH_MP
      sol[i] = mpz_get_si(*listPtr->vector->the_vector);
#else
      sol[i] = (int64_t)*listPtr->vector->the_vector;
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

/// Solve these constraints for lexmin. solution.
//  TODO: GLPK-based path doesn't exist here.
int64_t *pluto_constraints_lexmin(const PlutoConstraints *cst, int negvar) {
  if (cst->context->options->islsolve)
    return pluto_constraints_lexmin_isl(cst, negvar);
  // Use PIP.
  return pluto_constraints_lexmin_pip(cst, negvar);
}

#ifdef GLPK
/* Constructs constraints for glpk problem in-memory.  Assumes that there
 * are no rows or cols in glp_prob lp.*/
void set_glpk_constraints_from_pluto_constraints(glp_prob *lp,
                                                 const PlutoConstraints *cst) {
  int i, j, k, nrows, ncols, non_zero;
  int *row, *col;
  double *coeff;

  nrows = cst->nrows;
  ncols = cst->ncols;
  non_zero = pluto_constraints_get_num_non_zero_coeffs(cst);

  assert(lp != NULL);

  glp_add_rows(lp, cst->nrows);
  glp_add_cols(lp, cst->ncols - 1);

  /* These are indexed from 1 */
  row = (int *)malloc((non_zero + 1) * sizeof(int));
  col = (int *)malloc((non_zero + 1) * sizeof(int));
  coeff = (double *)malloc((non_zero + 1) * sizeof(double));

  k = 1;
  for (i = 0; i < nrows; i++) {
    if (cst->is_eq[i]) {
      /* For equality constraints the bounds are fixed */
      glp_set_row_bnds(lp, i + 1, GLP_FX, -cst->val[i][cst->ncols - 1], 0.0);
    } else {
      /* Add a lower bound on each row of the inequality */
      glp_set_row_bnds(lp, i + 1, GLP_LO, -cst->val[i][cst->ncols - 1], 0.0);
    }

    for (j = 0; j < ncols - 1; j++) {
      if (cst->val[i][j] != 0) {
        row[k] = i + 1;
        col[k] = j + 1;
        coeff[k] = (double)cst->val[i][j];
        k++;
      }
    }
  }
  glp_load_matrix(lp, non_zero, row, col, coeff);
  free(row);
  free(col);
  free(coeff);
}

/* Retrives ilp solution from the input glpk problem.
 * Assumes that the optimal solution exists and has been found*/
double *get_lp_solution_from_glpk_problem(glp_prob *lp) {
  int ncols = glp_get_num_cols(lp);
  double *sol = (double *)malloc(sizeof(double) * ncols);
  for (int i = 0; i < glp_get_num_cols(lp); i++) {
    double x = glp_mip_col_val(lp, i + 1);
    sol[i] = x;
  }
  return sol;
}

/* Retrives ilp solution from the input glpk problem.
 * Assumes that the optimal solution exists and has been found*/
int64_t *get_ilp_solution_from_glpk_problem(glp_prob *lp,
                                            PlutoContext *context) {
  int i, ncols;
  int64_t *sol;

  ncols = glp_get_num_cols(lp);
  sol = (int64_t *)malloc(sizeof(int64_t) * ncols);
  for (i = 0; i < glp_get_num_cols(lp); i++) {
    double x = glp_mip_col_val(lp, i + 1);
    IF_DEBUG(printf("c%d = %ld, ", i, (int64_t)round(x)););
    sol[i] = (int64_t)round(x);
  }
  IF_DEBUG(printf("\n"););
  return sol;
}

/* Set glpk problem parameters.
 * Checks feasibility of the LP problem using simplex */
void set_glpk_problem_params(glp_prob *lp, PlutoContext *context) {
  PlutoOptions *options = context->options;
  if (!options->moredebug) {
    glp_term_out(GLP_OFF);
  }

  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.presolve = GLP_ON;

  parm.msg_lev = GLP_MSG_OFF;
  IF_MORE_DEBUG(parm.msg_lev = GLP_MSG_ON;);
  IF_MORE_DEBUG(parm.msg_lev = GLP_MSG_ALL;);

  glp_scale_prob(lp, GLP_SF_AUTO);
  glp_adv_basis(lp, 0);
  glp_simplex(lp, &parm);
}

void find_optimal_solution_glpk(glp_prob *lp, double tol,
                                PlutoContext *context) {
  glp_iocp iocp;
  glp_init_iocp(&iocp);
  /* The default is 1e-5; one may need to reduce it even further
   * depending on how large a coefficient we might see */
  iocp.tol_int = tol;
  IF_MORE_DEBUG(printf("Setting GLPK integer tolerance to %e\n", iocp.tol_int));

  iocp.msg_lev = GLP_MSG_OFF;
  IF_MORE_DEBUG(iocp.msg_lev = GLP_MSG_ON;);
  IF_MORE_DEBUG(iocp.msg_lev = GLP_MSG_ALL;);

  /* Find optimal solution */
  glp_intopt(lp, &iocp);
}

/* Returns 0 if a solution was found else returns 1 */
int pluto_constraints_solve_glpk(glp_prob *lp, PlutoContext *context) {
  set_glpk_problem_params(lp, context);
  int lp_status = glp_get_status(lp);

  if (lp_status == GLP_INFEAS || lp_status == GLP_UNDEF) {
    glp_delete_prob(lp);
    return 1;
  }

  find_optimal_solution_glpk(lp, 1e-7, context);

  int ilp_status = glp_mip_status(lp);

  /* ILP optimization can return any of these exit codes representing cases
   * where the problem is infeasible (or not feasible) undefined or unbounded.
   * All these cases are undesireable. If the solution is read when the lp
   * problem returns with any of these exit codes, then a zero solution will be
   * obtained */
  if (ilp_status == GLP_INFEAS || ilp_status == GLP_UNDEF ||
      ilp_status == GLP_NOFEAS || ilp_status == GLP_UNBND) {
    glp_delete_prob(lp);
    return 1;
  }

  double z = glp_mip_obj_val(lp);
  IF_DEBUG(printf("Objective value: z = %lf\n", z););

  return 0;
}

double *pluto_mip_scale_solutions_glpk(glp_prob *ilp, PlutoContext *context) {
  double *scale_sols;

  set_glpk_problem_params(ilp, context);

  int lp_status = glp_get_status(ilp);

  if (lp_status == GLP_INFEAS || lp_status == GLP_UNDEF) {
    glp_delete_prob(ilp);
    return NULL;
  }

  find_optimal_solution_glpk(ilp, 1e-2, context);

  int ilp_status = glp_mip_status(ilp);

  if (ilp_status == GLP_NOFEAS) {
    glp_delete_prob(ilp);
    return NULL;
  }

  double z = glp_mip_obj_val(ilp);
  IF_DEBUG(printf("z = %lf\n", z););

  scale_sols = get_lp_solution_from_glpk_problem(ilp);
  return scale_sols;
}

glp_prob *get_scaling_lp_glpk(double *fpsol, int num_sols, double **val,
                              int **index, int npar, int num_ccs,
                              PlutoContext *context) {
  PlutoOptions *options = context->options;
  int i, num_rows, num_sols_to_be_scaled, col_num;
  glp_prob *lp;

  /* first npar+1 elements in the LP solution (corresponding to u and w) need
   * not be scaled.
   * The first num_ccs coeffs correspond to the scaling factors.
   * The remaining num_elements_to_be_scaled need to be scaled up if
   * they are higher than the INT_TOL limit set. */
  num_sols_to_be_scaled = num_sols - npar - 1;

  IF_DEBUG(printf("[Pluto]: Number of connected components: %d\n", num_ccs););
  IF_DEBUG(printf("[Pluto]: Number of solutions to be scaled: %d\n",
                  num_sols_to_be_scaled););

  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);
  glp_add_cols(lp, num_ccs + num_sols_to_be_scaled);

  num_rows = 0;

  for (i = npar + 1; i < num_sols; i++) {
    col_num = i - (npar + 1) + num_ccs + 1;
    if (fabs(fpsol[i]) > 1e-7) {
      glp_add_rows(lp, 1);
      num_rows++;
      val[i - npar - 1][1] = fpsol[i];
      glp_set_row_bnds(lp, num_rows, GLP_FX, 0.0, 0.0);
      glp_set_col_bnds(lp, col_num, GLP_LO, 1.0, 0.0);
      glp_set_col_kind(lp, col_num, GLP_IV);
      glp_set_mat_row(lp, num_rows, 2, index[i - npar - 1], val[i - npar - 1]);
    } else {
      glp_set_col_bnds(lp, col_num, GLP_FX, 0.0, 0.0);
    }
  }

  /* The lower bound of the scaling factor for each CC has to be one */
  for (i = 0; i < num_ccs; i++) {
    glp_set_col_bnds(lp, i + 1, GLP_LO, 1.0, 0.0);
    glp_set_obj_coef(lp, i + 1, 1.0);
  }
  if (options->debug) {
    glp_write_lp(lp, NULL, "pluto-scaling-mip.lp");
  }
  return lp;
}

void set_glpk_objective_from_pluto_matrix(glp_prob *lp, PlutoMatrix *obj) {
  int j;
  for (j = 0; j < obj->ncols; j++) {
    glp_set_obj_coef(lp, j + 1, (double)obj->val[0][j]);
  }
}

/* Construct ILP in cplex format. The last four parameters are used for
 * scaling solutions to integers */
int64_t *pluto_prog_constraints_lexmin_glpk(const PlutoConstraints *cst,
                                            PlutoMatrix *obj, double **val,
                                            int **index, int npar,
                                            int num_ccs) {
  int i, j, is_unsat, num_sols;
  int64_t *sol;
  double *fpsol, *scale_sols;
  PlutoContext *context = cst->context;
  PlutoOptions *options = context->options;

  IF_DEBUG(printf("[pluto] pluto_prog_constraints_lexmin_glpk (%d variables, "
                  "%d constraints)\n",
                  cst->ncols - 1, cst->nrows););

  glp_prob *lp;
  lp = glp_create_prob();

  glp_set_obj_dir(lp, GLP_MIN);

  set_glpk_constraints_from_pluto_constraints(lp, cst);

  set_glpk_objective_from_pluto_matrix(lp, obj);

  for (i = 0; i < cst->ncols - 1; i++) {
    glp_set_col_bnds(lp, i + 1, GLP_LO, 0.0, 0.0);
  }

  if (!options->lp) {
    for (i = 0; i < cst->ncols - 1; i++) {
      glp_set_col_kind(lp, i + 1, GLP_IV);
    }
  }
  IF_DEBUG(glp_write_lp(lp, NULL, "pluto-debug-glpk.lp"););

  is_unsat = pluto_constraints_solve_glpk(lp, context);

  if (is_unsat) {
    return NULL;
  }

  if (options->lp || options->dfp) {

    fpsol = get_lp_solution_from_glpk_problem(lp);
    num_sols = glp_get_num_cols(lp);
    glp_delete_prob(lp);

    lp = get_scaling_lp_glpk(fpsol, num_sols, val, index, npar, num_ccs,
                             cst->context);

    scale_sols = pluto_mip_scale_solutions_glpk(lp, context);
    int64_t *sol = (int64_t *)malloc(sizeof(int64_t) * (num_sols));

    /* Ideally u and w have to be set to be computed on a per CC basis. Since it
     * is not
     * used further down the tool chain, it is set to the maximum scaling
     * factor.*/
    int64_t max_scale_factor = (int64_t)round(glp_mip_col_val(lp, 1));
    IF_MORE_DEBUG(printf("Scaling Factor for CC 1:%ld\n", max_scale_factor););
    for (j = 0; j < num_ccs; j++) {
      IF_MORE_DEBUG(printf("Scaling Factor for CC %d:%ld\n", j + 1,
                           (int64_t)round(glp_mip_col_val(lp, j + 1))););
      if (scale_sols[j] > max_scale_factor) {
        max_scale_factor = (int64_t)round(scale_sols[j]);
      }
    }

    for (int j = 0; j < npar + 1; j++) {
      sol[j] = max_scale_factor * (int64_t)round(fpsol[j]);
    }

    int col_iter = num_ccs;
    for (j = npar + 1; j < cst->ncols - 1; j++) {
      double x = scale_sols[col_iter++];
      IF_DEBUG(printf("c%d = %ld, ", j, (int64_t)round(x)););
      sol[j] = (int64_t)round(x);
    }
    IF_DEBUG(printf("\n"););

    glp_delete_prob(lp);
    free(fpsol);
    free(scale_sols);
    return sol;
  } else {
    sol = get_ilp_solution_from_glpk_problem(lp, cst->context);
    glp_delete_prob(lp);
    return sol;
  }
}

double *pluto_fcg_constraints_lexmin_glpk(const PlutoConstraints *cst,
                                          PlutoMatrix *obj) {
  PlutoContext *context = cst->context;
  PlutoOptions *options = cst->context->options;
  glp_prob *lp;
  double *sol;
  int i, is_unsat;

  lp = glp_create_prob();

  glp_set_obj_dir(lp, GLP_MIN);

  set_glpk_constraints_from_pluto_constraints(lp, cst);

  set_glpk_objective_from_pluto_matrix(lp, obj);

  for (i = 0; i < cst->ncols - 1; i++) {
    glp_set_col_bnds(lp, i + 1, GLP_LO, 0.0, 0.0);
  }

  if (options->ilp) {
    for (i = 0; i < cst->ncols - 1; i++) {
      glp_set_col_kind(lp, i + 1, GLP_IV);
    }
  }
  IF_DEBUG(glp_write_lp(lp, NULL, "pluto-pairwise-constraints-glpk.lp"););

  is_unsat = pluto_constraints_solve_glpk(lp, cst->context);

  if (is_unsat) {
    return NULL;
  } else {
    sol = get_lp_solution_from_glpk_problem(lp);
    glp_delete_prob(lp);
    return sol;
  }
}

#endif

/* All equalities */
PlutoConstraints *pluto_constraints_from_equalities(const PlutoMatrix *mat) {
  PlutoConstraints *cst =
      pluto_constraints_alloc(mat->nrows, mat->ncols, mat->context);

  for (unsigned i = 0; i < mat->nrows; i++) {
    cst->is_eq[i] = 1;
    for (unsigned j = 0; j < mat->ncols; j++) {
      cst->val[i][j] = mat->val[i][j];
    }
  }
  cst->nrows = mat->nrows;

  return cst;
}

/* Add an inequality (>= 0); initialize it to all zero; will be added
 * as the last row */
void pluto_constraints_add_inequality(PlutoConstraints *cst) {
  if (cst->nrows == cst->alloc_nrows) {
    pluto_constraints_resize_single(cst, cst->nrows + 1, cst->ncols);
  } else {
    cst->nrows++;
  }

  for (unsigned j = 0; j < cst->ncols; j++) {
    cst->val[cst->nrows - 1][j] = 0;
  }
  cst->is_eq[cst->nrows - 1] = 0;

  if (cst->next) {
    pluto_constraints_add_inequality(cst->next);
  }
}

/* Add an equality (== 0); INITIALIZE IT TO ALL ZERO; will be added
 * as the last row */
void pluto_constraints_add_equality(PlutoConstraints *cst) {
  assert(cst->nrows <= cst->alloc_nrows);

  if (cst->nrows == cst->alloc_nrows) {
    pluto_constraints_resize_single(cst, cst->nrows + 1, cst->ncols);
  } else {
    cst->nrows++;
  }

  for (unsigned j = 0; j < cst->ncols; j++) {
    cst->val[cst->nrows - 1][j] = 0;
  }
  cst->is_eq[cst->nrows - 1] = 1;

  if (cst->next != NULL) {
    pluto_constraints_add_equality(cst->next);
  }
}

/* Add a constraint; initialize it to all zero */
void pluto_constraints_add_constraint(PlutoConstraints *cst, int is_eq) {
  if (is_eq)
    pluto_constraints_add_equality(cst);
  else
    pluto_constraints_add_inequality(cst);
}

/* Remove a row; pos is 0-indexed */
void pluto_constraints_remove_row(PlutoConstraints *cst, unsigned pos) {
  assert(pos <= (int)cst->nrows - 1);

  for (int i = pos; i < (int)cst->nrows - 1; i++) {
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = cst->val[i + 1][j];
    }
    cst->is_eq[i] = cst->is_eq[i + 1];
  }
  cst->nrows--;
}

/* Remove a variable */
void pluto_constraints_remove_dim(PlutoConstraints *cst, int pos) {
  assert(pos >= 0 && pos <= (int)cst->ncols - 2);

  if (cst->names)
    free(cst->names[pos]);

  for (unsigned i = 0; i < cst->nrows; i++) {
    for (int j = pos; j < (int)cst->ncols - 1; j++) {
      cst->val[i][j] = cst->val[i][j + 1];
    }
  }

  for (int j = pos; j < (int)cst->ncols - 1; j++) {
    if (cst->names && j < (int)cst->ncols - 2 && cst->names[j + 1]) {
      cst->names[j] = cst->names[j + 1];
    }
  }
  if (cst->names) {
    cst->names =
        (char **)realloc(cst->names, (cst->ncols - 2) * sizeof(char *));
  }

  cst->ncols--;

  if (cst->next != NULL)
    pluto_constraints_remove_dim(cst->next, pos);
}

/* Blank 'pos' th row */
void pluto_constraints_zero_row(PlutoConstraints *cst, int pos) {
  assert(cst->next == NULL);
  assert(pos >= 0 && pos <= cst->nrows - 1);

  for (unsigned j = 0; j < cst->ncols; j++) {
    cst->val[pos][j] = 0;
  }
}

/* Normalize row by its gcd */
void pluto_constraints_normalize_row(PlutoConstraints *cst, int pos) {
  int i, j, k;

  assert(cst->next == NULL);
  assert(pos >= 0 && pos <= cst->nrows - 1);

  /* Normalize cst first */
  for (i = 0; i < cst->nrows; i++) {
    if (cst->val[i][0] == 0)
      continue;
    int rowgcd = abs(cst->val[i][0]);
    for (j = 1; j < cst->ncols; j++) {
      if (cst->val[i][j] == 0)
        break;
      rowgcd = gcd(rowgcd, abs(cst->val[i][j]));
    }
    if (i == cst->nrows) {
      if (rowgcd >= 2) {
        for (k = 0; k < cst->ncols; k++) {
          cst->val[i][k] /= rowgcd;
        }
      }
    }
  }
}

/* Add a variable; resize if necessary; initialize to zero */
void pluto_constraints_add_dim(PlutoConstraints *cst, int pos,
                               const char *name) {
  int i, j;

  /* Has to be a new variable */
  assert(pos >= 0 && pos <= cst->ncols - 1);

  if (cst->ncols == cst->alloc_ncols) {
    pluto_constraints_resize_single(cst, cst->nrows, cst->ncols + 1);
  } else {
    cst->ncols++;
  }
  if (cst->names)
    cst->names =
        (char **)realloc(cst->names, (cst->ncols - 1) * sizeof(char *));

  for (j = cst->ncols - 2; j >= pos; j--) {
    for (i = 0; i < cst->nrows; i++) {
      cst->val[i][j + 1] = cst->val[i][j];
    }
    if (cst->names && j <= cst->ncols - 3) {
      cst->names[j + 1] = cst->names[j];
    }
  }

  /* Initialize to zero */
  for (i = 0; i < cst->nrows; i++) {
    cst->val[i][pos] = 0;
  }
  if (cst->names)
    cst->names[pos] = strdup(name ? name : "_u");

  if (cst->next != NULL) {
    pluto_constraints_add_dim(cst->next, pos, name);
  }
}

/// Adds dimensions to the constraints at the beginning. The corresponding
/// coefficients of the new variables are nitialized to zero. This is a very
/// inefficient implemetation of the routine. The constraints can be resized at
/// one go and the initialization of these dimensions can be done to zero
/// instead of calling pluto_constraints_add_dim num_dims times
void pluto_constraints_add_leading_dims(PlutoConstraints *cst, int num_dims) {
  for (int i = 0; i < num_dims; i++) {
    pluto_constraints_add_dim(cst, 0, NULL);
  }
}
/* Add a lower bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_lb(PlutoConstraints *cst, int varnum, int64_t lb) {
  assert(varnum >= 0 && varnum <= cst->ncols - 2);

  pluto_constraints_add_inequality(cst);

  while (cst != NULL) {
    cst->val[cst->nrows - 1][varnum] = 1;
    cst->val[cst->nrows - 1][cst->ncols - 1] = -lb;
    cst = cst->next;
  }
}

/* Add an upper bound for 'varnum' variable: varnum: 0-indexed */
void pluto_constraints_add_ub(PlutoConstraints *cst, int varnum, int64_t ub) {
  assert(varnum >= 0 && varnum <= cst->ncols - 2);

  pluto_constraints_add_inequality(cst);

  while (cst != NULL) {
    cst->val[cst->nrows - 1][varnum] = -1;
    cst->val[cst->nrows - 1][cst->ncols - 1] = ub;
    cst = cst->next;
  }
}

/* Set a value for a variable: varnum: 0-indexed */
void pluto_constraints_set_var(PlutoConstraints *cst, int varnum, int64_t val) {
  assert(varnum >= 0 && varnum <= cst->ncols - 2);

  pluto_constraints_add_equality(cst);

  while (cst != NULL) {
    cst->val[cst->nrows - 1][varnum] = 1;
    cst->val[cst->nrows - 1][cst->ncols - 1] = -val;
    cst = cst->next;
  }
}

/* Populate a PIP matrix */
PipMatrix *pip_matrix_populate(int64_t **cst, int nrows, int ncols) {
  int i, j;
  PipMatrix *matrix;
  Entier *p;

  matrix = pip_matrix_alloc(nrows, ncols);

  p = matrix->p_Init;
  for (i = 0; i < matrix->NbRows; i++) {
    for (j = 0; j < matrix->NbColumns; j++) {
      *(p++) = cst[i][j];
    }
  }
  return matrix;
}

PlutoConstraints *pluto_constraints_select_row(const PlutoConstraints *cst,
                                               int pos) {
  int j;

  PlutoConstraints *row = pluto_constraints_alloc(1, cst->ncols, cst->context);
  row->is_eq[0] = cst->is_eq[pos];
  for (j = 0; j < cst->ncols; j++) {
    row->val[0][j] = cst->val[pos][j];
  }
  row->nrows = 1;
  return row;
}

/*
 * Negate a single row in cst.
 */
void pluto_constraints_negate_row(PlutoConstraints *cst, int pos) {
  assert(cst->next == NULL);

  int j;

  for (j = 0; j < cst->ncols; j++) {
    cst->val[pos][j] = -cst->val[pos][j];
  }
}

/*
 * Negation of a single constraint in cst.
 */
void pluto_constraints_negate_constraint(PlutoConstraints *cst, int pos) {
  int j;

  assert(cst->next == NULL);

  for (j = 0; j < cst->ncols - 1; j++) {
    cst->val[pos][j] = -cst->val[pos][j];
  }
  cst->val[pos][cst->ncols - 1]--;
}

/* Convert everything to >= 0 form; for all equalities, add a single extra
 * constraint that puts the sum of those to be <= 0; it is only performed on a
 * single PlutoConstraints in a list */
PlutoConstraints *
pluto_constraints_to_pure_inequalities_single(const PlutoConstraints *cst) {
  int i, j;

  PlutoConstraints *ineq = pluto_constraints_dup_single(cst);

  for (i = 0; i < ineq->nrows; i++) {
    ineq->is_eq[i] = 0;
  }

  int has_eq = 0;
  for (i = 0; i < cst->nrows; i++) {
    if (cst->is_eq[i]) {
      has_eq = 1;
      break;
    }
  }

  if (has_eq) {
    /* Add a constraint to make sum of all equalities <= 0 */
    PlutoConstraints *neg_eq =
        pluto_constraints_alloc(1, cst->ncols, cst->context);

    for (i = 0; i < cst->nrows; i++) {
      if (cst->is_eq[i]) {
        for (j = 0; j < cst->ncols; j++) {
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
void pluto_constraints_project_out_single(PlutoConstraints *cst, int start,
                                          int num) {
  int i, end;

  assert(num >= 0);

  if (num == 0)
    return;

  end = start + num - 1;

  assert(start >= 0 && end <= cst->ncols - 2);

  PlutoMatrix *func =
      pluto_matrix_alloc(cst->ncols - num, cst->ncols, cst->context);
  pluto_matrix_set(func, 0);
  for (i = 0; i < start; i++) {
    func->val[i][i] = 1;
  }
  for (i = end + 1; i < cst->ncols; i++) {
    func->val[i - num][i] = 1;
  }

  PlutoConstraints *img = pluto_constraints_image(cst, func);
  pluto_constraints_copy_single(cst, img);

  pluto_constraints_free(img);
  pluto_matrix_free(func);
}

/* start: 0-indexed */
void pluto_constraints_project_out(PlutoConstraints *cst, int start, int num) {
  pluto_constraints_project_out_single(cst, start, num);

  if (cst->next != NULL) {
    pluto_constraints_project_out(cst->next, start, num);
  }
}

void pluto_constraints_interchange_cols(PlutoConstraints *cst, int col1,
                                        int col2) {
  if (col1 == col2)
    return;

  for (unsigned r = 0; r < cst->nrows; r++) {
    int64_t tmp = cst->val[r][col1];
    cst->val[r][col1] = cst->val[r][col2];
    cst->val[r][col2] = tmp;
  }

  if (cst->next != NULL) {
    pluto_constraints_interchange_cols(cst->next, col1, col2);
  }
}

void check_redundancy(PlutoConstraints *cst) {
  PlutoContext *context = cst->context;
  PlutoConstraints *check =
      pluto_constraints_alloc(cst->nrows, cst->ncols, context);
  int count;

  assert(cst->next == NULL);

  count = 0;

  for (int i = 0; i < cst->nrows; i++) {
    pluto_constraints_copy(check, cst);
    PlutoConstraints *row = pluto_constraints_select_row(cst, i);
    pluto_constraints_remove_row(check, i);
    pluto_constraints_negate_constraint(row, 0);
    pluto_constraints_add(check, row);
    if (pluto_constraints_is_empty(check)) {
      count++;
    }
    pluto_constraints_free(row);
  }
  IF_DEBUG(printf("%d constraints redundant\n", count););
}

int pluto_constraints_num_in_list(const PlutoConstraints *const cst) {
  if (cst == NULL)
    return 0;

  return 1 + pluto_constraints_num_in_list(cst->next);
}

/* In-place union: first argument is modified */
PlutoConstraints *pluto_constraints_unionize(PlutoConstraints *cst1,
                                             const PlutoConstraints *cst2) {
  PlutoConstraints *ucst = pluto_constraints_union(cst1, cst2);
  pluto_constraints_copy(cst1, ucst);
  pluto_constraints_free(ucst);
  return cst1;
}

/* In-place union: first argument is modified */
PlutoConstraints *pluto_constraints_unionize_isl(PlutoConstraints *cst1,
                                                 const PlutoConstraints *cst2) {
  PlutoConstraints *ucst = pluto_constraints_union_isl(cst1, cst2);
  pluto_constraints_copy(cst1, ucst);
  pluto_constraints_free(ucst);
  return cst1;
}

PlutoConstraints *pluto_constraints_universe(int ncols, PlutoContext *context) {
  PlutoConstraints *universe = pluto_constraints_alloc(1, ncols, context);
  return universe;
}

/* Return an empty polyhedron */
PlutoConstraints *pluto_constraints_empty(int ncols, PlutoContext *context) {
  PlutoConstraints *empty = pluto_constraints_alloc(1, ncols, context);
  pluto_constraints_add_equality(empty);
  empty->val[0][ncols - 1] = 1;
  return empty;
}

int pluto_constraints_is_empty(const PlutoConstraints *cst) {
  int64_t *sol;
  bool is_empty;

  if (cst->context->options->islsolve) {
    isl_ctx *ctx = isl_ctx_alloc();
    isl_set *iset = isl_set_from_pluto_constraints(cst, ctx);
    is_empty = isl_set_is_empty(iset);
    isl_set_free(iset);
    isl_ctx_free(ctx);
  } else {
    sol = pluto_constraints_lexmin_pip(cst, ALLOW_NEGATIVE_COEFF);
    is_empty = (sol == NULL);
    free(sol);
  }

  if (!is_empty)
    return 0;

  if (cst->next != NULL)
    return pluto_constraints_is_empty(cst->next);

  /* This one is empty and there are no more */
  return 1;
}

/* Append the constraints together */
PlutoConstraints *
pluto_constraints_unionize_simple(PlutoConstraints *cst1,
                                  const PlutoConstraints *cst2) {
  while (cst1->next != NULL) {
    cst1 = cst1->next;
  }
  cst1->next = pluto_constraints_dup(cst2);
  return cst1;
}

/* Get lower bound for pos^th variable if it's a (single) constant: return 1
 * if constant, 0 otherwise */
int pluto_constraints_get_const_lb(const PlutoConstraints *cnst, int pos,
                                   int64_t *lb) {
  int i, retval;

  PlutoConstraints *cst = pluto_constraints_dup(cnst);

  pluto_constraints_project_out_single(cst, 0, pos);
  pluto_constraints_project_out_single(cst, 1, cst->ncols - 2);

  retval = 0;
  *lb = INT_MIN;
  for (i = 0; i < cst->nrows; i++) {
    if (cst->is_eq[i]) {
      *lb = cst->val[i][cst->ncols - 1];
      retval = 1;
      break;
    }
    if (!cst->is_eq[i] && cst->val[i][0] >= 1) {
      retval = 1;
      *lb = PLMAX(*lb, -cst->val[i][cst->ncols - 1]);
    }
  }
  pluto_constraints_free(cst);

  if (retval && cnst->next != NULL) {
    int64_t next_lb;
    retval = pluto_constraints_get_const_lb(cnst->next, pos, &next_lb);
    if (*lb != next_lb)
      retval = 0;
  }

  return retval;
}

/* Get upper bound for the pos^th variable if it's a (single) constant: return
 * 0 if not constant,  1 otherwise */
int pluto_constraints_get_const_ub(const PlutoConstraints *cnst, int pos,
                                   int64_t *ub) {
  int i, retval;

  PlutoConstraints *cst = pluto_constraints_dup(cnst);

  pluto_constraints_project_out_single(cst, 0, pos);
  pluto_constraints_project_out_single(cst, 1, cst->ncols - 2);

  retval = 0;
  *ub = LLONG_MAX;
  for (i = 0; i < cst->nrows; i++) {
    if (cst->is_eq[i]) {
      *ub = cst->val[i][cst->ncols - 1];
      retval = 1;
      break;
    }
    if (!cst->is_eq[i] && cst->val[i][0] <= -1) {
      *ub = PLMIN(*ub, cst->val[i][cst->ncols - 1]);
      retval = 1;
    }
  }
  pluto_constraints_free(cst);

  if (retval && cnst->next != NULL) {
    int64_t next_ub;
    retval = pluto_constraints_get_const_ub(cnst->next, pos, &next_ub);
    if (*ub != next_ub)
      retval = 0;
  }

  return retval;
}

void print_polylib_visual_sets_internal(char *str, int k,
                                        PlutoConstraints *cst) {
  int i, j, first = 0;
  char name[100];

  if (!cst->context->options->moredebug)
    return;

  if (cst == NULL)
    return;

  sprintf(name, "%si%d", str, k);

  printf("%s := { ", name);
  for (i = 0; i < cst->ncols - 1; i++) {
    if (i == cst->ncols - 2)
      printf("t%d | ", i);
    else
      printf("t%d,", i);
  }

  for (i = 0; i < cst->nrows; i++) {
    first = 0;
    for (j = 0; j < cst->ncols - 1; j++) {
      if (cst->val[i][j] == 0)
        continue;

      if (first == 0) {
        first = 1;
        if (cst->val[i][j] == 1)
          printf("t%d", j);
        else if (cst->val[i][j] == -1)
          printf("-t%d", j);
        else if (cst->val[i][j] != 0)
          printf("%ld*t%d", cst->val[i][j], j);
      } else {
        first = 1;
        if (cst->val[i][j] == 1)
          printf("+t%d", j);
        else if (cst->val[i][j] == -1)
          printf("-t%d", j);
        else if (cst->val[i][j] != 0)
          printf("%+ld*t%d", cst->val[i][j], j);
      }
    }

    if (cst->val[i][j] != 0) {
      if (first != 0)
        printf("%+ld", cst->val[i][j]);
      else
        printf("%ld", cst->val[i][j]);
    }

    printf("%s0 ", cst->is_eq[i] ? "=" : ">=");
    if (i != cst->nrows - 1)
      printf(",");
  }

  printf("};\n\n");

  if (cst->next != NULL)
    print_polylib_visual_sets_internal(str, k + 1, cst->next);
}

void print_polylib_visual_sets_internal_new(char *str, int k,
                                            PlutoConstraints *cst) {

  int i, j, first = 0;
  char name[100];

  sprintf(name, "%s%d", str, k);

  printf("%s := { ", name);
  for (i = 0; i < cst->ncols - 1; i++) {
    if (i == cst->ncols - 2)
      printf("t%d | ", i);
    else
      printf("t%d,", i);
  }

  for (i = 0; i < cst->nrows; i++) {
    first = 0;
    for (j = 0; j < cst->ncols - 1; j++) {
      if (cst->val[i][j] == 0)
        continue;

      if (first == 0) {
        first = 1;
        if (cst->val[i][j] == 1)
          printf("t%d", j);
        else if (cst->val[i][j] == -1)
          printf("-t%d", j);
        else if (cst->val[i][j] != 0)
          printf("%ld*t%d", cst->val[i][j], j);
      } else {
        first = 1;
        if (cst->val[i][j] == 1)
          printf("+t%d", j);
        else if (cst->val[i][j] == -1)
          printf("-t%d", j);
        else if (cst->val[i][j] != 0)
          printf("%+ld*t%d", cst->val[i][j], j);
      }
    }

    if (cst->val[i][j] != 0) {
      if (first != 0)
        printf("%+ld", cst->val[i][j]);
      else
        printf("%ld", cst->val[i][j]);
    }

    printf("%s0 ", cst->is_eq[i] ? "=" : ">=");
    if (i != cst->nrows - 1)
      printf(",");
  }

  printf("};\n\n");

  if (cst->next != NULL)
    print_polylib_visual_sets_internal(str, k + 1, cst->next);
}

void print_polylib_visual_sets_new(char *name, PlutoConstraints *cst) {
  print_polylib_visual_sets_internal_new(name, 0, cst);
  return;
}

void print_polylib_visual_sets(char *name, PlutoConstraints *cst) {
  print_polylib_visual_sets_internal(name, 0, cst);
  return;
}

PlutoDepList *pluto_dep_list_alloc(Dep *dep) {
  PlutoDepList *list = (PlutoDepList *)malloc(sizeof(PlutoDepList));
  assert(dep != NULL);
  list->dep = dep;
  list->next = NULL;
  return list;
}

void pluto_deps_list_free(PlutoDepList *deplist) {
  if (deplist == NULL)
    return;
  pluto_deps_list_free(deplist->next);
}

PlutoDepList *pluto_deps_list_dup(PlutoDepList *src) {

  assert(src != NULL);

  PlutoDepList *newdeps = pluto_dep_list_alloc(src->dep);

  if (src->next != NULL) {
    newdeps->next = pluto_deps_list_dup(src->next);
  }

  return newdeps;
}

void pluto_deps_list_append(PlutoDepList *list, Dep *dep) {

  assert(list != NULL);
  assert(dep != NULL);

  PlutoDepList *newdeps = pluto_dep_list_alloc(dep);
  newdeps->next = list->next;
  list->next = newdeps;

  return;
}

PlutoConstraintsList *pluto_constraints_list_alloc(PlutoConstraints *cst) {

  PlutoConstraintsList *list =
      (PlutoConstraintsList *)malloc(sizeof(PlutoConstraintsList));
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
void pluto_constraints_list_add(PlutoConstraintsList *list,
                                const PlutoConstraints *cst, Dep *dep,
                                int copyDep) {

  assert(list != NULL);
  assert(cst != NULL);

  PlutoConstraintsList *newdeps =
      pluto_constraints_list_alloc(pluto_constraints_dup(cst));
  newdeps->next = list->next;
  list->next = newdeps;

  if (copyDep) {
    newdeps->deps = pluto_deps_list_dup(list->deps);
    pluto_deps_list_append(newdeps->deps, dep);
  } else {
    newdeps->deps = pluto_dep_list_alloc(dep);
  }
  return;
}

void pluto_constraints_list_replace(PlutoConstraintsList *list,
                                    PlutoConstraints *cst) {

  assert(list != NULL);
  assert(cst != NULL);

  pluto_constraints_free(list->constraints);
  list->constraints = pluto_constraints_dup(cst);

  return;
}

int pluto_constraints_get_num_equalities(const PlutoConstraints *cst) {
  int r, count;

  count = 0;
  for (r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r])
      count++;
  }
  return count;
}

PlutoMatrix *pluto_constraints_extract_equalities(const PlutoConstraints *cst) {
  int r, c, count;

  PlutoMatrix *mat = pluto_matrix_alloc(
      pluto_constraints_get_num_equalities(cst), cst->ncols, cst->context);

  count = 0;
  for (r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r]) {
      for (c = 0; c < cst->ncols; c++) {
        mat->val[count][c] = cst->val[r][c];
      }
      count++;
    }
  }

  return mat;
}

void pluto_constraints_remove_names_single(PlutoConstraints *cst) {
  int i;

  if (cst->names) {
    for (i = 0; i < cst->ncols - 1; i++) {
      free(cst->names[i]);
    }
  }
  free(cst->names);
  cst->names = NULL;
}

/*
 * Existing names if any are freed
 */
void pluto_constraints_set_names(PlutoConstraints *cst, char **names) {
  int i;

  if (cst->names) {
    for (i = 0; i < cst->ncols - 1; i++) {
      free(cst->names[i]);
    }
  }

  if (!cst->names) {
    cst->names = (char **)malloc((cst->ncols - 1) * sizeof(char *));
  }

  assert(names);

  for (i = 0; i < cst->ncols - 1; i++) {
    cst->names[i] = names[i] ? strdup(names[i]) : NULL;
  }
}

void pluto_constraints_set_names_range(PlutoConstraints *cst, char **names,
                                       int dest_offset, int src_offset,
                                       int num) {
  int i;

  if (num == 0)
    return;

  assert(names);
  assert(dest_offset + num <= cst->ncols - 1);

  if (!cst->names) {
    cst->names = (char **)malloc((cst->ncols - 1) * sizeof(char *));
    for (i = 0; i < cst->ncols - 1; i++) {
      cst->names[i] = NULL;
    }
  } else {
    for (i = 0; i < num; i++) {
      free(cst->names[dest_offset + i]);
    }
  }

  for (i = 0; i < num; i++) {
    if (names[src_offset + i] != NULL) {
      cst->names[dest_offset + i] = strdup(names[src_offset + i]);
    } else {
      cst->names[dest_offset + i] = NULL;
    }
  }
}

/*
 * Returns the best candidate to eliminate (exact index in cst)
 * max_elim: maximum number of variables to eliminate (from the right)
 *
 * FIXME: update it for constraints with equalities
 */
int pluto_constraints_best_elim_candidate(const PlutoConstraints *cst,
                                          int max_elim) {
  int64_t **csm;
  int64_t i, j, ub, lb, nb, num_eq, cost;

  int min_cost = cst->nrows * cst->nrows / 4;
  int best_candidate = cst->ncols - 2;

  csm = cst->val;

  for (j = cst->ncols - 2; j > cst->ncols - 2 - max_elim; j--) {
    ub = 0;
    lb = 0;
    nb = 0;
    num_eq = 0;
    for (i = 0; i < cst->nrows; i++) {
      if (cst->is_eq[i] && csm[i][j] != 0)
        num_eq++;
      else if (csm[i][j] > 0)
        ub++;
      else if (csm[i][j] < 0)
        lb++;
      else
        nb++;
    }
    if (num_eq >= 1)
      cost = cst->nrows - 1;
    else
      cost = lb * ub + nb;
    if (cost < min_cost) {
      min_cost = cost;
      best_candidate = j;
    }
  }

  return best_candidate;
}

PlutoConstraints *pluto_hyperplane_get_non_negative_half_space(Hyperplane *h) {
  assert(h->nrows == 1);
  assert(h->is_eq[0]);

  PlutoConstraints *pos_h = pluto_constraints_dup(h);
  pos_h->is_eq[0] = 0;
  return pos_h;
}

PlutoConstraints *pluto_hyperplane_get_negative_half_space(Hyperplane *h) {
  assert(h->nrows == 1);
  assert(h->is_eq[0]);

  PlutoConstraints *neg_h = pluto_constraints_dup(h);
  neg_h->is_eq[0] = 0;
  pluto_constraints_negate_row(neg_h, 0);
  neg_h->val[neg_h->nrows - 1][neg_h->ncols - 1] -= 1;
  return neg_h;
}

/* Shift a particular dimension by an affine function of other dimensions */
void pluto_constraints_shift_dim(PlutoConstraints *cst, int pos,
                                 PlutoMatrix *func) {
  int i, j;
  assert(func->ncols == cst->ncols);
  assert(func->nrows == 1);
  assert(func->val[0][pos] == 0);

  for (i = 0; i < cst->nrows; i++) {
    for (j = 0; j < cst->ncols; j++) {
      if (j != pos) {
        cst->val[i][j] -= cst->val[i][pos] * func->val[0][j];
      }
    }
  }
}

int pluto_constraints_is_ub(PlutoConstraints *cst, int row, int pos) {
  if (cst->val[row][pos] <= -1)
    return 1;

  return 0;
}

void pluto_constraints_remove_const_ub(PlutoConstraints *cst, int pos) {
  int i, j, r;
  assert(pos >= 0 && pos <= cst->ncols - 2);

  for (i = 0, r = 0; r < cst->nrows; r++) {
    if (pluto_constraints_is_ub(cst, i, pos)) {
      int sum = 0;
      for (j = 0; j < cst->ncols - 1 && j != pos; j++) {
        sum += abs(cst->val[i][j]);
      }
      if (sum == 0) {
        pluto_constraints_remove_row(cst, i);
      } else
        i++;
    }
  }
}

/* Returns the number of non-zero elements in the
 * constraint matrix, excluding the constant term */
int pluto_constraints_get_num_non_zero_coeffs(const PlutoConstraints *cst) {
  int nrows, ncols, non_zeros, i, j;

  nrows = cst->nrows;
  ncols = cst->ncols;
  non_zeros = 0;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols - 1; j++) {
      if (cst->val[i][j] != 0) {
        non_zeros++;
      }
    }
  }
  return non_zeros;
}

/*
 * Eliminates the pos^th variable, where pos has to be between 0 and
 *cst->ncols-2;
 * Remember that the last column is for the constant. The implementation does
 *not
 * have a complex redundancy check; it just uses pluto_constraints_simplify
 *which
 * eliminates duplicates after a gcd normalization and eliminates all zero
 * constraints
 *
 * Uses Gaussian elimination if there is an equality involving the variable
 */
void fourier_motzkin_eliminate_smart(PlutoConstraints *cst, unsigned pos) {
  int i, r, k, l, p, q;
  int64_t lb, ub, nb;
  int *bound;

  // At least one variable
  assert(cst->ncols >= 2);
  assert((int)pos <= (int)cst->ncols - 2);

  for (i = 0; i < cst->nrows; i++) {
    if (cst->is_eq[i] && cst->val[i][pos] != 0) {
      pluto_constraints_gaussian_eliminate(cst, pos);
      pluto_constraints_simplify(cst);
      return;
    }
  }

  PlutoConstraints *newcst;

  for (i = 0; i < cst->nrows; i++) {
    if (cst->val[i][pos] != 0)
      break;
  }

  if (i == cst->nrows) {
    newcst = pluto_constraints_dup_single(cst);
    pluto_constraints_remove_dim(newcst, pos);
  } else {
    bound = (int *)malloc(cst->nrows * sizeof(int));

    lb = 0;
    ub = 0;
    nb = 0;
    /* Variable does appear */
    for (r = 0; r < cst->nrows; r++) {
      if (cst->val[r][pos] == 0) {
        bound[r] = NB;
        nb++;
      } else if (cst->val[r][pos] >= 1) {
        bound[r] = LB;
        lb++;
      } else {
        bound[r] = UB;
        ub++;
      }
    }
    newcst =
        pluto_constraints_alloc(lb * ub + nb, cst->ncols - 1, cst->context);
    newcst->nrows = 0;

    p = 0;
    for (r = 0; r < cst->nrows; r++) {
      if (bound[r] == UB) {
        for (k = 0; k < cst->nrows; k++) {
          if (bound[k] == LB) {
            q = 0;
            for (l = 0; l < cst->ncols; l++) {
              if (l != pos) {
                newcst->val[p][q] =
                    cst->val[r][l] * (lcm(cst->val[k][pos], -cst->val[r][pos]) /
                                      (-cst->val[r][pos])) +
                    cst->val[k][l] * (lcm(-cst->val[r][pos], cst->val[k][pos]) /
                                      cst->val[k][pos]);
                q++;
              }
            }
            newcst->is_eq[p] = 0;
            p++;
          }
        }
      } else if (bound[r] == NB) {
        q = 0;
        for (l = 0; l < cst->ncols; l++) {
          if (l != pos) {
            newcst->val[p][q] = cst->val[r][l];
            q++;
          }
        }
        newcst->is_eq[p] = cst->is_eq[r];
        p++;
      }
    }
    assert(p == lb * ub + nb);
    newcst->nrows = p;
    free(bound);
  }

  pluto_constraints_simplify(newcst);
  pluto_constraints_copy_single(cst, newcst);
  pluto_constraints_free(newcst);

  if (cst->next != NULL)
    fourier_motzkin_eliminate(cst->next, pos);
}

void pluto_constraints_gaussian_eliminate(PlutoConstraints *cst, int pos) {
  int r, r2, c;
  int factor1, factor2;

  assert(pos >= 0);
  assert(pos <= cst->ncols - 2);

  for (r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r] && cst->val[r][pos] != 0) {
      break;
    }
  }

  if (r == cst->nrows) {
    printf("Can't eliminate dimension via GE\n");
    assert(0);
  }

  /* cst->val[r] is an equality */
  for (r2 = 0; r2 < cst->nrows; r2++) {
    if (r2 == r || cst->val[r2][pos] == 0)
      continue;
    if (cst->val[r2][pos] >= 1) {
      factor1 = lcm(llabs(cst->val[r][pos]), llabs(cst->val[r2][pos])) /
                cst->val[r2][pos];
      factor2 = lcm(llabs(cst->val[r][pos]), llabs(cst->val[r2][pos])) /
                cst->val[r][pos];
    } else if (cst->val[r2][pos] <= -1) {
      factor1 = -lcm(llabs(cst->val[r][pos]), llabs(cst->val[r2][pos])) /
                cst->val[r2][pos];
      factor2 = -lcm(llabs(cst->val[r][pos]), llabs(cst->val[r2][pos])) /
                cst->val[r][pos];
    }
    for (c = 0; c < cst->ncols; c++) {
      cst->val[r2][c] = cst->val[r2][c] * factor1 - cst->val[r][c] * factor2;
    }
  }
  pluto_constraints_remove_row(cst, r);

  pluto_constraints_remove_dim(cst, pos);

  if (cst->next != NULL)
    pluto_constraints_gaussian_eliminate(cst->next, pos);
}
