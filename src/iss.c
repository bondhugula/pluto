/* Index set splitting */

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constraints.h"
#include "math_support.h"
#include "pluto.h"
#include "program.h"

PlutoConstraints **get_lin_ind_constraints(PlutoMatrix *mat, int *orthonum) {
  int i, j, k;
  PlutoConstraints **orthcst;
  isl_ctx *ctx;
  isl_mat *h;

  assert(mat != NULL);

  int ndim = mat->ncols;

  ctx = isl_ctx_alloc();
  assert(ctx);

  // printf("Input matrix\n");
  // pluto_matrix_print(stdout, mat);

  h = isl_mat_alloc(ctx, mat->nrows, ndim);

  for (i = 0; i < mat->ncols; i++) {
    for (j = 0; j < mat->nrows; j++) {
      h = isl_mat_set_element_si(h, j, i, mat->val[j][i]);
    }
  }

  h = isl_mat_right_kernel(h);

  PlutoMatrix *ortho = pluto_matrix_from_isl_mat(h);

  isl_mat_free(h);

  orthcst =
      (PlutoConstraints **)malloc((ndim + 1) * sizeof(PlutoConstraints *));

  for (i = 0; i < ndim + 1; i++) {
    orthcst[i] = pluto_constraints_alloc(1, ndim + 1);
    orthcst[i]->ncols = ndim + 1;
  }

  /* All non-negative orthant only */
  /* An optimized version where the constraints are added as
   * c_1 >= 0, c_2 >= 0, ..., c_n >= 0, c_1+c_2+..+c_n >= 1
   *
   * basically only look in the orthogonal space where everything is
   * non-negative
   */

  /* Normalize ortho first */
  for (j = 0; j < ortho->ncols; j++) {
    if (ortho->val[0][j] == 0)
      continue;
    int colgcd = abs(ortho->val[0][j]);
    for (i = 1; i < ortho->nrows; i++) {
      if (ortho->val[i][j] == 0)
        break;
      colgcd = gcd(colgcd, abs(ortho->val[i][j]));
    }
    if (i == ortho->nrows) {
      if (colgcd > 1) {
        for (k = 0; k < ortho->nrows; k++) {
          ortho->val[k][j] /= colgcd;
        }
      }
    }
  }
  // printf("Ortho matrix\n");
  // pluto_matrix_print(stdout, ortho);

  for (i = 0; i < ortho->ncols; i++) {
    for (j = 0; j < ndim; j++) {
      orthcst[i]->val[0][j] = ortho->val[j][i];
    }
    orthcst[i]->nrows = 1;
    orthcst[i]->val[0][ndim] = -1;
    orthcst[i]->val[0][ndim] = 0;
  }

  // pluto_matrix_print(stdout, stmt->trans);

  if (ortho->ncols >= 1) {
    /* Sum of all of the above is the last constraint */
    for (j = 0; j < ndim + 1; j++) {
      for (i = 0; i < ortho->ncols; i++) {
        orthcst[ortho->ncols]->val[0][j] += orthcst[i]->val[0][j];
      }
    }
    orthcst[ortho->ncols]->nrows = 1;
    orthcst[ortho->ncols]->val[0][ndim] = -1;
    *orthonum = ortho->ncols + 1;
  } else
    *orthonum = 0;

  // printf("Ortho constraints: %d set(s)\n", *orthonum);
  // for (i=0; i<*orthonum; i++) {
  // print_polylib_visual_sets("li", orthcst[i]);
  // pluto_constraints_print(stdout, orthcst[i]);
  // }

  /* Free the unnecessary ones */
  for (i = *orthonum; i < ndim + 1; i++) {
    pluto_constraints_free(orthcst[i]);
  }

  pluto_matrix_free(ortho);
  isl_ctx_free(ctx);

  return orthcst;
}

/*
 * Index Set Splitting for close-to mid-point cutting
 *
 * Refer to PACT'14 paper on tiling periodic domains for
 * the formulation
 */
PlutoConstraints *pluto_find_iss(const PlutoConstraints **doms, int ndoms,
                                 int npar, PlutoConstraints *indcst) {
  int i, j, k, ndim;

  if (ndoms == 0)
    return NULL;

  assert(ndoms >= 0);

  const PlutoConstraints *dom0 = doms[0];

  for (i = 0; i < ndoms; i++) {
    assert((doms[i]->ncols - 1 - npar) % 2 == 0);
  }
  ndim = (doms[0]->ncols - 1 - npar) / 2;

  /*
   * [m | sigma(h) | h | P r| const ]
   *
   * The ISS is given by h.i = P.p + r, i \in dom
   * m = bound on distance from mid-point
   *
   * */
  int iss_cst_width = 2 + ndim + npar + 1 + 1;
  PlutoConstraints *cst = pluto_constraints_alloc(10, iss_cst_width);

  for (k = 0; k < ndoms; k++) {
    // printf("[iss] Domain %d\n", k);
    const PlutoConstraints *dom = doms[k];

    /* Linearize m + 2v(p) - h.s - h.t >= 0 */
    PlutoMatrix *mat = pluto_matrix_alloc(dom->ncols, iss_cst_width);
    pluto_matrix_set(mat, 0);

    for (i = 0; i < ndim; i++) {
      mat->val[i][2 + i] = -1;
    }
    for (i = 0; i < ndim; i++) {
      mat->val[ndim + i][2 + i] = -1;
    }
    for (i = 0; i < npar; i++) {
      mat->val[2 * ndim + i][2 + ndim + i] = 2;
    }
    /* m + 2*r for the const part */
    mat->val[2 * ndim + npar][2 + ndim + npar] = 2;
    mat->val[2 * ndim + npar][0] = 1;

    /*
     * cst is of the following form
     * [m | sigma(h) | h | P r| const ]
     */
    PlutoConstraints *cst1 = farkas_lemma_affine(dom, mat);

    /* Linearize: m - 2v(p) + h.s + h.t >= 0 */
    pluto_matrix_set(mat, 0);

    for (i = 0; i < ndim; i++) {
      mat->val[i][2 + i] = 1;
    }
    for (i = 0; i < ndim; i++) {
      mat->val[ndim + i][2 + i] = 1;
    }
    for (i = 0; i < npar; i++) {
      mat->val[2 * ndim + i][2 + ndim + i] = -2;
    }
    /* m + 2*r for the const part */
    mat->val[2 * ndim + npar][2 + ndim + npar] = -2;
    mat->val[2 * ndim + npar][0] = 1;

    PlutoConstraints *cst2 = farkas_lemma_affine(dom, mat);
    pluto_matrix_free(mat);

    pluto_constraints_add(cst, cst1);
    pluto_constraints_add(cst, cst2);
    pluto_constraints_free(cst1);
    pluto_constraints_free(cst2);
  }

  /* m <= 5 */
  pluto_constraints_add_ub(cst, 0, 5);

  /* Set sigma h */
  pluto_constraints_add_equality(cst);
  /* \sigma h_i */
  for (j = 2; j < 2 + ndim; j++) {
    cst->val[cst->nrows - 1][j] = -1;
  }
  cst->val[cst->nrows - 1][1] = 1;

  /* Avoid trivial zero solution of h */
  PlutoConstraints *nz = pluto_constraints_alloc(1, iss_cst_width);
  nz->nrows = 1;
  for (j = 0; j < ndim; j++) {
    nz->val[0][2 + j] = 1;
  }
  nz->val[0][nz->ncols - 1] = -1;
  pluto_constraints_add(cst, nz);

  /* Add linear independence constraints */
  if (indcst) {
    /* for m and sigma h */
    pluto_constraints_add_dim(indcst, 0, NULL);
    pluto_constraints_add_dim(indcst, 0, NULL);

    for (j = 0; j < npar + 1; j++) {
      pluto_constraints_add_dim(indcst, ndim + 2, NULL);
    }
    pluto_constraints_add(cst, indcst);
  }

  int64 *sol = pluto_constraints_lexmin(cst, DO_NOT_ALLOW_NEGATIVE_COEFF);

  pluto_constraints_free(cst);
  pluto_constraints_free(nz);

  if (sol) {
    PlutoConstraints *h = pluto_constraints_alloc(1, ndim + npar + 1);
    h->nrows = 1;

    h->is_eq[0] = 1;
    for (j = 0; j < ndim; j++) {
      h->val[0][j] = sol[2 + j];
    }
    for (j = 0; j < npar; j++) {
      h->val[0][ndim + j] = -sol[2 + ndim + j];
    }
    h->val[0][ndim + npar] = -sol[2 + ndim + npar];
    pluto_constraints_set_names_range(h, dom0->names, 0, 0, ndim);
    pluto_constraints_set_names_range(h, dom0->names, ndim, 2 * ndim, npar);

    PLUTO_MESSAGE(printf("[iss] m = %lld\n", sol[0]););
    PLUTO_MESSAGE(printf("[iss] h (cut) is "););
    PLUTO_MESSAGE(pluto_constraints_compact_print(stdout, h););
    free(sol);
    return h;
  } else {
    PLUTO_MESSAGE(printf("[iss] No solution to close to mid-point cut)\n"););
    return NULL;
  }
}

int is_long_bidirectional_dep(const Dep *dep, int dim, int npar) {
  assert(dep->src == dep->dest);

  // printf("Dep %d, dim; %d\n", dep->id+1, dim);
  // pluto_constraints_compact_print(stdout, dep->dpolytope);

  int ndim = (dep->dpolytope->ncols - 1 - npar) / 2;

  assert(dim >= 0);
  assert(dim <= ndim - 1);

  PlutoConstraints *dpolyc = pluto_constraints_dup(dep->dpolytope);
  /* Add new variable: t = i' - i  */
  pluto_constraints_add_dim(dpolyc, 0, NULL);
  pluto_constraints_add_equality(dpolyc);
  dpolyc->val[dpolyc->nrows - 1][0] = 1;
  dpolyc->val[dpolyc->nrows - 1][1 + dim] = 1;
  dpolyc->val[dpolyc->nrows - 1][1 + ndim + dim] = -1;

  int64 lb, ub;
  int retval1, retval2;

  retval1 = pluto_constraints_get_const_lb(dpolyc, 0, &lb);
  retval2 = pluto_constraints_get_const_ub(dpolyc, 0, &ub);

  pluto_constraints_free(dpolyc);

  // printf("lb = %lld, ub = %lld\n", lb, ub);

  return !(retval1 && retval2 && abs(ub) <= 5 && abs(lb) <= 5);
}

/*
 * Update dependences after ISS
 */
void pluto_update_deps_after_iss(PlutoProg *prog, PlutoConstraints **cuts,
                                 int num_cuts, int iss_stmt_id,
                                 int base_stmt_id) {
  int i, k, s, t, num_iss_deps;

  Dep **iss_deps = NULL;
  num_iss_deps = 0;

  if (num_cuts == 0)
    return;

  for (i = 0; i < prog->ndeps; i++) {
    int num_s_cuts, num_d_cuts;

    Dep *dep = prog->deps[i];

    if (dep->src != iss_stmt_id && dep->dest != iss_stmt_id) {
      num_iss_deps++;
      iss_deps = realloc(iss_deps, num_iss_deps * sizeof(Dep *));
      iss_deps[num_iss_deps - 1] = dep;
      continue;
    }
    if (dep->src == iss_stmt_id)
      num_s_cuts = num_cuts;
    else
      num_s_cuts = 1;
    if (dep->dest == iss_stmt_id)
      num_d_cuts = num_cuts;
    else
      num_d_cuts = 1;
    for (s = 0; s < num_s_cuts; s++) {
      for (t = 0; t < num_d_cuts; t++) {
        PlutoConstraints *dpolytope = pluto_constraints_dup(dep->dpolytope);

        Stmt *dest_stmt = prog->stmts[dep->dest];
        Stmt *src_stmt = prog->stmts[dep->src];
        if (dep->src == iss_stmt_id) {
          PlutoConstraints *scut = pluto_constraints_dup(cuts[s]);
          for (k = 0; k < dest_stmt->dim; k++) {
            pluto_constraints_add_dim(scut, src_stmt->dim, NULL);
          }
          pluto_constraints_add(dpolytope, scut);
          pluto_constraints_free(scut);
        }

        if (dep->dest == iss_stmt_id) {
          PlutoConstraints *dcut = pluto_constraints_dup(cuts[t]);
          for (k = 0; k < src_stmt->dim; k++) {
            pluto_constraints_add_dim(dcut, 0, NULL);
          }
          pluto_constraints_add(dpolytope, dcut);
          pluto_constraints_free(dcut);
        }

        if (!pluto_constraints_is_empty(dpolytope)) {
          num_iss_deps++;
          iss_deps = realloc(iss_deps, num_iss_deps * sizeof(Dep *));
          iss_deps[num_iss_deps - 1] = pluto_dep_dup(dep);

          Dep *iss_dep = iss_deps[num_iss_deps - 1];
          pluto_constraints_free(iss_dep->dpolytope);
          iss_dep->dpolytope = dpolytope;

          /* Update the source and target of the dependence */
          if (dep->src == iss_stmt_id) {
            iss_dep->src = base_stmt_id + s;
            iss_dep->src_acc = NULL;
          }
          if (dep->dest == iss_stmt_id) {
            iss_dep->dest = base_stmt_id + t;
            iss_dep->dest_acc = NULL;
          }
        } else {
          pluto_constraints_free(dpolytope);
        }
      }
    }
    pluto_dep_free(dep);
  }

  if (num_iss_deps >= 1) {
    free(prog->deps);

    prog->deps = iss_deps;
    prog->ndeps = num_iss_deps;
  }
}

/*
 * Perform Index Set Splitting: add new statements to the end
 */
void pluto_iss(Stmt *stmt, PlutoConstraints **cuts, int num_cuts,
               PlutoProg *prog) {
  int i;

  int prev_num_stmts = prog->nstmts;

  PLUTO_MESSAGE(printf("[iss] Splitting S%d into %d statements\n", stmt->id + 1,
                       num_cuts););

  for (i = 0; i < num_cuts; i++) {
    Stmt *nstmt = pluto_stmt_dup(stmt);
    pluto_constraints_add(nstmt->domain, cuts[i]);
    pluto_add_given_stmt(prog, nstmt);
  }

  pluto_update_deps_after_iss(prog, cuts, num_cuts, stmt->id, prev_num_stmts);
}

/*
 * Index set splitting based on near mid-point cutting of dependences
 */
void pluto_iss_dep(PlutoProg *prog) {
  int ndeps = prog->ndeps;
  int npar = prog->npar;
  int i, j;

  if (prog->nstmts == 0)
    return;

  /* TEMP restriction: all statements should have same dimensionality
   * and should be fused till the innermost level
   * TODO: change this to ISS for every statement with long self dependences
   * (potentially transitive)
   */
  for (i = 1; i < prog->nstmts; i++) {
    if (prog->stmts[i]->dim != prog->stmts[0]->dim) {
      return;
    }
  }
  if (!pluto_are_stmts_fused(prog->stmts, prog->nstmts, prog))
    return;

  if (ndeps == 0)
    return;

  int ndim = prog->stmts[0]->dim;

  int is_long[prog->ndeps][ndim];
  int num_long_deps[ndim];

  bzero(num_long_deps, ndim * sizeof(int));

  for (i = 0; i < ndeps; i++) {
    if (prog->deps[i]->src != prog->deps[i]->dest)
      continue;
    int ndim = (prog->deps[i]->dpolytope->ncols - 1 - npar) / 2;
    for (j = 0; j < ndim; j++) {
      is_long[i][j] = is_long_bidirectional_dep(prog->deps[i], j, npar);
      num_long_deps[j] += is_long[i][j];
    }
  }

  for (j = 0; j < ndim; j++) {
    PLUTO_MESSAGE(
        printf("[iss] Dimension %d has %d long deps\n", j, num_long_deps[j]););
  }

  PlutoConstraints ***long_dep_doms =
      (PlutoConstraints ***)malloc(sizeof(PlutoConstraints **) * ndim);
  for (i = 0; i < ndim; i++) {
    if (num_long_deps[i] >= 1) {
      long_dep_doms[i] = malloc(num_long_deps[i] * sizeof(PlutoConstraints *));
    } else
      long_dep_doms[i] = NULL;
  }

  for (j = 0; j < ndim; j++) {
    int q = 0;
    for (i = 0; i < ndeps; i++) {
      if (prog->deps[i]->src != prog->deps[i]->dest)
        continue;
      if (is_long[i][j]) {
        assert(q <= num_long_deps[j] - 1);
        long_dep_doms[j][q++] = prog->deps[i]->dpolytope;
      }
    }
  }

  int num_cuts;

  PlutoConstraints **cuts = NULL;
  int *pos = NULL;
  num_cuts = 0;

  for (i = 0; i < ndim; i++) {
    if (num_long_deps[i] == 0)
      continue;
    PLUTO_MESSAGE(printf("[iss] Dimension %d\n", i););
    PlutoConstraints *h =
        pluto_find_iss((const PlutoConstraints **)long_dep_doms[i],
                       num_long_deps[i], npar, NULL);
    if (h) {
      PlutoConstraints *negh = pluto_hyperplane_get_negative_half_space(h);
      PlutoConstraints *posh = pluto_hyperplane_get_non_negative_half_space(h);

      if (num_cuts == 0) {
        cuts = (PlutoConstraints **)malloc(2 * sizeof(PlutoConstraints *));
        pos = (int *)malloc(2 * sizeof(sizeof(int)));
        cuts[0] = negh;
        cuts[1] = posh;
        num_cuts = 2;
      } else if (num_cuts <= 2) {
        cuts = (PlutoConstraints **)realloc(
            cuts, 2 * num_cuts * sizeof(PlutoConstraints *));
        for (i = 0; i < num_cuts; i++) {
          cuts[num_cuts + i] = pluto_constraints_dup(cuts[i]);
          pluto_constraints_add(cuts[i], negh);
          pluto_constraints_add(cuts[num_cuts + i], posh);
        }
        num_cuts *= 2;
        pluto_constraints_free(negh);
        pluto_constraints_free(posh);
      }
      pluto_constraints_free(h);
    }
  }

  if (num_cuts >= 2) {
    int nstmts = prog->nstmts;
    /* Split the old statements */
    for (i = 0; i < nstmts; i++) {
      pluto_iss(prog->stmts[i], cuts, num_cuts, prog);
      // pluto_stmts_print(stdout, prog->stmts, prog->nstmts);
    }
    /* Remove old statements */
    for (i = 0; i < nstmts; i++) {
      pluto_remove_stmt(prog, 0);
    }
  }

  for (i = 0; i < num_cuts; i++) {
    pluto_constraints_free(cuts[i]);
  }
  free(cuts);
  free(pos);

  for (i = 0; i < ndim; i++) {
    free(long_dep_doms[i]);
  }
  free(long_dep_doms);

  IF_DEBUG(printf("After ISS\n"));
  IF_DEBUG(pluto_prog_print(stdout, prog));
}
