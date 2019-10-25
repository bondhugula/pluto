/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software: you can redistribute it and/or modify
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
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "constraints.h"
#include "ddg.h"
#include "isl_support.h"
#include "math_support.h"
#include "pluto.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "program.h"

#include "isl/constraint.h"
#include "isl/mat.h"
#include "isl/set.h"

#define CONSTRAINTS_SIMPLIFY_THRESHOLD 10000U
#define MAX_FARKAS_CST 2000

// This value represents the upper bound on w in case the dependence is short.
#define SHORT_DEP_DISTANCE 10

static int pluto_dep_satisfies_instance(const Dep *dep, const PlutoProg *prog,
                                        unsigned level);
static int pluto_dep_remove_satisfied_instances(Dep *dep, PlutoProg *prog,
                                                unsigned level);
static bool is_stmt_in_stmt_list(int stmt_id, std::vector<Stmt *> stmts);
static PlutoConstraints *
get_feautrier_schedule_constraints(PlutoProg *prog, std::vector<Stmt *> stmts);
/**
 *
 * Each constraint row has the following format
 *
 *      [dep distance bound | mapping coeff.s for S1, S2,... |constant]
 * Size:[       npar+1      | (nvar+1)*nstmts                | 1      ]
 *
 * npar - number of parameters in whole program
 * nvar - number of parameters in whole program
 */

/* Builds validity and bounding function constraints for a dependence */
static void compute_permutability_constraints_dep(Dep *dep, PlutoProg *prog) {
  PlutoConstraints *cst, *tiling_valid_cst, *bounding_func_cst;
  int j, r;
  PlutoMatrix *phi;
  PlutoContext *context = prog->context;
  PlutoOptions *options = prog->context->options;

  int nvar = prog->nvar;
  int npar = prog->npar;
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  /* IMPORTANT: It's assumed that all statements are of dimensionality nvar */

  IF_DEBUG(printf("[pluto] compute permutability constraints: Dep %d\n",
                  dep->id + 1););

  int dest_stmt = dep->dest;
  int src_stmt = dep->src;

  PlutoConstraints *dpoly = pluto_constraints_dup(dep->dpolytope);

  if (src_stmt != dest_stmt) {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, 2 * (nvar + 1) + 1, context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][r] = -1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][(nvar + 1) + (r - nvar)] = 1;
    }
    /* No parametric shifts: all zero for 2*nvar to 2*nvar+npar */

    /* Translation coefficients */
    phi->val[2 * nvar + npar][(nvar + 1) + nvar] = 1;
    phi->val[2 * nvar + npar][nvar] = -1;
  } else {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, (nvar + 1) + 1, context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][r] = -1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][r - nvar] = 1;
    }
    /* No parametric shifts: so all zero for 2*nvar to 2*nvar+npar-1 */

    /* Translation coefficients cancel out;
     * so nothing for 2*nvar+npar */
  }

  /* Apply Farkas lemma for tiling validity constraints */
  tiling_valid_cst = farkas_lemma_affine(dpoly, phi);

  pluto_matrix_free(phi);

  if (src_stmt != dest_stmt) {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, npar + 1 + 2 * (nvar + 1) + 1,
                             context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][npar + 1 + r] = 1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][npar + 1 + (nvar + 1) + (r - nvar)] = -1;
    }
    for (r = 2 * nvar; r < 2 * nvar + npar; r++) {
      /* for \vec{u} - parametric bounding function */
      phi->val[r][r - 2 * nvar] = 1;
    }

    /* Translation coefficients of statements */
    phi->val[2 * nvar + npar][npar + 1 + nvar] = 1;
    phi->val[2 * nvar + npar][npar + 1 + (nvar + 1) + nvar] = -1;
    /* for w */
    phi->val[2 * nvar + npar][npar] = 1;
  } else {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, npar + 1 + (nvar + 1) + 1,
                             context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][npar + 1 + r] = 1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][npar + 1 + (r - nvar)] = -1;
    }
    for (r = 2 * nvar; r < 2 * nvar + npar; r++) {
      /* for u */
      phi->val[r][r - 2 * nvar] = 1;
      /* No parametric shift coefficients */
    }
    /* Statement's translation coefficients cancel out */

    /* for w */
    phi->val[2 * nvar + npar][npar] = 1;
  }

  /* Apply Farkas lemma for bounding function constraints */
  bounding_func_cst = farkas_lemma_affine(dep->bounding_poly, phi);

  pluto_matrix_free(phi);
  pluto_constraints_free(dpoly);

  /* Aggregate permutability and bounding function constraints together in
   * global format; note that tiling_valid_cst and bounding_func_cst are
   * local to a  dependence/statements pertaining to it) */

  /* Everything initialized to zero during allocation */
  cst = pluto_constraints_alloc(
      tiling_valid_cst->nrows + bounding_func_cst->nrows, CST_WIDTH, context);
  cst->nrows = 0;
  cst->ncols = CST_WIDTH;

  int src_offset = npar + 1 + src_stmt * (nvar + 1);
  int dest_offset = npar + 1 + dest_stmt * (nvar + 1);

  /* Permutability constraints */
  if (!IS_RAR(dep->type)) {
    /* Permutability constraints only for non-RAR deps */
    for (unsigned k = 0; k < tiling_valid_cst->nrows; k++) {
      pluto_constraints_add_constraint(cst, tiling_valid_cst->is_eq[k]);
      for (j = 0; j < nvar + 1; j++) {
        cst->val[cst->nrows - 1][src_offset + j] = tiling_valid_cst->val[k][j];
        if (src_stmt != dest_stmt) {
          cst->val[cst->nrows - 1][dest_offset + j] =
              tiling_valid_cst->val[k][nvar + 1 + j];
        }
      }
      /* constant part */
      if (src_stmt == dest_stmt) {
        cst->val[cst->nrows - 1][cst->ncols - 1] =
            tiling_valid_cst->val[k][nvar + 1];
      } else {
        cst->val[cst->nrows - 1][cst->ncols - 1] =
            tiling_valid_cst->val[k][2 * nvar + 2];
      }
    }
  }

  /* Add bounding function constraints */
  if (!options->nodepbound) {
    /* Bounding function constraints in global format */
    PlutoConstraints *bcst_g;

    int src_offset = npar + 1 + src_stmt * (nvar + 1);
    int dest_offset = npar + 1 + dest_stmt * (nvar + 1);

    bcst_g =
        pluto_constraints_alloc(bounding_func_cst->nrows, CST_WIDTH, context);

    for (unsigned k = 0; k < bounding_func_cst->nrows; k++) {
      pluto_constraints_add_constraint(bcst_g, bounding_func_cst->is_eq[k]);
      for (j = 0; j < npar + 1; j++) {
        bcst_g->val[bcst_g->nrows - 1][j] = bounding_func_cst->val[k][j];
      }
      for (j = 0; j < nvar + 1; j++) {
        bcst_g->val[bcst_g->nrows - 1][src_offset + j] =
            bounding_func_cst->val[k][npar + 1 + j];
        if (src_stmt != dest_stmt) {
          bcst_g->val[bcst_g->nrows - 1][dest_offset + j] =
              bounding_func_cst->val[k][npar + 1 + nvar + 1 + j];
        }
      }
      /* constant part */
      if (src_stmt == dest_stmt) {
        bcst_g->val[bcst_g->nrows - 1][bcst_g->ncols - 1] =
            bounding_func_cst->val[k][npar + 1 + nvar + 1];
      } else {
        bcst_g->val[bcst_g->nrows - 1][bcst_g->ncols - 1] =
            bounding_func_cst->val[k][npar + 1 + 2 * nvar + 2];
      }
    }
    pluto_constraints_add(cst, bcst_g);

    pluto_constraints_free(dep->bounding_cst);
    dep->bounding_cst = bcst_g;
  }

  /* Coefficients of those dimensions that were added for padding
   * are of no utility */
  for (int k = 0; k < nvar; k++) {
    if (!stmts[src_stmt]->is_orig_loop[k]) {
      for (unsigned j = 0; j < cst->nrows; j++) {
        cst->val[j][src_offset + k] = 0;
      }
    }
    if (src_stmt != dest_offset && !stmts[dest_stmt]->is_orig_loop[k]) {
      for (unsigned j = 0; j < cst->nrows; j++) {
        cst->val[j][dest_offset + k] = 0;
      }
    }
  }

  pluto_constraints_free(dep->cst);
  dep->cst = cst;

  pluto_constraints_free(tiling_valid_cst);
  pluto_constraints_free(bounding_func_cst);
}

/* Computes permutability constraints for a dependence.
 * An interfacing routine to expose compute_permutability_constraints_dep
 */
void compute_pairwise_permutability(Dep *dep, PlutoProg *prog) {
  compute_permutability_constraints_dep(dep, prog);
}

/* Computes permutatbility constraints for all the dependences within the input
 * SCC */
PlutoConstraints *get_scc_permutability_constraints(int scc_id,
                                                    PlutoProg *prog) {
  int i, ndeps, src, dest;
  Dep *dep;
  PlutoContext *context = prog->context;

  ndeps = prog->ndeps;

  PlutoConstraints *scc_dep_cst = NULL;

  for (i = 0; i < ndeps; i++) {
    dep = prog->deps[i];
    src = dep->src;
    dest = dep->dest;
    if (dep_is_satisfied(dep)) {
      continue;
    }
    if (prog->stmts[src]->scc_id == scc_id &&
        prog->stmts[src]->scc_id == prog->stmts[dest]->scc_id) {
      if (dep->cst == NULL) {
        compute_permutability_constraints_dep(dep, prog);
      }
      if (scc_dep_cst == NULL) {
        scc_dep_cst = pluto_constraints_alloc(
            (prog->ddg->sccs[scc_id].size * dep->cst->nrows), dep->cst->ncols,
            context);
        scc_dep_cst->nrows = 0;
        scc_dep_cst->ncols = dep->cst->ncols;
      }
      pluto_constraints_add(scc_dep_cst, dep->cst);
    }
  }
  return scc_dep_cst;
}

/// Computes permutatbility constraints for all the dependences within the input
/// connected component.
PlutoConstraints *get_cc_permutability_constraints(int cc_id, PlutoProg *prog) {
  int ndeps = prog->ndeps;

  PlutoConstraints *cc_dep_cst = NULL;

  for (int i = 0; i < ndeps; i++) {
    Dep *dep = prog->deps[i];
    int src = dep->src;
    int dest = dep->dest;
    if (dep_is_satisfied(dep)) {
      continue;
    }
    if (prog->stmts[src]->cc_id == cc_id &&
        prog->stmts[src]->cc_id == prog->stmts[dest]->cc_id) {
      if (dep->cst == NULL) {
        compute_permutability_constraints_dep(dep, prog);
      }
      if (cc_dep_cst == NULL) {
        cc_dep_cst = pluto_constraints_alloc((ndeps * dep->cst->nrows),
                                             dep->cst->ncols, prog->context);
      }
      pluto_constraints_add(cc_dep_cst, dep->cst);
    }
  }
  return cc_dep_cst;
}

/* This function itself is NOT thread-safe for the same PlutoProg */
PlutoConstraints *get_permutability_constraints(PlutoProg *prog) {
  int nstmts, nvar, npar, ndeps, total_cst_rows;
  PlutoConstraints *globcst;
  Dep **deps;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  nstmts = prog->nstmts;
  ndeps = prog->ndeps;
  deps = prog->deps;
  nvar = prog->nvar;
  npar = prog->npar;

  FILE *skipfp = fopen("skipdeps.txt", "r");
  int *skipdeps = (int *)malloc(ndeps * sizeof(int));
  bzero(skipdeps, ndeps * sizeof(int));

  /* For debugging (skip deps listed here) */
  if (skipfp) {
    int num;
    while (!feof(skipfp)) {
      fscanf(skipfp, "%d", &num);
      skipdeps[num - 1] = 1;
      printf("\tskipping dep %d\n", num);
    }
  }

  total_cst_rows = 0;

  /* Compute the constraints and store them in dep->cst */
  for (int i = 0; i < ndeps; i++) {
    Dep *dep = deps[i];

    if (skipdeps[i])
      continue;

    if (options->rar == 0 && IS_RAR(dep->type)) {
      continue;
    }

    if (dep->cst == NULL) {
      /* First time, compute the constraints */
      compute_permutability_constraints_dep(dep, prog);

      IF_DEBUG(fprintf(stdout, "\tFor dep %d; num_constraints: %d\n", i + 1,
                       dep->cst->nrows));
    }
    total_cst_rows += dep->cst->nrows;
  }

  int ncols = CST_WIDTH;
  if (options->per_cc_obj) {
    ncols += (npar + 1) * prog->ddg->num_ccs;
  }
  if (!prog->globcst) {
    prog->globcst = pluto_constraints_alloc(total_cst_rows, ncols, context);
  }

  globcst = prog->globcst;

  globcst->ncols = ncols;
  globcst->nrows = 0;

  /* Add constraints to globcst */
  int i;
  unsigned inc;
  for (i = 0, inc = 0; i < ndeps; i++) {
    Dep *dep = deps[i];

    if (skipdeps[i])
      continue;

    /* print_polylib_visual_sets("BB_cst", dep->bounding_cst); */

    if (options->rar == 0 && IS_RAR(dep->type))
      continue;

    /* Note that dependences would be marked satisfied (in
     * pluto_auto_transform) only after all possible independent solutions
     * are found at a depth
     */
    if (dep_is_satisfied(dep)) {
      continue;
    }

    /* Subsequent calls can just use the old ones */
    if (options->per_cc_obj) {
      PlutoConstraints *newcst = pluto_constraints_dup(dep->cst);
      int num_cols_to_add = (npar + 1) * (prog->ddg->num_ccs);
      pluto_constraints_add_leading_dims(newcst, num_cols_to_add);
      int cc_id = prog->stmts[dep->src]->cc_id;
      unsigned source_col_offset = (npar + 1) * prog->ddg->num_ccs;
      unsigned target_col_offset = npar + 1 + (npar + 1) * cc_id;
      for (int j = 0; j < npar + 1; j++) {
        pluto_constraints_interchange_cols(newcst, source_col_offset + j,
                                           target_col_offset + j);
      }
      pluto_constraints_add(globcst, newcst);
      pluto_constraints_free(newcst);
    } else {
      pluto_constraints_add(globcst, dep->cst);
    }
    /* print_polylib_visual_sets("global", dep->cst); */

    IF_DEBUG(fprintf(stdout, "\tAfter dep %d; num_constraints: %d\n", i + 1,
                     globcst->nrows));
    /* This is for optimization as opposed to for correctness. We will
     * simplify constraints only if it crosses the threshold: at the time
     * it crosses the threshold or at 1000 increments thereon. Simplifying
     * each time slows down Pluto whenever there are few hundreds of
     * dependences. Not simplifying at all also leads to a slow down
     * because it leads to a large globcst and a number of constraits in
     * it are redundant */
    assert(globcst->nrows >= dep->cst->nrows);
    if (globcst->nrows >= CONSTRAINTS_SIMPLIFY_THRESHOLD + (3000 * inc) &&
        globcst->nrows - dep->cst->nrows <
            CONSTRAINTS_SIMPLIFY_THRESHOLD + (3000U * inc)) {
      pluto_constraints_simplify(globcst);
      inc++;
      IF_DEBUG(fprintf(stdout,
                       "\tAfter dep %d; num_constraints_simplified: %d\n",
                       i + 1, globcst->nrows));
    }
  }

  pluto_constraints_simplify(globcst);

  free(skipdeps);
  if (skipfp)
    fclose(skipfp);

  IF_DEBUG(fprintf(
      stdout,
      "\tAfter all dependences: num constraints: %d, num variables: %d\n",
      globcst->nrows, globcst->ncols - 1));

  return globcst;
}

/* Builds 1-dimensional schedule constraints (sequential schedule) for a
 * dependence */
PlutoConstraints *get_feautrier_schedule_constraints_dep(Dep *dep,
                                                         PlutoProg *prog) {
  PlutoConstraints *cst, *sched_valid_cst;
  int nstmts, nvar, npar, src_stmt, dest_stmt, j, k, r;
  int src_offset, dest_offset;
  PlutoMatrix *phi;
  Stmt **stmts;
  PlutoContext *context = prog->context;

  nvar = prog->nvar;
  npar = prog->npar;
  stmts = prog->stmts;
  nstmts = prog->nstmts;

  /* IMPORTANT: It's assumed that all statements are of dimensionality nvar */

  IF_DEBUG(
      printf("[pluto] get 1-d scheduling constraints: Dep %d\n", dep->id + 1););

  dest_stmt = dep->dest;
  src_stmt = dep->src;

  PlutoConstraints *dpoly = pluto_constraints_dup(dep->dpolytope);

  if (src_stmt != dest_stmt) {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, 2 * (nvar + 1) + 1, context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][r] = -1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][(nvar + 1) + (r - nvar)] = 1;
    }
    /* No parametric shifts: all zero for 2*nvar to 2*nvar+npar */

    /* Translation coefficients */
    phi->val[2 * nvar + npar][(nvar + 1) + nvar] = 1;
    phi->val[2 * nvar + npar][nvar] = -1;
    /* \phi(t) - \phi(s) - 1 >= 0: this is for the -1 */
    phi->val[2 * nvar + npar][2 * (nvar + 1)] = -1;
  } else {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, (nvar + 1) + 1, context);
    pluto_matrix_set(phi, 0);

    for (r = 0; r < nvar; r++) {
      /* Source stmt */
      phi->val[r][r] = -1;
    }
    for (r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][r - nvar] = 1;
    }
    /* No parametric shifts: so all zero for 2*nvar to 2*nvar+npar-1 */

    /* Translation coefficients cancel out;
     * so nothing for 2*nvar+npar */

    /* \phi(t) - \phi(s) - 1 >= 0: this is for the -1 */
    phi->val[2 * nvar + npar][nvar + 1] = -1;
  }

  /* Apply Farkas lemma */
  sched_valid_cst = farkas_lemma_affine(dpoly, phi);

  pluto_matrix_free(phi);
  pluto_constraints_free(dpoly);

  /* Put scheduling constraints in global format */

  /* Everything initialized to zero during allocation */
  cst = pluto_constraints_alloc(sched_valid_cst->nrows, CST_WIDTH, context);
  cst->nrows = 0;
  cst->ncols = CST_WIDTH;

  src_offset = npar + 1 + src_stmt * (nvar + 1);
  dest_offset = npar + 1 + dest_stmt * (nvar + 1);

  /* Permutability constraints */
  if (!IS_RAR(dep->type)) {
    /* Permutability constraints only for non-RAR deps */
    for (unsigned k = 0; k < sched_valid_cst->nrows; k++) {
      pluto_constraints_add_constraint(cst, sched_valid_cst->is_eq[k]);
      for (j = 0; j < nvar + 1; j++) {
        cst->val[cst->nrows - 1][src_offset + j] = sched_valid_cst->val[k][j];
        if (src_stmt != dest_stmt) {
          cst->val[cst->nrows - 1][dest_offset + j] =
              sched_valid_cst->val[k][nvar + 1 + j];
        }
      }
      /* constant part */
      if (src_stmt == dest_stmt) {
        cst->val[cst->nrows - 1][cst->ncols - 1] =
            sched_valid_cst->val[k][nvar + 1];
      } else {
        cst->val[cst->nrows - 1][cst->ncols - 1] =
            sched_valid_cst->val[k][2 * nvar + 2];
      }
    }
  }

  /* Coefficients of those dimensions that were added for padding
   * are of no utility */
  for (k = 0; k < nvar; k++) {
    if (!stmts[src_stmt]->is_orig_loop[k]) {
      for (unsigned j = 0; j < cst->nrows; j++) {
        cst->val[j][src_offset + k] = 0;
      }
    }
    if (src_stmt != dest_offset && !stmts[dest_stmt]->is_orig_loop[k]) {
      for (unsigned j = 0; j < cst->nrows; j++) {
        cst->val[j][dest_offset + k] = 0;
      }
    }
  }

  pluto_constraints_free(sched_valid_cst);

  return cst;
}

/*
 * 1-d affine schedule for a set of statements.
 */
PlutoConstraints *
get_feautrier_schedule_constraints(PlutoProg *prog, std::vector<Stmt *> stmts) {
  int ndeps = prog->ndeps;
  Dep **deps = prog->deps;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;

  PlutoConstraints *fcst = pluto_constraints_alloc(
      ndeps * nvar, (npar + 1 + prog->nstmts * (nvar + 1) + 1), context);

  /* Compute the constraints and store them */
  for (int i = 0, inc = 0; i < ndeps; i++) {
    Dep *dep = deps[i];

    if (context->options->rar == 0 && IS_RAR(dep->type)) {
      continue;
    }

    if (!is_stmt_in_stmt_list(dep->src, stmts) ||
        !is_stmt_in_stmt_list(dep->dest, stmts)) {
      continue;
    }

    PlutoConstraints *fcst_d =
        get_feautrier_schedule_constraints_dep(dep, prog);

    IF_DEBUG(fprintf(stdout, "\tFor dep %d; num_constraints: %d\n", i + 1,
                     fcst_d->nrows));
    IF_MORE_DEBUG(fprintf(stdout, "Constraints for dep %d\n", i + 1));
    IF_MORE_DEBUG(pluto_constraints_pretty_print(stdout, fcst_d));

    pluto_constraints_add(fcst, fcst_d);
    assert(fcst->nrows >= fcst_d->nrows);
    pluto_constraints_free(fcst_d);

    if (fcst->nrows >= CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000 * inc) &&
        fcst->nrows - fcst_d->nrows <
            CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000 * inc)) {
      pluto_constraints_simplify(fcst);
      IF_DEBUG(fprintf(stdout,
                       "\tAfter dep %d; num_constraints_simplified: %d\n",
                       i + 1, fcst->nrows));
      if (fcst->nrows >= CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000 * inc)) {
        inc++;
      }
    }
  }
  pluto_constraints_simplify(fcst);

  IF_DEBUG(fprintf(
      stdout,
      "\tAfter all dependences: num constraints: %d, num variables: %d\n",
      fcst->nrows, fcst->ncols - 1));
  IF_DEBUG2(pluto_constraints_pretty_print(stdout, fcst));

  return fcst;
}

/*
 * Returns linear independence constraints for a single statement.
 *
 * In particular, if H contains the first rows of an affine transformation,
 * then return a constraint on the coefficients of the next row that
 * ensures that this next row is linearly independent of the first rows.
 * Furthermore, the constraint is constructed in such a way that it allows
 * for a solution when combined with the other constraints on the coefficients
 * (currcst), provided any such constraint can be constructed.
 *
 * We do this by computing a basis for the null space of H and returning
 * a constraint that enforces the sum of these linear expressions
 * over the coefficients to be strictly greater than zero.
 * In this sum, some of the linear expressions may be negated to ensure
 * that a solution exists.
 *
 * The return value is a list of constraints, the first *orthonum corresponding
 * to the linear expressions that form a basis of the null space
 * and the final constraint the actual linear independence constraint.
 *
 * If the null space is 0-dimensional, *orthonum will be zero and the return
 * value is NULL
 */
PlutoConstraints **get_stmt_lin_ind_constraints(Stmt *stmt,
                                                const PlutoProg *prog,
                                                const PlutoConstraints *currcst,
                                                int *orthonum) {
  int nvar = prog->nvar;
  int npar = prog->npar;
  int nstmts = prog->nstmts;
  HyperplaneProperties *hProps = prog->hProps;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  IF_DEBUG(printf("[pluto] get_stmt_ortho constraints S%d\n", stmt->id + 1););

  /* Transformation has full column rank already */
  if (pluto_stmt_get_num_ind_hyps(stmt) >= stmt->dim_orig) {
    *orthonum = 0;
    return NULL;
  }

  /* Get rid of the variables that don't appear in the domain of this
   * statement and also beta rows */
  unsigned p = 0;
  for (int i = 0; i < nvar; i++) {
    if (stmt->is_orig_loop[i]) {
      p++;
    }
  }

  assert(stmt->trans != NULL);

  unsigned q = 0;
  for (unsigned j = 0; j < stmt->trans->nrows; j++) {
    if (hProps[j].type != H_SCALAR) {
      q++;
    }
  }

  isl_ctx *ctx = isl_ctx_alloc();
  assert(ctx);

  isl_mat *h = isl_mat_alloc(ctx, q, p);

  p = 0;
  for (int i = 0; i < nvar; i++) {
    if (stmt->is_orig_loop[i]) {
      unsigned q = 0;
      for (unsigned j = 0; j < stmt->trans->nrows; j++) {
        /* Skip rows of h that are zero */
        if (hProps[j].type != H_SCALAR) {
          h = isl_mat_set_element_si(h, q, p, stmt->trans->val[j][i]);
          q++;
        }
      }
      p++;
    }
  }

  h = isl_mat_right_kernel(h);

  PlutoMatrix *ortho = pluto_matrix_from_isl_mat(h, context);

  isl_mat_free(h);

  PlutoConstraints **orthcst =
      (PlutoConstraints **)malloc((nvar + 1) * sizeof(PlutoConstraints *));

  int ncols = CST_WIDTH;
  if (options->per_cc_obj) {
    ncols += (npar + 1) * prog->ddg->num_ccs;
  }

  for (int i = 0; i < nvar + 1; i++) {
    orthcst[i] = pluto_constraints_alloc(1, ncols, context);
    orthcst[i]->ncols = ncols;
  }

  /* All non-negative orthant only */
  /* An optimized version where the constraints are added as
   * c_1 >= 0, c_2 >= 0, ..., c_n >= 0, c_1+c_2+..+c_n >= 1,
   * basically, only looking in the orthogonal space where everything is
   * non-negative
   */

  /* Normalize ortho first */
  for (unsigned j = 0; j < ortho->ncols; j++) {
    if (ortho->val[0][j] == 0)
      continue;
    int colgcd = abs(ortho->val[0][j]);
    unsigned i;
    for (i = 1; i < ortho->nrows; i++) {
      if (ortho->val[i][j] == 0)
        break;
      colgcd = gcd(colgcd, abs(ortho->val[i][j]));
    }
    if (i == ortho->nrows) {
      if (colgcd > 1) {
        for (unsigned k = 0; k < ortho->nrows; k++) {
          ortho->val[k][j] /= colgcd;
        }
      }
    }
  }

  /* Fast linear independence check */
  isl_basic_set *isl_currcst;
  if (options->flic)
    isl_currcst = NULL;
  else
    isl_currcst = isl_basic_set_from_pluto_constraints(ctx, currcst);

  assert(p == ortho->nrows);
  p = 0;
  int coeff_offset = npar + 1 + ncols - CST_WIDTH;
  for (unsigned i = 0; i < ortho->ncols; i++) {
    isl_basic_set *orthcst_i;

    unsigned j = 0;
    for (int q = 0; q < nvar; q++) {
      if (stmt->is_orig_loop[q]) {
        orthcst[p]->val[0][coeff_offset + (stmt->id) * (nvar + 1) + q] =
            ortho->val[j][i];
        j++;
      }
    }
    orthcst[p]->nrows = 1;
    orthcst[p]->val[0][ncols - 1] = -1;
    if (!options->flic) {
      orthcst_i = isl_basic_set_from_pluto_constraints(ctx, orthcst[p]);
    }
    orthcst[p]->val[0][ncols - 1] = 0;

    if (!options->flic) {
      orthcst_i =
          isl_basic_set_intersect(orthcst_i, isl_basic_set_copy(isl_currcst));
      if (isl_basic_set_plain_is_empty(orthcst_i) ||
          isl_basic_set_is_empty(orthcst_i)) {
        pluto_constraints_negate_row(orthcst[p], 0);
      }
      isl_basic_set_free(orthcst_i);
    }
    p++;
  }

  if (p >= 1) {
    /* Sum of all of the above is the last constraint */
    for (int j = 0; j < ncols; j++) {
      for (unsigned i = 0; i < p; i++) {
        orthcst[p]->val[0][j] += orthcst[i]->val[0][j];
      }
    }
    orthcst[p]->nrows = 1;
    orthcst[p]->val[0][ncols - 1] = -1;
    p++;
  }

  *orthonum = p;

  IF_DEBUG2(printf("Ortho constraints for S%d; %d disjuncts\n", stmt->id + 1,
                   *orthonum - 1));
  for (int i = 0; i < *orthonum; i++) {
    // print_polylib_visual_sets("li", orthcst[i]);
    IF_DEBUG2(pluto_constraints_compact_print(stdout, orthcst[i]));
  }

  /* Free the unnecessary ones */
  for (int i = p; i < nvar + 1; i++) {
    pluto_constraints_free(orthcst[i]);
  }

  pluto_matrix_free(ortho);
  isl_basic_set_free(isl_currcst);
  isl_ctx_free(ctx);

  return orthcst;
}

/* PlutoConstraints to avoid trivial solutions (all zeros)
 *
 * hyp_search_mode = EAGER: If a statement's transformation is not full-ranked,
 * a hyperplane, if found, will be a loop hyperplane.
 *                 = LAZY: at least one of the hyperplanes for non-full
 *statements
 *  should be a loop hyperplane as opposed to all
 */
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *prog,
                                                  bool hyp_search_mode) {
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  int ncols = CST_WIDTH;
  if (options->per_cc_obj) {
    ncols += (npar + 1) * prog->ddg->num_ccs;
  }
  PlutoConstraints *nzcst = pluto_constraints_alloc(nstmts, ncols, context);
  nzcst->ncols = ncols;

  int coeff_offset = npar + 1 + ncols - CST_WIDTH;

  if (hyp_search_mode == EAGER) {
    for (int i = 0; i < nstmts; i++) {
      /* Don't add the constraint if enough solutions have been found */
      if (pluto_stmt_get_num_ind_hyps(stmts[i]) >= stmts[i]->dim_orig) {
        IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
        continue;
      }
      int stmt_offset = coeff_offset + i * (nvar + 1);
      for (int j = 0; j < nvar; j++) {
        if (stmts[i]->is_orig_loop[j] == 1) {
          nzcst->val[nzcst->nrows][stmt_offset + j] = 1;
        }
      }
      nzcst->val[nzcst->nrows][ncols - 1] = -1;
      nzcst->nrows++;
    }
  } else {
    assert(hyp_search_mode == LAZY);
    for (int i = 0; i < nstmts; i++) {
      /* Don't add the constraint if enough solutions have been found */
      if (pluto_stmt_get_num_ind_hyps(stmts[i]) >= stmts[i]->dim_orig) {
        IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
        continue;
      }
      int stmt_offset = coeff_offset + i * (nvar + 1);
      for (int j = 0; j < nvar; j++) {
        if (stmts[i]->is_orig_loop[j] == 1) {
          nzcst->val[0][stmt_offset + j] = 1;
        }
      }
      nzcst->val[0][ncols - 1] = -1;
    }
    nzcst->nrows = 1;
  }

  return nzcst;
}

/// Bounds for Pluto ILP variables.
PlutoConstraints *get_coeff_bounding_constraints(const PlutoProg *prog) {
  int npar = prog->npar;
  int nstmts = prog->nstmts;
  int nvar = prog->nvar;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  unsigned ncols = CST_WIDTH;
  if (options->per_cc_obj) {
    ncols += (npar + 1) * prog->ddg->num_ccs;
  }

  // Allocate more assuming both upper and lower bounds for each coefficient.
  PlutoConstraints *cst = pluto_constraints_alloc(2 * ncols, ncols, context);

  unsigned trans_coeff_offset = npar + 1 + ncols - CST_WIDTH;

  /* Lower bound for bounding coefficients (all non-negative) */
  for (unsigned i = 0; i < trans_coeff_offset; i++) {
    pluto_constraints_add_lb(cst, i, 0);
  }
  /* Lower bound for transformation coefficients (all non-negative) */
  for (int i = trans_coeff_offset; i < (int)cst->ncols - 1; i++) {
    IF_DEBUG2(
        printf("Adding lower bound %d for transformation coefficients\n", 0););
    pluto_constraints_add_lb(cst, i, 0);
  }

  if (options->coeff_bound != -1) {
    for (int i = trans_coeff_offset; i < (int)cst->ncols - 1; i++) {
      IF_DEBUG2(
          printf("Adding upper bound %d for transformation coefficients\n",
                 options->coeff_bound););
      pluto_constraints_add_ub(cst, i, options->coeff_bound);
    }
  } else {
    /* Add upper bounds for transformation coefficients */
    int ub = pluto_prog_get_largest_const_in_domains(prog);

    /* Putting too small an upper bound can prevent useful transformations;
     * also, note that an upper bound is added for all statements globally due
     * to the lack of an easy way to determine bounds for each coefficient to
     * prevent spurious transformations that involve shifts proportional to
     * loop bounds
     */
    if (ub >= 10) {
      for (int i = trans_coeff_offset; i < (int)cst->ncols - 1; i++) {
        IF_DEBUG2(
            printf("Adding upper bound %d for transformation coefficients\n",
                   ub););
        pluto_constraints_add_ub(cst, i, ub);
      }
    }
  }

  return cst;
}

/*
 * Check whether the dependence is satisfied at level 'level'
 * (works whether the dep is const or non-const, inter-stmt or
 * self edge
 */
bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, unsigned level) {
  PlutoConstraints *cst;
  int j, src_dim, dest_dim, npar;
  bool is_empty;
  PlutoContext *context = prog->context;

  npar = prog->npar;

  Stmt *src_stmt = prog->stmts[dep->src];
  Stmt *dest_stmt = prog->stmts[dep->dest];

  src_dim = src_stmt->dim;
  dest_dim = dest_stmt->dim;

  assert(src_stmt->trans != NULL);
  assert(dest_stmt->trans != NULL);
  assert(level < src_stmt->trans->nrows);
  assert(level < dest_stmt->trans->nrows);

  cst = pluto_constraints_alloc(2 * (1 + dep->dpolytope->nrows),
                                src_dim + dest_dim + npar + 1, context);

  /*
   * constraint format
   * \phi(src) - \phi (dest) >= 0
   * (reverse of satisfaction)
   */

  cst->is_eq[0] = 0;
  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = src_stmt->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = -dest_stmt->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar + 1; j++) {
    cst->val[0][j] = src_stmt->trans->val[level][j - dest_dim] -
                     dest_stmt->trans->val[level][j - src_dim];
  }

  cst->nrows = 1;

  pluto_constraints_add(cst, dep->dpolytope);

  /* If no solution exists, the dependence is satisfied, i.e., no points
   * satisfy \phi(src) - \phi(dest) <= 0 */
  is_empty = pluto_constraints_is_empty(cst);
  pluto_constraints_free(cst);

  return is_empty;
}

/*
 * Remove dependence instances from 'dep' that have been satisified at 'level'
 *
 * A dependence instance is satisfied if \phi(t) - \phi(s) >= 1; hence, those
 * <s,t> with \phi(t) - \phi(s) <= 0 are the unsatisfied ones. In fact, there
 * will be a violation if \phi(t) - \phi(s) <= -1
 *
 * Retval: true if at least one dependence instance was satisfied
 */
static int pluto_dep_remove_satisfied_instances(Dep *dep, PlutoProg *prog,
                                                unsigned level) {
  PlutoConstraints *cst;
  int j, src, dest, src_dim, dest_dim, retval;

  int npar = prog->npar;

  Stmt **stmts = prog->stmts;
  PlutoContext *context = prog->context;

  src = dep->src;
  dest = dep->dest;

  src_dim = prog->stmts[src]->dim;
  dest_dim = prog->stmts[dest]->dim;

  assert(level < stmts[src]->trans->nrows);
  assert(level < stmts[dest]->trans->nrows);

  cst = pluto_constraints_alloc(1 + dep->dpolytope->nrows,
                                src_dim + dest_dim + npar + 1, context);

  /*
   * constraint format
   * \phi(dest) - \phi (src) <= 0
   */

  cst->is_eq[0] = 0;
  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = stmts[src]->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = -stmts[dest]->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
    cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] -
                     stmts[dest]->trans->val[level][j - src_dim];
  }
  j = src_dim + dest_dim + npar;
  cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] -
                   stmts[dest]->trans->val[level][j - src_dim];

  cst->nrows = 1;

  pluto_constraints_intersect_isl(dep->depsat_poly, cst);

  retval = !pluto_constraints_is_empty(cst);

  pluto_constraints_free(cst);

  return retval;
}

/*
 * A precise dep satisfaction computation.
 * Returns: number of dependences satisfied
 */
int pluto_compute_dep_satisfaction_precise(PlutoProg *prog) {
  int i, num_satisfied;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] computing_dep_satisfaction_precise\n"););

  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i]->satisfied = false;
    prog->deps[i]->satisfaction_level = -1;
  }

  num_satisfied = 0;

  /* Piplib is not thread-safe (use multiple threads only with --islsolve) */
  /* #pragma omp parallel for if (options->islsolve) */
  for (i = 0; i < prog->ndeps; i++) {
    int level;
    Dep *dep = prog->deps[i];

    if (dep->depsat_poly != NULL) {
      pluto_constraints_free(dep->depsat_poly);
    }
    dep->depsat_poly = pluto_constraints_dup(dep->dpolytope);

    if (dep->satvec != NULL)
      free(dep->satvec);
    dep->satvec = (int *)malloc(prog->num_hyperplanes * sizeof(int));

    dep->satisfaction_level = -1;

    for (level = 0; level < prog->num_hyperplanes; level++) {
      dep->satvec[level] = pluto_dep_satisfies_instance(dep, prog, level);
      pluto_dep_remove_satisfied_instances(dep, prog, level);
      if (dep->satvec[level]) {
        dep->satisfaction_level = PLMAX(dep->satisfaction_level, level);
      }
      if (pluto_constraints_is_empty(dep->depsat_poly)) {
        dep->satisfied = true;
        IF_MORE_DEBUG(printf("\tdep %d satisfied\n", dep->id + 1););
        if (!IS_RAR(dep->type))
          num_satisfied++;
        break;
      }
    }
    if (level == prog->num_hyperplanes && !IS_RAR(dep->type)) {
      /* Dep has not been satisfied fully */
    }
    level++;
    for (; level < prog->num_hyperplanes; level++) {
      dep->satvec[level] = 0;
    }
  }
  IF_DEBUG(printf("\t %d (out of %d) dep(s) satisfied\n", num_satisfied,
                  prog->ndeps););
  return num_satisfied;
}

/* Retval: true if some iterations are satisfied */
static int pluto_dep_satisfies_instance(const Dep *dep, const PlutoProg *prog,
                                        unsigned level) {
  PlutoConstraints *cst;
  int j, src, dest, src_dim, dest_dim, retval;
  PlutoContext *context = prog->context;

  int npar = prog->npar;

  Stmt **stmts = prog->stmts;

  src = dep->src;
  dest = dep->dest;

  src_dim = prog->stmts[src]->dim;
  dest_dim = prog->stmts[dest]->dim;

  assert(level < stmts[src]->trans->nrows);
  assert(level < stmts[dest]->trans->nrows);

  cst = pluto_constraints_alloc(1 + dep->dpolytope->nrows,
                                src_dim + dest_dim + npar + 1, context);

  /* Constraint format: \phi(dest) - \phi (src) >= 1 */
  cst->is_eq[0] = 0;
  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = stmts[dest]->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] +
                     stmts[dest]->trans->val[level][j - src_dim];
  }
  j = src_dim + dest_dim + npar;
  cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] +
                   stmts[dest]->trans->val[level][j - src_dim] - 1;

  cst->nrows = 1;

  pluto_constraints_intersect_isl(cst, dep->depsat_poly);

  retval = !pluto_constraints_is_empty(cst);

  pluto_constraints_free(cst);

  return retval;
}

/* Direction vector component at level 'level'
 */
DepDir get_dep_direction(const Dep *dep, const PlutoProg *prog,
                         unsigned level) {
  PlutoConstraints *cst;
  int j, src, dest;
  PlutoContext *context = prog->context;

  int npar = prog->npar;
  Stmt **stmts = prog->stmts;

  src = dep->src;
  dest = dep->dest;

  Stmt *src_stmt = stmts[dep->src];
  Stmt *dest_stmt = stmts[dep->dest];

  int src_dim = src_stmt->dim;
  int dest_dim = dest_stmt->dim;

  assert(level < stmts[src]->trans->nrows);
  assert(level < stmts[dest]->trans->nrows);

  cst = pluto_constraints_alloc(2 * (2 + dep->dpolytope->nrows),
                                (src_dim + dest_dim) + npar + 1, context);

  /*
   * Check for zero
   *
   * To test \phi (dest) - \phi(src) = 0, we try
   *
   * \phi(dest) - \phi(src) >= 1
   */
  cst->is_eq[0] = 0;
  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = stmts[dest]->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] +
                     stmts[dest]->trans->val[level][j - src_dim];
  }
  cst->val[0][src_dim + dest_dim + npar] =
      -stmts[src]->trans->val[level][src_dim + npar] +
      stmts[dest]->trans->val[level][dest_dim + npar] - 1;
  cst->nrows = 1;

  pluto_constraints_add(cst, dep->dpolytope);

  bool is_empty = pluto_constraints_is_empty(cst);

  if (is_empty) {
    for (j = 0; j < src_dim; j++) {
      cst->val[0][j] = stmts[src]->trans->val[level][j];
    }
    for (j = src_dim; j < src_dim + dest_dim; j++) {
      cst->val[0][j] = -stmts[dest]->trans->val[level][j - src_dim];
    }
    for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
      cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] -
                       stmts[dest]->trans->val[level][j - src_dim];
    }
    cst->val[0][src_dim + dest_dim + npar] =
        stmts[src]->trans->val[level][src_dim + npar] -
        stmts[dest]->trans->val[level][dest_dim + npar] - 1;
    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    is_empty = pluto_constraints_is_empty(cst);

    /* If no solution exists, all points satisfy \phi (dest) - \phi (src) = 0 */
    if (is_empty) {
      pluto_constraints_free(cst);
      return DEP_ZERO;
    }
  }

  /*
   * Check for PLUS
   * Constraint format
   * \phi(dest) - \phi (src) <= -1
   * (reverse of plus)
   */

  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = stmts[src]->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = -stmts[dest]->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
    cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] -
                     stmts[dest]->trans->val[level][j - src_dim];
  }
  cst->val[0][src_dim + dest_dim + npar] =
      stmts[src]->trans->val[level][src_dim + npar] -
      stmts[dest]->trans->val[level][dest_dim + npar] - 1;

  cst->nrows = 1;

  pluto_constraints_add(cst, dep->dpolytope);

  is_empty = pluto_constraints_is_empty(cst);

  if (is_empty) {
    pluto_constraints_free(cst);
    return DEP_PLUS;
  }

  /*
   * Check for MINUS
   *
   * Constraint format
   * \phi(dest) - \phi (src) >= 1
   * reverse of minus, we alraedy know that it's not zero
   */

  for (j = 0; j < src_dim; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j];
  }
  for (j = src_dim; j < src_dim + dest_dim; j++) {
    cst->val[0][j] = stmts[dest]->trans->val[level][j - src_dim];
  }
  for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
    cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] +
                     stmts[dest]->trans->val[level][j - src_dim];
  }
  cst->val[0][src_dim + dest_dim + npar] =
      -stmts[src]->trans->val[level][src_dim + npar] +
      stmts[dest]->trans->val[level][dest_dim + npar] - 1;
  cst->nrows = 1;

  pluto_constraints_add(cst, dep->dpolytope);

  is_empty = pluto_constraints_is_empty(cst);
  pluto_constraints_free(cst);

  if (is_empty) {
    return DEP_MINUS;
  }

  /* Neither ZERO, nor PLUS, nor MINUS, has to be STAR */
  return DEP_STAR;
}

/* The routine is used to populate the csr matrices for scaling rational
 * solutions of pluto-lp. These matrices are used to construct constraints for
 * scaling MIP */
void populate_scaling_csr_matrices_for_pluto_program(int ***index,
                                                     double ***val, int nrows,
                                                     PlutoProg *prog) {
  int i, num_ccs, num_rows, stmt_offset, nstmts, cc_id;
  Stmt **stmts;

  nstmts = prog->nstmts;
  stmts = prog->stmts;

  ddg_compute_cc(prog);
  num_ccs = prog->ddg->num_ccs;

  *val = (double **)malloc(sizeof(double *) * nrows);
  *index = (int **)malloc(sizeof(int *) * nrows);
  for (i = 0; i < nrows; i++) {
    (*val)[i] = (double *)malloc(sizeof(double) * 3);
    (*index)[i] = (int *)malloc(sizeof(int) * 3);
  }

  num_rows = 0;
  stmt_offset = 0;
  for (i = 0; i < nstmts; i++) {
    cc_id = stmts[i]->cc_id;
    for (unsigned j = 0; j < stmts[i]->dim_orig + 1; j++) {
      (*index)[num_rows][0] = 0;
      (*val)[num_rows][0] = 0.0f;

      (*index)[num_rows][1] = cc_id + 1;
      /* val[num_row][i] will be set once the mip solution is found */

      (*index)[num_rows][2] = num_ccs + stmt_offset + j + 1;
      (*val)[num_rows][2] = -1.0;
      num_rows++;
    }
    stmt_offset += stmts[i]->dim_orig + 1;
  }

  /* Safety check. */
  assert(nrows == num_rows);
}

// Adds a hyperplane for each statement from the ILP solution.
void pluto_add_hyperplane_from_ilp_solution(int64_t *sol, PlutoProg *prog) {
  int nvar = prog->nvar;
  int npar = prog->npar;
  int nstmts = prog->nstmts;
  Stmt **stmts = prog->stmts;
  PlutoOptions *options = prog->context->options;

  pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_LOOP);

  int coeff_offset = npar + 1;
  if (options->per_cc_obj) {
    coeff_offset += (npar + 1) * prog->ddg->num_ccs;
  }

  for (int j = 0; j < nstmts; j++) {
    Stmt *stmt = stmts[j];
    pluto_stmt_add_hyperplane(stmt, H_UNKNOWN, stmt->trans->nrows);
    for (int k = 0; k < nvar; k++) {
      stmt->trans->val[stmt->trans->nrows - 1][k] =
          sol[coeff_offset + j * (nvar + 1) + k];
    }
    /* No parameteric shifts */
    for (int k = nvar; k < nvar + npar; k++) {
      stmt->trans->val[stmt->trans->nrows - 1][k] = 0;
    }
    /* Constant loop shift */
    stmt->trans->val[stmt->trans->nrows - 1][nvar + npar] =
        sol[coeff_offset + j * (nvar + 1) + nvar];

    stmt->hyp_types[stmt->trans->nrows - 1] =
        pluto_is_hyperplane_scalar(stmt, stmt->trans->nrows - 1) ? H_SCALAR
                                                                 : H_LOOP;
  }
}

/// The routine checks whether the statement given by stat_id is present in
/// the vector of input statements.
static bool is_stmt_in_stmt_list(int stmt_id, std::vector<Stmt *> stmts) {
  for (auto itr = stmts.begin(); itr != stmts.end(); itr++) {
    Stmt *stmt = *itr;
    if (stmt->id == stmt_id)
      return true;
  }
  return false;
}

/// Finds domain face that allows concurrent start (for diamond tiling)
/// FIXME: iteration space boundaries are assumed to be corresponding to
/// rectangular ones
/// Returns: matrix with row i being the concurrent start face for Stmt i
static PlutoMatrix *
get_face_with_concurrent_start_for_stmts(PlutoProg *prog,
                                         std::vector<Stmt *> stmts) {
  PlutoConstraints *bcst;
  PlutoMatrix *conc_start_faces;
  PlutoContext *context = prog->context;

  int npar = prog->npar;
  int nvar = prog->nvar;

  PlutoConstraints *fcst = get_feautrier_schedule_constraints(prog, stmts);

  bcst = get_coeff_bounding_constraints(prog);
  pluto_constraints_add(fcst, bcst);
  pluto_constraints_free(bcst);

  int64_t *sol = pluto_prog_constraints_lexmin(fcst, prog);
  pluto_constraints_free(fcst);

  if (!sol) {
    IF_DEBUG(printf("[pluto] get_face_with_concurrent_start: no valid 1-d "
                    "schedules \n"););
    return NULL;
  }

  conc_start_faces = pluto_matrix_alloc(stmts.size(), nvar + npar + 1, context);
  pluto_matrix_set(conc_start_faces, 0);

  for (int s = 0, _s = 0; s < prog->nstmts; s++) {
    if (!is_stmt_in_stmt_list(s, stmts))
      continue;
    assert(_s <= (int)stmts.size() - 1);
    for (int j = 0; j < nvar; j++) {
      conc_start_faces->val[_s][j] = sol[npar + 1 + s * (nvar + 1) + j];
    }
    conc_start_faces->val[_s][nvar + npar] =
        sol[npar + 1 + s * (nvar + 1) + nvar];
    _s++;
  }
  free(sol);

  IF_DEBUG(printf("[pluto] get_face_with_concurrent_start: 1-d schedules\n"););
  unsigned s = 0;
  for (auto stmt_itr = stmts.begin(); stmt_itr != stmts.end(); stmt_itr++) {
    Stmt *stmt = *stmt_itr;
    IF_DEBUG(printf("\tf(S%d) = ", stmt->id + 1););
    IF_DEBUG(pluto_affine_function_print(stdout, conc_start_faces->val[s++],
                                         nvar + npar,
                                         (const char **)stmt->domain->names););
    IF_DEBUG(printf("\n"););
  }

  /* 1-d schedule should be parallel to an iteration space boundary
   * FIXME: assuming canonical boundaries */
  for (s = 0; s < stmts.size(); s++) {
    unsigned nz = 0;
    for (int j = 0; j < nvar; j++) {
      if (conc_start_faces->val[s][j])
        nz++;
    }
    if (nz != 1)
      break;
  }

  if (s < stmts.size()) {
    pluto_matrix_free(conc_start_faces);
    IF_DEBUG(printf("[pluto] No iteration space faces with concurrent start "
                    "for all statements\n"););
    return NULL;
  }

  IF_DEBUG(
      printf(
          "[pluto] faces with concurrent start found for all statements\n"););

  return conc_start_faces;
}

/// Return the face with concurrent start for a band. This routine is the
/// interface between the routine get_face_with_concurrent_start_for_stmts and
/// pluto.c for diamond tiling.
PlutoMatrix *get_face_with_concurrent_start(PlutoProg *prog, Band *band) {
  PlutoContext *context = prog->context;
  IF_DEBUG(printf("[pluto] get_face_with_concurrent_start for band\n\t"););
  IF_DEBUG(pluto_band_print(band););

  std::vector<Stmt *> stmts;
  for (unsigned i = 0; i < band->loop->nstmts; i++) {
    stmts.push_back(band->loop->stmts[i]);
  }
  return get_face_with_concurrent_start_for_stmts(prog, stmts);
}

/// Returns the dependence distance bounding constraints for each dependence in
/// the scc.
static PlutoConstraints *
get_dep_distance_bounding_constraints_for_scc(int scc_id, PlutoProg *prog) {
  PlutoConstraints *dep_dist_bound_cst = NULL;
  PlutoContext *context = prog->context;
  Dep **deps = prog->deps;
  Stmt **stmts = prog->stmts;
  for (int i = 0; i < prog->ndeps; i++) {
    Dep *dep = deps[i];
    if (IS_RAR(dep->type))
      continue;
    int src_stmt = dep->src;
    int dest_stmt = dep->dest;
    if (stmts[src_stmt]->scc_id != scc_id || stmts[dest_stmt]->scc_id != scc_id)
      continue;

    if (!dep->bounding_cst)
      compute_permutability_constraints_dep(dep, prog);
    if (!dep_dist_bound_cst) {
      dep_dist_bound_cst = pluto_constraints_alloc(
          dep->bounding_cst->nrows, dep->bounding_cst->ncols, context);
      dep_dist_bound_cst->nrows = 0;
      dep_dist_bound_cst->ncols = dep->bounding_cst->ncols;
      pluto_constraints_add(dep_dist_bound_cst, dep->bounding_cst);
    } else
      dep_dist_bound_cst =
          pluto_constraints_add(dep_dist_bound_cst, dep->bounding_cst);
  }
  return dep_dist_bound_cst;
}

/// Return constraints of the form \phi(t) - \phi(s) >= -(u+w) for the input
/// dependence. This is used to lower bound the dependence distances to a small
/// value (Dep distances may even be negative). The callee has to set the value
/// of u to be zero if the dependence distance has to be short.
static PlutoConstraints *
compute_dep_dist_lb_constraints_for_stencil_check(Dep *dep, PlutoProg *prog) {
  int src_stmt = dep->src;
  int dest_stmt = dep->dest;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoMatrix *phi = NULL;
  PlutoContext *context = prog->context;

  if (src_stmt != dest_stmt) {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, npar + 1 + 2 * (nvar + 1) + 1,
                             context);
    pluto_matrix_set(phi, 0);

    for (int r = 0; r < nvar; r++) {
      /* Source stmt. */
      phi->val[r][npar + 1 + r] = -1;
    }
    for (int r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt. */
      phi->val[r][npar + 1 + (nvar + 1) + (r - nvar)] = 1;
    }
    for (int r = 2 * nvar; r < 2 * nvar + npar; r++) {
      /* for \vec{u} - parametric bounding function */
      phi->val[r][r - 2 * nvar] = 1;
    }

    /* Translation coefficients of statements */
    phi->val[2 * nvar + npar][npar + 1 + nvar] = -1;
    phi->val[2 * nvar + npar][npar + 1 + (nvar + 1) + nvar] = 1;
    /* for w. */
    phi->val[2 * nvar + npar][npar] = 1;
  } else {
    phi = pluto_matrix_alloc(2 * nvar + npar + 1, npar + 1 + (nvar + 1) + 1,
                             context);
    pluto_matrix_set(phi, 0);

    for (int r = 0; r < nvar; r++) {
      /* Source stmt. */
      phi->val[r][npar + 1 + r] = -1;
    }
    for (int r = nvar; r < 2 * nvar; r++) {
      /* Dest stmt */
      phi->val[r][npar + 1 + (r - nvar)] = 1;
    }
    for (int r = 2 * nvar; r < 2 * nvar + npar; r++) {
      /* for u */
      phi->val[r][r - 2 * nvar] = 1;
      /* No parametric shift coefficients */
    }
    /* Statement's translation coefficients cancel out */

    /* for w */
    phi->val[2 * nvar + npar][npar] = 1;
  }

  /* Apply Farkas lemma for bounding function constraints */
  PlutoConstraints *bounding_func_cst =
      farkas_lemma_affine(dep->bounding_poly, phi);

  pluto_matrix_free(phi);
  // Bring this to global format
  PlutoConstraints *bcst_g;

  int src_offset = npar + 1 + src_stmt * (nvar + 1);
  int dest_offset = npar + 1 + dest_stmt * (nvar + 1);
  int nstmts = prog->nstmts;

  bcst_g =
      pluto_constraints_alloc(bounding_func_cst->nrows, CST_WIDTH, context);

  for (unsigned k = 0; k < bounding_func_cst->nrows; k++) {
    pluto_constraints_add_constraint(bcst_g, bounding_func_cst->is_eq[k]);
    for (int j = 0; j < npar + 1; j++) {
      bcst_g->val[bcst_g->nrows - 1][j] = bounding_func_cst->val[k][j];
    }
    for (int j = 0; j < nvar + 1; j++) {
      bcst_g->val[bcst_g->nrows - 1][src_offset + j] =
          bounding_func_cst->val[k][npar + 1 + j];
      if (src_stmt != dest_stmt) {
        bcst_g->val[bcst_g->nrows - 1][dest_offset + j] =
            bounding_func_cst->val[k][npar + 1 + nvar + 1 + j];
      }
    }
    /* constant part */
    if (src_stmt == dest_stmt) {
      bcst_g->val[bcst_g->nrows - 1][bcst_g->ncols - 1] =
          bounding_func_cst->val[k][npar + 1 + nvar + 1];
    } else {
      bcst_g->val[bcst_g->nrows - 1][bcst_g->ncols - 1] =
          bounding_func_cst->val[k][npar + 1 + 2 * nvar + 2];
    }
  }
  pluto_constraints_free(bounding_func_cst);
  return bcst_g;
}

/// Get a lower bound for dependence distances this is used to check if the
/// dependences in the SCC are short. See routine
/// compute_dep_dist_lb_constraints_for_stencil_check for more details.
static PlutoConstraints *
get_dep_dist_lower_bound_constraints_for_scc(int scc_id, PlutoProg *prog) {
  PlutoContext *context = prog->context;
  Dep **deps = prog->deps;
  Stmt **stmts = prog->stmts;
  PlutoConstraints *dist_lb_csts = NULL;
  for (int i = 0; i < prog->ndeps; i++) {
    Dep *dep = deps[i];
    if (IS_RAR(dep->type))
      continue;
    int src_stmt = dep->src;
    int dest_stmt = dep->dest;
    if (stmts[src_stmt]->scc_id != scc_id || stmts[dest_stmt]->scc_id != scc_id)
      continue;

    PlutoConstraints *dep_dist_lb_cst =
        compute_dep_dist_lb_constraints_for_stencil_check(dep, prog);
    if (dist_lb_csts == NULL) {
      dist_lb_csts =
          pluto_constraints_alloc(prog->ndeps * dep_dist_lb_cst->nrows,
                                  dep_dist_lb_cst->ncols, context);
      dist_lb_csts->nrows = 1;
      dist_lb_csts->ncols = dep_dist_lb_cst->ncols;
    }
    dist_lb_csts = pluto_constraints_add(dist_lb_csts, dep_dist_lb_cst);
    pluto_constraints_free(dep_dist_lb_cst);
  }
  return dist_lb_csts;
}

/// Returns true if the scc in the ddg given by scc_id has the stencil pattern
/// of dependences. That is, the routine checks if there a face with concurrent
/// start and there are short dependences on either sides of the face with
/// concurrent start.
bool is_scc_stencil(int scc_id, PlutoProg *prog) {
  std::vector<Stmt *> stmts;
  PlutoContext *context = prog->context;
  for (int i = 0; i < prog->nstmts; i++) {
    if (prog->stmts[i]->scc_id == scc_id) {
      stmts.push_back(prog->stmts[i]);
    }
  }
  // The face with concurrent start is guaranteed to satisfy some dependence.
  PlutoMatrix *conc_start_faces =
      get_face_with_concurrent_start_for_stmts(prog, stmts);
  if (!conc_start_faces)
    return false;
  // printf("Concurrent start faces\n");
  // pluto_matrix_print(stdout, conc_start_faces);
  // Get dependence distance bounding constraints.
  PlutoConstraints *dep_bound_cst =
      get_dep_distance_bounding_constraints_for_scc(scc_id, prog);

  PlutoConstraints *dep_dist_lower_bound_cst =
      get_dep_dist_lower_bound_constraints_for_scc(scc_id, prog);

  // Bounding Constraints for coefficients.
  PlutoConstraints *coeff_bound_cst = pluto_constraints_alloc(
      dep_bound_cst->ncols, dep_bound_cst->ncols, context);
  coeff_bound_cst->nrows = 0;
  coeff_bound_cst->ncols = dep_bound_cst->ncols;
  unsigned npar = prog->npar;

  for (unsigned i = npar; i < dep_bound_cst->ncols - 1; i++) {
    pluto_constraints_add_lb(coeff_bound_cst, i, 0);
  }
  // Set \vec{u} to zero.
  for (unsigned i = 0; i < npar; i++) {
    pluto_constraints_add_equality(coeff_bound_cst);
    coeff_bound_cst->val[coeff_bound_cst->nrows - 1][i] = 1;
  }

  int short_dist_bound = SHORT_DEP_DISTANCE;
  // Upper bound on w for short dependence distance.
  pluto_constraints_add_ub(coeff_bound_cst, npar, short_dist_bound);
  // Lower bound on w for short dependence distance. Note that in case of
  // stencils this can be negative.
  // pluto_constraints_add_lb(coeff_bound_cst, npar, -short_dist_bound);

  PlutoConstraints *orthocst = pluto_constraints_alloc(
      conc_start_faces->nrows, dep_bound_cst->ncols, context);
  orthocst->nrows = conc_start_faces->nrows;
  orthocst->ncols = dep_bound_cst->ncols;
  unsigned count = 0;
  unsigned nvar = prog->nvar;
  for (auto stmt_itr = stmts.begin(); stmt_itr != stmts.end(); stmt_itr++) {
    Stmt *stmt = *stmt_itr;
    unsigned offset = npar + 1 + stmt->id * (nvar + 1);
    for (int i = 0; i < prog->nvar; i++) {
      if (conc_start_faces->val[count][i] == 0) {
        orthocst->val[count][offset + i] = 1;
      } else {
        // The coefficient corresponding to the dimension of the concurrent
        // start is zero.
        pluto_constraints_add_equality(orthocst);
        orthocst->val[orthocst->nrows - 1][offset + i] = 1;
      }
    }
    orthocst->val[count][orthocst->ncols - 1] = -1;
    count++;
  }

  dep_bound_cst = pluto_constraints_add(dep_bound_cst, orthocst);
  dep_bound_cst =
      pluto_constraints_add(dep_bound_cst, dep_dist_lower_bound_cst);
  dep_bound_cst = pluto_constraints_add(dep_bound_cst, coeff_bound_cst);

  bool is_stencil = true;
  int64_t *sol = pluto_constraints_lexmin(dep_bound_cst, ALLOW_NEGATIVE_COEFF);
  if (!sol) {
    is_stencil = false;
  } else if (sol[npar] == 0) {
    // w should be non zero.
    is_stencil = false;
  }
  free(sol);

  pluto_constraints_free(dep_dist_lower_bound_cst);
  pluto_constraints_free(dep_bound_cst);
  pluto_constraints_free(coeff_bound_cst);
  pluto_matrix_free(conc_start_faces);
  pluto_constraints_free(orthocst);
  return is_stencil;
}
