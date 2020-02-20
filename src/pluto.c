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
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>

#include "constraints.h"
#include "ddg.h"
#include "math_support.h"
#include "pluto.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"
#include "version.h"

void pluto_print_colours(int *colour, PlutoProg *prog);
bool *innermost_dep_satisfaction_dims(PlutoProg *prog,
                                      bool *tile_preventing_deps);
bool colour_scc(int scc_id, int *colour, int c, int stmt_pos, int pv,
                PlutoProg *prog);

int get_num_unsatisfied_deps(Dep **deps, int ndeps);
int get_num_unsatisfied_inter_stmt_deps(Dep **deps, int ndeps);
int get_num_unsatisfied_inter_scc_deps(PlutoProg *prog);

bool pluto_diamond_tile(PlutoProg *prog);

static double rtclock() {
  struct timeval Tp;
  int stat;
  stat = gettimeofday(&Tp, NULL);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

/*
 * Returns the number of (new) satisfied dependences at this level
 *
 * NOTE: for every unsatisfied dependence, this function tests if the entire
 * dependence has been satisfied at 'level'
 *
 */
int dep_satisfaction_update(PlutoProg *prog, int level) {
  int i;
  int num_new_carried;
  PlutoContext *context = prog->context;

  int ndeps = prog->ndeps;
  Dep **deps = prog->deps;

  num_new_carried = 0;

  IF_DEBUG(printf("[pluto] dep_satisfaction_update (level %d)\n", level););

  for (i = 0; i < ndeps; i++) {
    Dep *dep = deps[i];
    if (!dep_is_satisfied(dep)) {
      dep->satisfied = dep_satisfaction_test(dep, prog, level);
      if (dep->satisfied) {
        IF_MORE_DEBUG(
            printf("[pluto] dep_satisfaction_update: dep %d satisfied\n",
                   i + 1););
        if (!IS_RAR(dep->type))
          num_new_carried++;
        dep->satisfaction_level = level;
      }
    }
  }
  IF_DEBUG(printf("\t %d dep(s) satisfied\n", num_new_carried););

  return num_new_carried;
}

/* Check whether all deps are satisfied */
int deps_satisfaction_check(PlutoProg *prog) {
  int i;

  for (i = 0; i < prog->ndeps; i++) {
    if (IS_RAR(prog->deps[i]->type))
      continue;
    if (!dep_is_satisfied(prog->deps[i])) {
      return false;
    }
  }
  return true;
}

void pluto_dep_satisfaction_reset(PlutoProg *prog) {
  int i;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] pluto_dep_satisfaction_reset\n"););

  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i]->satisfied = false;
    prog->deps[i]->satisfaction_level = -1;
  }
}

/*
 * Conservative but powerful enough: until a dependence has been completely
 * satisfied (a level at which it is completely satisifed), a non-zero
 * dependence component would set satvec for that level to one
 */
void pluto_compute_dep_satisfaction(PlutoProg *prog) {
  int i, level;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] pluto_compute_dep_satisfaction\n"););

  pluto_dep_satisfaction_reset(prog);

  for (i = 0; i < prog->num_hyperplanes; i++) {
    dep_satisfaction_update(prog, i);
  }

  /* Create and set satisfaction vectors */
  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (IS_RAR(dep->type))
      continue;

    if (dep->satvec)
      free(dep->satvec);
    dep->satvec = (int *)malloc(prog->num_hyperplanes * sizeof(int));

    /* Direction vectors should be available */
    assert(dep->dirvec != NULL);

    for (level = 0; level < prog->num_hyperplanes; level++) {
      if (dep->dirvec[level] != DEP_ZERO &&
          (dep->satisfaction_level >= level || dep->satisfaction_level == -1)) {
        dep->satvec[level] = 1;
      } else {
        dep->satvec[level] = 0;
      }
    }
  }
}

bool dep_is_satisfied(Dep *dep) { return dep->satisfied; }

int num_satisfied_deps(Dep *deps, int ndeps) {
  int i;

  int num_satisfied = 0;
  for (i = 0; i < ndeps; i++) {
    if (IS_RAR(deps[i].type))
      continue;
    if (dep_is_satisfied(&deps[i]))
      num_satisfied++;
  }

  return num_satisfied;
}

int num_inter_stmt_deps(Dep *deps, int ndeps) {
  int i;
  int count;

  count = 0;
  for (i = 0; i < ndeps; i++) {
    if (IS_RAR(deps[i].type))
      continue;
    if (deps[i].src != deps[i].dest) {
      count++;
    }
  }
  return count;
}

int num_inter_scc_deps(Stmt *stmts, Dep *deps, int ndeps) {
  int i, count;

  count = 0;
  for (i = 0; i < ndeps; i++) {
    if (IS_RAR(deps[i].type))
      continue;
    if (dep_is_satisfied(&deps[i]))
      continue;
    if (stmts[deps[i].src].scc_id != stmts[deps[i].dest].scc_id)
      count++;
  }
  return count;
}

/**
 * Coefficient bounds when finding the cone complement; the cone complement
 * could have (and always has in the case of Pluto as opposed to Pluto+)
 * negative coefficients. So, we can't assume non-negative coefficients as in
 * the remaining Pluto hyperplanes
 */
static PlutoConstraints *
get_coeff_bounding_constraints_for_cone_complement(PlutoProg *prog) {
  int i, npar, nstmts, nvar;
  PlutoContext *context = prog->context;

  npar = prog->npar;
  nstmts = prog->nstmts;
  nvar = prog->nvar;

  PlutoConstraints *cst = pluto_constraints_alloc(1, CST_WIDTH, context);

  /* Lower bound for bounding coefficients */
  for (i = 0; i < npar + 1; i++) {
    pluto_constraints_add_lb(cst, i, 0);
  }
  /* Lower bound for transformation coefficients */
  for (int s = 0; s < nstmts; s++) {
    for (i = 0; i < nvar; i++) {
      /* Set this to -4 (is enough) */
      IF_DEBUG2(printf("[pluto_get_coeff_bound_for_cone_complement] Adding "
                       "lower bound %d for stmt dim coefficients\n",
                       -4););
      pluto_constraints_add_lb(cst, npar + 1 + s * (nvar + 1) + i, -4);
    }
    /* Translation coefficients need not be negative */
    IF_DEBUG2(printf("[pluto_get_coeff_bound_for_cone_complement] Adding lower "
                     "bound %d for stmt translation coefficient\n",
                     0););
    pluto_constraints_add_lb(cst, npar + 1 + s * (nvar + 1) + nvar, 0);
  }
  return cst;
}

PlutoMatrix *construct_cplex_objective(const PlutoConstraints *cst,
                                       const PlutoProg *prog) {
  int npar = prog->npar;
  int nvar = prog->nvar;
  PlutoMatrix *obj = pluto_matrix_alloc(1, cst->ncols - 1, prog->context);
  pluto_matrix_set(obj, 0);

  /* u */
  for (int j = 0; j < npar; j++) {
    obj->val[0][j] = 5 * 5 * nvar * prog->nstmts;
  }
  /* w */
  obj->val[0][npar] = 5 * nvar * prog->nstmts;

  for (int i = 0, j = npar + 1; i < prog->nstmts; i++) {
    unsigned k;
    for (k = j; k < j + prog->stmts[i]->dim_orig; k++) {
      obj->val[0][k] = (nvar + 2) * (prog->stmts[i]->dim_orig - (k - j));
    }
    /* constant shift */
    obj->val[0][k] = 1;
    j += prog->stmts[i]->dim_orig + 1;
  }
  return obj;
}

/*
 * This calls pluto_constraints_lexmin, but before doing that does some
 * preprocessing
 * - removes variables that we know will be assigned 0 - also do some
 *   permutation/substitution of variables
 */
int64_t *pluto_prog_constraints_lexmin(PlutoConstraints *cst, PlutoProg *prog) {
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;
  int nvar = prog->nvar;
  int npar = prog->npar;
  int num_ccs = prog->ddg->num_ccs;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  int ncols = CST_WIDTH;
  if (options->per_cc_obj) {
    ncols += (npar + 1) * num_ccs;
  }
  assert((int)cst->ncols - 1 == ncols - 1);

  unsigned coeff_offset = npar + 1 + ncols - CST_WIDTH;
  /* Remove redundant variables - that don't appear in your outer loops */
  int redun[coeff_offset + nstmts * (nvar + 1) + 1];
  for (unsigned i = 0; i < coeff_offset; i++) {
    redun[i] = 0;
  }

  for (int i = 0; i < nstmts; i++) {
    for (int j = 0; j < nvar; j++) {
      redun[coeff_offset + i * (nvar + 1) + j] = !stmts[i]->is_orig_loop[j];
    }
    redun[coeff_offset + i * (nvar + 1) + nvar] = 0;
  }
  redun[coeff_offset + nstmts * (nvar + 1)] = 0;

  /* TODO: Write a new in-place copy routine where some columns from the source
   * can be removed during copy using a mask. This will be helpful when a large
   * number of columns need to be removed. */
  unsigned newcols = 0;
  PlutoConstraints *newcst = pluto_constraints_dup(cst);
  for (int i = 0; i < cst->nrows; i++) {
    unsigned count = 0;
    for (int j = 0; j < cst->ncols; j++) {
      if (!redun[j]) {
        newcst->val[i][count++] = newcst->val[i][j];
      }
    }
    if (newcols == 0)
      newcols = count;
  }
  newcst->ncols = newcols;

  /* Permute the constraints so that if all else is the same, the original
   * hyperplane order is preserved (no strong reason to do this) */
  /* We do not need to permute in case of pluto-lp-dfp */
  if (!options->dfp) {
    unsigned j = coeff_offset;
    for (int i = 0; i < nstmts; i++) {
      for (unsigned k = j; k < j + (stmts[i]->dim_orig) / 2; k++) {
        pluto_constraints_interchange_cols(
            newcst, k, j + (stmts[i]->dim_orig - 1 - (k - j)));
      }
      j += stmts[i]->dim_orig + 1;
    }
  }
  IF_DEBUG(printf("[pluto] pluto_prog_constraints_lexmin (%d variables, %d "
                  "constraints)\n",
                  cst->ncols - 1, cst->nrows););

  int64_t *sol = NULL;
  /* Solve the constraints using the chosen solver. */
  if (options->islsolve) {
    double t_start = rtclock();
    sol = pluto_constraints_lexmin_isl(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
    prog->mipTime += rtclock() - t_start;
  } else if (options->glpk || options->lp || options->dfp || options->gurobi) {
    double **val = NULL;
    int **index = NULL;
    int nrows;

    nrows = 0;

    PlutoMatrix *obj = construct_cplex_objective(newcst, prog);

#if defined(GLPK) || defined(GUROBI)
    int num_ccs;

    num_ccs = 0;
#endif

    if (options->lp) {
      nrows = newcst->ncols - 1 - npar - 1;
      populate_scaling_csr_matrices_for_pluto_program(&index, &val, nrows,
                                                      prog);
#if defined(GLPK) || defined(GUROBI)
      num_ccs = prog->ddg->num_ccs;
#endif
    }

    double t_start = rtclock();
    if (options->glpk) {
#ifdef GLPK
      sol = pluto_prog_constraints_lexmin_glpk(newcst, obj, val, index, npar,
                                               num_ccs);
#endif
    } else if (options->gurobi) {
#ifdef GUROBI
      sol = pluto_prog_constraints_lexmin_gurobi(newcst, obj, val, index, npar,
                                                 num_ccs);
#endif
    }
    prog->mipTime += rtclock() - t_start;

    pluto_matrix_free(obj);
    if (options->lp) {
      for (int i = 0; i < nrows; i++) {
        free(val[i]);
        free(index[i]);
      }
      free(val);
      free(index);
    }
  } else {
    /* Use PIP */
    double t_start = rtclock();
    sol = pluto_constraints_lexmin_pip(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
    prog->mipTime += rtclock() - t_start;
  }

  int64_t *fsol = NULL;
  if (sol) {
    int k1, k2, q;
    int64_t tmp;
    /* Permute the solution in line with the permuted cst */
    if (!options->dfp) {
      unsigned j = coeff_offset;
      for (int i = 0; i < nstmts; i++) {
        for (unsigned k = j; k < j + (stmts[i]->dim_orig) / 2; k++) {
          k1 = k;
          k2 = j + (stmts[i]->dim_orig - 1 - (k - j));
          tmp = sol[k1];
          sol[k1] = sol[k2];
          sol[k2] = tmp;
        }
        j += stmts[i]->dim_orig + 1;
      }
    }

    fsol = (int64_t *)malloc((cst->ncols - 1) * sizeof(int64_t));

    /* Fill the soln with zeros for the redundant variables */
    q = 0;
    for (int j = 0; j < (int)cst->ncols - 1; j++) {
      fsol[j] = redun[j] ? 0 : sol[q++];
    }
    free(sol);
  }

  pluto_constraints_free(newcst);

  return fsol;
}

/* Is there an edge between some vertex of SCC1 and some vertex of SCC2? */
int ddg_sccs_direct_connected(Graph *g, PlutoProg *prog, int scc1, int scc2) {
  int i, j;

  for (i = 0; i < prog->nstmts; i++) {
    if (prog->stmts[i]->scc_id == scc1) {
      for (j = 0; j < prog->nstmts; j++) {
        if (prog->stmts[j]->scc_id == scc2) {
          if (g->adj->val[i][j] > 0) {
            return 1;
          }
        }
      }
    }
  }

  return 0;
}

/* Cut dependences between two SCCs
 * Returns: number of dependences cut  */
int cut_between_sccs(PlutoProg *prog, Graph *ddg, int scc1, int scc2) {
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;
  PlutoContext *context = prog->context;

  int nvar = prog->nvar;
  int npar = prog->npar;

  int i, j, num_satisfied;

  if (!ddg_sccs_direct_connected(ddg, prog, scc1, scc2)) {
    return 0;
  }

  IF_DEBUG(printf("[pluto] Cutting between SCC id %d and id %d\n", scc1, scc2));

  pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
    for (j = 0; j < nvar + npar; j++) {
      stmts[i]->trans->val[stmts[i]->trans->nrows - 1][j] = 0;
    }
    if (stmts[i]->scc_id < scc2) {
      stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + npar] = 0;
    } else {
      stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + npar] = 1;
    }
  }
  num_satisfied = dep_satisfaction_update(prog, stmts[0]->trans->nrows - 1);
  if (num_satisfied >= 1) {
    IF_DEBUG(
        pluto_transformation_print_level(prog, prog->num_hyperplanes - 1););
    ddg_update(ddg, prog);
  } else {
    for (i = 0; i < nstmts; i++) {
      stmts[i]->trans->nrows--;
    }
    prog->num_hyperplanes--;
  }

  return num_satisfied;
}

/*
 * Cut dependences between all SCCs
 */
int cut_all_sccs(PlutoProg *prog, Graph *ddg) {
  int i, j, num_satisfied;
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] Cutting between all SCCs\n"));

  if (ddg->num_sccs == 1) {
    IF_DEBUG(printf("\tonly one SCC\n"));
    return 0;
  }

  pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
    for (j = 0; j < nvar + npar; j++) {
      stmts[i]->trans->val[stmts[i]->trans->nrows - 1][j] = 0;
    }
    stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + npar] =
        stmts[i]->scc_id;
  }
  IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes - 1););
  num_satisfied = dep_satisfaction_update(prog, stmts[0]->trans->nrows - 1);
  ddg_update(ddg, prog);

  return num_satisfied;
}

/*
 * Cut based on dimensionalities of SCCs; if two SCCs are of different
 * dimensionalities; separate them
 * SCC1 -> SCC2 -> SCC3 ... ->SCCn
 * Two neighboring SCCs won't be cut if they are of the same
 * dimensionality
 *
 * (Due to the algorithm used for computing SCCs, there are no edges
 * between SCC<i> and SCC<j> where i > j)
 */
int cut_scc_dim_based(PlutoProg *prog, Graph *ddg) {
  int i, j, k, count;
  Stmt **stmts = prog->stmts;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;

  if (ddg->num_sccs == 1)
    return 0;

  IF_DEBUG(printf("Cutting based on SCC dimensionalities\n"));

  count = 0;

  int cur_max_dim = ddg->sccs[0].max_dim;

  pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

  for (k = 0; k < ddg->num_sccs; k++) {
    if (cur_max_dim != ddg->sccs[k].max_dim) {
      cur_max_dim = ddg->sccs[k].max_dim;
      count++;
    }

    for (i = 0; i < prog->nstmts; i++) {
      if (stmts[i]->scc_id == k) {
        pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
        for (j = 0; j < nvar; j++) {
          stmts[i]->trans->val[stmts[i]->trans->nrows - 1][j] = 0;
        }
        stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + npar] = count;
      }
    }
  }

  int num_new_carried =
      dep_satisfaction_update(prog, stmts[0]->trans->nrows - 1);

  if (num_new_carried >= 1) {
    IF_DEBUG(
        pluto_transformation_print_level(prog, prog->num_hyperplanes - 1););
    ddg_update(ddg, prog);
  } else {
    for (i = 0; i < prog->nstmts; i++) {
      stmts[i]->trans->nrows--;
    }
    prog->num_hyperplanes--;
  }

  return num_new_carried;
}

/* Heuristic cut */
void cut_smart(PlutoProg *prog, Graph *ddg) {
  if (ddg->num_sccs == 0)
    return;

  if (pluto_transformations_full_ranked(prog)) {
    /* Enough linearly independent solutions have been found */
    cut_all_sccs(prog, ddg);
    return;
  }

  int i, j;

  int num_new_carried = 0;

  /* First time, cut between SCCs of different dimensionalities */
  if (cut_scc_dim_based(prog, ddg)) {
    return;
  }

  /* Cut in the center */
  if (cut_between_sccs(prog, ddg, ceil(ddg->num_sccs / 2.0) - 1,
                       ceil(ddg->num_sccs / 2.0))) {
    return;
  }

  /* Cut between SCCs that are far away */
  for (i = 0; i < ddg->num_sccs - 1; i++) {
    for (j = ddg->num_sccs - 1; j >= i + 1; j--) {
      if ((int)prog->stmts[0]->trans->nrows <= 4 * prog->nvar + 2) {
        if (ddg_sccs_direct_connected(ddg, prog, i, j)) {
          // if (ddg->sccs[i].max_dim == ddg->sccs[j].max_dim) {
          num_new_carried += cut_between_sccs(prog, ddg, i, j);
          // }
        }
      } else {
        cut_all_sccs(prog, ddg);
        return;
      }
    }
  }
}

/* Distribute conservatively to maximize (rather random) fusion chance */
void cut_conservative(PlutoProg *prog, Graph *ddg) {
  int i, j;

  if (cut_scc_dim_based(prog, ddg)) {
    return;
  }

  /* Cut in the center */
  if (cut_between_sccs(prog, ddg, ceil(ddg->num_sccs / 2.0) - 1,
                       ceil(ddg->num_sccs / 2.0))) {
    return;
  }

  /* Cut between SCCs that are far away */
  for (i = 0; i < ddg->num_sccs - 1; i++) {
    for (j = ddg->num_sccs - 1; j >= i + 1; j--) {
      if ((int)prog->stmts[0]->trans->nrows <= 4 * prog->nvar + 2) {
        if (cut_between_sccs(prog, ddg, i, j)) {
          return;
        }
      } else {
        cut_all_sccs(prog, ddg);
        return;
      }
    }
  }
}

/*
 * Determine constraints to ensure linear independence of hyperplanes
 *
 * lin_ind_mode = EAGER: all statement hyperplanes have to be linearly
 *independent
 * w.r.t existing ones (ignoring stmts that already have enough lin ind solns)
 *              = LAZY: at least one statement that does not have enough
 * linearly independent solutions will get a new linearly independent
 * hyperplane (this is enough to make progress)
 */
PlutoConstraints *get_linear_ind_constraints(const PlutoProg *prog,
                                             const PlutoConstraints *cst,
                                             bool lin_ind_mode) {
  int npar, nvar, nstmts, i, j, k, orthosum;
  int orthonum[prog->nstmts];
  PlutoConstraints ***orthcst;
  Stmt **stmts;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] get_linear_ind_constraints\n"););

  npar = prog->npar;
  nvar = prog->nvar;
  nstmts = prog->nstmts;
  stmts = prog->stmts;

  orthcst = (PlutoConstraints ***)malloc(nstmts * sizeof(PlutoConstraints **));

  orthosum = 0;

  /* Get orthogonality constraints for each statement */
  for (j = 0; j < nstmts; j++) {
    orthcst[j] =
        get_stmt_lin_ind_constraints(stmts[j], prog, cst, &orthonum[j]);
    orthosum += orthonum[j];
  }

  int ncols = CST_WIDTH;
  if (prog->context->options->per_cc_obj) {
    ncols += (npar + 1) * prog->ddg->num_ccs;
  }
  PlutoConstraints *indcst = pluto_constraints_alloc(1, ncols, context);

  if (orthosum >= 1) {
    if (lin_ind_mode == EAGER) {
      /* Look for linearly independent hyperplanes for all stmts */
      for (j = 0; j < nstmts; j++) {
        if (orthonum[j] >= 1) {
          IF_DEBUG2(printf("Added ortho constraints for S%d\n", j + 1););
          pluto_constraints_add(indcst, orthcst[j][orthonum[j] - 1]);
        }
      }
    } else {
      assert(lin_ind_mode == LAZY);
      /* At least one stmt should have a linearly independent hyperplane */
      for (i = 0; i < prog->nstmts; i++) {
        /* Everything was initialized to zero */
        if (orthonum[i] >= 1) {
          for (j = 0; j < ncols - 1; j++) {
            indcst->val[0][j] += orthcst[i][orthonum[i] - 1]->val[0][j];
          }
        }
      }
      indcst->val[0][ncols - 1] = -1;
      indcst->nrows = 1;
      IF_DEBUG2(printf("Added \"at least one\" linear ind constraints\n"););
      IF_DEBUG2(pluto_constraints_pretty_print(stdout, indcst););
    }
  }

  for (j = 0; j < nstmts; j++) {
    for (k = 0; k < orthonum[j]; k++) {
      pluto_constraints_free(orthcst[j][k]);
    }
    free(orthcst[j]);
  }
  free(orthcst);

  return indcst;
}

// Each CC has its own u and w. Lexmin has to minimize u and w for all ccs
// equally. Therefore we add new set of u and w which is the component wise sum
// of u and w of the individual connected component. This is given by the
// constraints \wedge\limits_{i=0}^{npar+1} u_i = \sigma\limits_{j=0}^{num_ccs}
// u_i^j
PlutoConstraints *get_per_cc_obj_sum_constraints(PlutoProg *prog) {
  int num_ccs = prog->ddg->num_ccs;
  int npar = prog->npar;
  int nrows = npar + 1;
  int nvar = prog->nvar;
  int nstmts = prog->nstmts;
  int ncols = (num_ccs * (npar + 1)) + CST_WIDTH;
  PlutoConstraints *obj_sum_cst =
      pluto_constraints_alloc(nrows, ncols, prog->context);
  obj_sum_cst->nrows = nrows;
  obj_sum_cst->ncols = ncols;

  for (int i = 0; i < npar + 1; i++) {
    obj_sum_cst->val[i][i] = 1;
    for (int j = 1; j <= num_ccs; j++) {
      obj_sum_cst->val[i][j + i] = -1;
    }
  }
  return obj_sum_cst;
}

/* Find all linearly independent permutable band of hyperplanes at a level.
 *
 * See sub-functions for hyp_search_mode and lin_ind_mode
 *
 * If all statements already have enough linearly independent solutions, no
 * independence constraints will be generated, and since no non-trivial
 * solution constraints are added in such a case, the trivial zero solution will
 * end up being returned.
 * */
int find_permutable_hyperplanes(PlutoProg *prog, bool hyp_search_mode,
                                int max_sols, int band_depth) {
  int nstmts = prog->nstmts;
  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  IF_DEBUG(fprintf(stdout,
                   "[pluto] find_permutable_hyperplanes: max "
                   "solution(s): %d; band depth: %d\n",
                   max_sols, band_depth));

  assert(max_sols >= 0);

  if (max_sols == 0)
    return 0;

  /* Don't free basecst */
  PlutoConstraints *basecst = get_permutability_constraints(prog);
  PlutoConstraints *boundcst = get_coeff_bounding_constraints(prog);
  pluto_constraints_add(basecst, boundcst);
  pluto_constraints_free(boundcst);
  // print_polylib_visual_sets("pluto", basecst);

  int num_sols_found = 0;
  /* We don't expect to add a lot to basecst - just ortho constraints
   * and trivial soln avoidance constraints; instead of duplicating basecst,
   * we will just allocate once and copy each time */
  PlutoConstraints *currcst = pluto_constraints_alloc(
      basecst->nrows + nstmts + nvar * nstmts, CST_WIDTH, context);

  int64_t *bestsol;
  do {
    pluto_constraints_copy(currcst, basecst);
    PlutoConstraints *nzcst =
        get_non_trivial_sol_constraints(prog, hyp_search_mode);
    pluto_constraints_add(currcst, nzcst);
    pluto_constraints_free(nzcst);

    PlutoConstraints *indcst =
        get_linear_ind_constraints(prog, currcst, hyp_search_mode);
    // print_polylib_visual_sets("ind", indcst);
    IF_DEBUG2(printf("linear independence constraints\n"));
    IF_DEBUG2(pluto_constraints_pretty_print(stdout, indcst););

    if (indcst->nrows == 0) {
      /* If you don't have any independence constraints, we would end
       * up finding the same solution that was found earlier; so we
       * won't find anything new */
      IF_DEBUG(printf("No linearly independent rows\n"););
      bestsol = NULL;
    } else {
      pluto_constraints_add(currcst, indcst);
      if (options->per_cc_obj) {
        PlutoConstraints *obj_sum_cst = get_per_cc_obj_sum_constraints(prog);
        pluto_constraints_add(currcst, obj_sum_cst);
        pluto_constraints_free(obj_sum_cst);
      }
      IF_DEBUG(printf("[pluto] (Band %d) Solving for hyperplane #%d\n",
                      band_depth + 1, num_sols_found + 1));
      // IF_DEBUG2(pluto_constraints_pretty_print(stdout, currcst));
      bestsol = pluto_prog_constraints_lexmin(currcst, prog);
    }
    pluto_constraints_free(indcst);

    if (bestsol != NULL) {
      IF_DEBUG(fprintf(
          stdout, "[pluto] find_permutable_hyperplanes: found a hyperplane\n"));
      num_sols_found++;

      pluto_add_hyperplane_from_ilp_solution(bestsol, prog);

      free(bestsol);
      IF_DEBUG(
          pluto_transformation_print_level(prog, prog->num_hyperplanes - 1););
    } else {
      IF_DEBUG(fprintf(
          stdout,
          "[pluto] find_permutable_hyperplanes: No hyperplane found\n"));
    }
  } while (num_sols_found < max_sols && bestsol != NULL);

  pluto_constraints_free(currcst);

  /* Same number of solutions are found for each stmt */
  return num_sols_found;
}

/*
 * Returns H_LOOP if this hyperplane is a real loop or H_SCALAR if it's a scalar
 * dimension (beta row or node splitter)
 */
int get_loop_type(Stmt *stmt, int level) {
  for (int j = 0; j < (int)stmt->trans->ncols - 1; j++) {
    if (stmt->trans->val[level][j] > 0) {
      return H_LOOP;
    }
  }

  return H_SCALAR;
}

/* Cut based on the .fst file; returns 0 if it fails  */
bool precut(PlutoProg *prog, Graph *ddg, int depth) {
  int ncomps;

  int nstmts = prog->nstmts;

  int stmtGrp[nstmts][nstmts];
  int grpCount[nstmts];

  int i, j, k;

  if (depth != 0)
    return false;

  Stmt **stmts = prog->stmts;
  int nvar = prog->nvar;
  int npar = prog->npar;

  FILE *cutFp = fopen(".fst", "r");

  if (cutFp) {
    int tile;

    fscanf(cutFp, "%d", &ncomps);

    if (ncomps > nstmts) {
      printf("You have an .fst in your directory that is invalid for this "
             "source\n");
      printf("No fusion/distribution forced\n");
      return false;
    }

    for (i = 0; i < ncomps; i++) {
      fscanf(cutFp, "%d", &grpCount[i]);
      assert(grpCount[i] <= nstmts);
      for (j = 0; j < grpCount[i]; j++) {
        fscanf(cutFp, "%d", &stmtGrp[i][j]);
        assert(stmtGrp[i][j] <= nstmts - 1);
      }
      fscanf(cutFp, "%d", &tile);
      for (j = 0; j < grpCount[i]; j++)
        for (unsigned k = 0; k < stmts[stmtGrp[i][j]]->dim_orig; k++)
          stmts[stmtGrp[i][j]]->tile = tile;
    }

    fclose(cutFp);

    /* Update transformation matrices */
    for (i = 0; i < nstmts; i++) {
      pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
    }

    for (i = 0; i < ncomps; i++) {
      for (j = 0; j < grpCount[i]; j++) {
        int id = stmtGrp[i][j];
        for (k = 0; k < nvar + npar; k++) {
          stmts[id]->trans->val[stmts[id]->trans->nrows - 1][k] = 0;
        }
        stmts[id]->trans->val[stmts[id]->trans->nrows - 1][nvar + npar] = i;
      }
    }

    pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

    dep_satisfaction_update(prog, prog->num_hyperplanes - 1);
    ddg_update(ddg, prog);

    return true;
  } else {
    FILE *precut = fopen(".precut", "r");
    int ignore, rows, tile, tiling_depth;
    unsigned cols;

    if (precut) {
      /* Num of statements */
      fscanf(precut, "%d", &ignore);

      assert(ignore == prog->nstmts);

      /* Tiling depth */
      fscanf(precut, "%d", &tiling_depth);

      for (i = 0; i < prog->nstmts; i++) {
        /* Read scatterings */
        fscanf(precut, "%d", &rows);
        fscanf(precut, "%d", &cols);

        for (k = 0; k < rows; k++) {
          /* Transformation is in polylib format
           * <stmt_orig_dim>+<npar>+1 (first column for equality) */
          assert(cols == 1 + stmts[i]->dim_orig + npar + 1);

          /* For equality - ignore the zero */
          fscanf(precut, "%d", &ignore);
          assert(ignore == 0);

          pluto_stmt_add_hyperplane(stmts[i], H_UNKNOWN,
                                    stmts[i]->trans->nrows);

          for (j = 0; j < nvar; j++) {
            if (stmts[i]->is_orig_loop[j]) {
              fscanf(precut, "%ld",
                     &stmts[i]->trans->val[stmts[i]->trans->nrows - 1][j]);
            } else {
              stmts[i]->trans->val[stmts[i]->trans->nrows - 1][j] = 0;
            }
          }
          for (j = 0; j < npar; j++) {
            fscanf(precut, "%d", &ignore);
            stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + j] = 0;
          }
          /* Constant part */
          fscanf(
              precut, "%ld",
              &stmts[i]->trans->val[stmts[i]->trans->nrows - 1][nvar + npar]);
          if (get_loop_type(stmts[i], stmts[i]->trans->nrows - 1) == H_LOOP) {
            stmts[i]->hyp_types[stmts[i]->trans->nrows - 1] = H_LOOP;
          } else {
            stmts[i]->hyp_types[stmts[i]->trans->nrows - 1] = H_SCALAR;
          }
        }

        /* Number of tiling levels */
        fscanf(precut, "%d", &ignore);

        /* FIXME: to tile or not is specified depth-wise, why? Just
         * specify once */
        for (j = 0; j < tiling_depth; j++) {
          fscanf(precut, "%d", &tile);
        }
        stmts[i]->tile = tile;
      }

      /* Set hProps correctly and update satisfied dependences */
      for (k = 0; k < rows; k++) {
        pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_UNKNOWN);
        for (i = 0; i < nstmts; i++) {
          if (get_loop_type(stmts[i], stmts[0]->trans->nrows - rows + k) ==
              H_LOOP) {
            prog->hProps[prog->num_hyperplanes - 1].type = H_LOOP;
          } else {
            prog->hProps[prog->num_hyperplanes - 1].type = H_SCALAR;
          }
        }
        dep_satisfaction_update(prog, prog->num_hyperplanes - 1);
        ddg_update(ddg, prog);
      }
      return true;
    }

    return false;
  }
}

void pluto_compute_dep_directions(PlutoProg *prog) {
  int i, level;

  Dep **deps = prog->deps;

  for (i = 0; i < prog->ndeps; i++) {
    if (deps[i]->dirvec != NULL) {
      free(deps[i]->dirvec);
    }
    deps[i]->dirvec = (DepDir *)malloc(prog->num_hyperplanes * sizeof(DepDir));
    for (level = 0; level < prog->num_hyperplanes; level++) {
      deps[i]->dirvec[level] = get_dep_direction(deps[i], prog, level);
    }
  }
}

void pluto_detect_hyperplane_types_stmtwise(PlutoProg *prog) {
  for (int s = 0; s < prog->nstmts; s++) {
    Stmt *stmt = prog->stmts[s];
    for (unsigned i = 0; i < stmt->trans->nrows; i++) {
      stmt->hyp_types[i] =
          pluto_is_hyperplane_loop(stmt, i) ? H_LOOP : H_SCALAR;
    }
  }
}

/* Detect H_LOOP or H_SCALAR from scratch */
void pluto_detect_hyperplane_types(PlutoProg *prog) {
  int i, depth;
  int nstmts = prog->nstmts;

  for (depth = 0; depth < prog->num_hyperplanes; depth++) {
    for (i = 0; i < nstmts; i++) {
      if (pluto_is_hyperplane_loop(prog->stmts[i], depth))
        break;
    }
    prog->hProps[depth].type = (i < nstmts) ? H_LOOP : H_SCALAR;
  }
}

void pluto_print_depsat_vectors(PlutoProg *prog, int levels) {
  int i, j;
  Dep **deps;

  deps = prog->deps;

  printf("\nSatisfaction vectors for transformed program\n");

  for (i = 0; i < prog->ndeps; i++) {
    assert(deps[i]->satvec != NULL);
    printf("Dep %d: S%d to S%d: ", i + 1, deps[i]->src + 1, deps[i]->dest + 1);
    printf("(");
    for (j = 0; j < levels; j++) {
      printf("%d, ", deps[i]->satvec[j]);
    }
    printf(")\n");
  }
}

void pluto_print_dep_directions(PlutoProg *prog) {
  int i, j, ndeps, nlevels;
  Dep **deps;

  deps = prog->deps;
  ndeps = prog->ndeps;
  nlevels = prog->num_hyperplanes;

  printf("\nDirection vectors for transformed program\n");

  for (i = 0; i < ndeps; i++) {
    printf("Dep %d: S%d to S%d: ", i + 1, deps[i]->src + 1, deps[i]->dest + 1);
    printf("(");
    for (j = 0; j < nlevels; j++) {
      printf("%c, ", deps[i]->dirvec[j]);
    }
    assert(deps[i]->satvec != NULL);
    printf(") satisfied: %s, satvec: (", deps[i]->satisfied ? "yes" : " no");
    for (j = 0; j < nlevels; j++) {
      printf("%d, ", deps[i]->satvec[j]);
    }
    printf("), sat level: %d", deps[i]->satisfaction_level);

    for (j = 0; j < nlevels; j++) {
      if (deps[i]->dirvec[j] > 0) {
        break;
      }
      /* Just because this is not printed does not mean a dependence has
       * not been violated; it could still be */
      if (deps[i]->dirvec[j] < 0) {
        printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src + 1,
               deps[i]->dest + 1);
        printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
      } else if (deps[i]->dirvec[j] < 0) {
        printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src + 1,
               deps[i]->dest + 1);
        printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
      }
    }

    printf("\n");
  }
}

/*
 * 1. Pad statement domains to maximum domain depth to make it easier to
 * construct scheduling constraints. These will be removed before the actual ILP
 * runs, and just before pluto_auto_transform returns.
 *
 * 2. Pre-process dependence domains to allow a single affine bounding
 * expression to be constructed. Fix this implementation (too brittle).
 */
void normalize_domains(PlutoProg *prog) {
  int i, j, k;

  int nvar = prog->nvar;
  int npar = prog->npar;
  PlutoContext *context = prog->context;

  /* if a dep distance <= N and another <= M, how do you bound it, no
   * way to express max(N,M) as a single affine function, when space is
   * built for each dependence, it doesn't know anything about a parameter
   * that does not appear in its dpolyhedron, and so it will assign the
   * coeff corresponding to that param in the bounding function constraints
   * local to the dependence to zero, and what if some other dep needs that
   * particular coeff to be >= 1 for bounding?
   *
   * Solution: global context should be available
   *
   * How to construct? Just put together all constraints on parameters alone in
   * the global context, i.e., eliminate iterators out of each domain and
   * aggregate constraints on the parameters, and add them to each
   * dependence polyhedron
   */
  int count = 0;
  if (npar >= 1) {
    PlutoConstraints *param_context =
        pluto_constraints_alloc(prog->nstmts * npar, npar + 1, prog->context);
    pluto_constraints_set_names(param_context, prog->params);
    for (i = 0; i < prog->nstmts; i++) {
      PlutoConstraints *copy = pluto_constraints_dup(prog->stmts[i]->domain);
      for (unsigned j = 0; j < prog->stmts[i]->dim_orig; j++) {
        fourier_motzkin_eliminate(copy, 0);
      }
      assert((int)copy->ncols == npar + 1);
      count += copy->nrows;

      if (count <= prog->nstmts * npar) {
        pluto_constraints_add(param_context, copy);
        pluto_constraints_free(copy);
      } else {
        pluto_constraints_free(copy);
        break;
      }
    }
    pluto_constraints_simplify(param_context);
    IF_DEBUG(printf("[pluto] parameter context from domains\n"););
    IF_DEBUG(pluto_constraints_compact_print(stdout, param_context););

    /* Add context to every dep polyhedron */
    for (i = 0; i < prog->ndeps; i++) {
      PlutoConstraints *bounding_poly =
          pluto_constraints_dup(prog->deps[i]->dpolytope);

      for (unsigned k = 0; k < param_context->nrows; k++) {
        pluto_constraints_add_inequality(bounding_poly);

        /* Already initialized to zero */

        for (j = 0; j < npar + 1; j++) {
          bounding_poly->val[bounding_poly->nrows - 1]
                            [j + bounding_poly->ncols - (npar + 1)] =
              param_context->val[k][j];
        }
      }
      /* Update reference, add_row can resize */
      pluto_constraints_free(prog->deps[i]->bounding_poly);
      prog->deps[i]->bounding_poly = bounding_poly;
    }
    pluto_constraints_free(param_context);
  } else {
    IF_DEBUG(printf("\tNo global context\n"));
  }

  /* Add padding dimensions to statement domains */
  for (i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    for (j = stmt->dim; j < nvar; j++) {
      pluto_sink_statement(stmt, stmt->dim, 0, prog);
    }
  }

  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    int src_dim = prog->stmts[dep->src]->dim;
    int target_dim = prog->stmts[dep->dest]->dim;
    assert(dep->dpolytope->ncols ==
           (unsigned)src_dim + target_dim + prog->npar + 1);
  }

  /* Normalize rows of dependence polyhedra */
  for (k = 0; k < prog->ndeps; k++) {
    /* Normalize by gcd */
    PlutoConstraints *dpoly = prog->deps[k]->dpolytope;

    for (unsigned i = 0; i < dpoly->nrows; i++) {
      pluto_constraints_normalize_row(dpoly, i);
    }
  }

  /* Avoid the need for bounding function coefficients to take negative
   * values */
  bool *neg = (bool *)malloc(sizeof(bool) * npar);
  for (k = 0; k < prog->ndeps; k++) {
    Dep *dep = prog->deps[k];
    PlutoConstraints *dpoly = dep->dpolytope;

    bzero(neg, npar * sizeof(bool));

    if (dpoly->nrows == 0)
      continue;

    for (j = 2 * nvar; j < 2 * nvar + npar; j++) {
      int min = dpoly->val[0][j];
      int max = dpoly->val[0][j];
      for (unsigned i = 1; i < dpoly->nrows; i++) {
        min = PLMIN(dpoly->val[i][j], min);
        max = PLMAX(dpoly->val[i][j], max);
      }

      if (min < 0 && max <= 0) {
        neg[j - 2 * nvar] = true;
        IF_DEBUG(printf("Dep %d has negative coeff's for parameter %s\n",
                        dep->id, prog->params[j - 2 * nvar]));
      }
    }

    /* For parameters appearing with negative coefficients in upper bounds */
    for (j = 0; j < npar; j++) {
      if (neg[j]) {
        pluto_constraints_add_inequality(dep->bounding_poly);
        dep->bounding_poly->val[dep->bounding_poly->nrows - 1][2 * nvar + j] =
            1;
      }
    }
  }
  free(neg);
}

/* Remove padding dimensions that were added earlier; transformation matrices
 * will have stmt->dim + npar + 1 after this function */
void denormalize_domains(PlutoProg *prog) {
  int i, j;

  int nvar = prog->nvar;
  int npar = prog->npar;

  for (i = 0; i < prog->nstmts; i++) {
    int del_count;
    Stmt *stmt = prog->stmts[i];
    del_count = 0;
    for (j = 0; j < nvar; j++) {
      if (!stmt->is_orig_loop[j - del_count]) {
        pluto_stmt_remove_dim(stmt, j - del_count, prog);
        if (stmt->evicted_hyp) {
          pluto_matrix_remove_col(stmt->evicted_hyp, j - del_count);
        }
        del_count++;
      }
    }

    assert(stmt->domain->ncols == stmt->dim + npar + 1);
    assert(stmt->trans->ncols == stmt->dim + npar + 1);

    for (unsigned j = 0; j < stmt->dim; j++) {
      stmt->is_orig_loop[j] = 1;
    }
  }
}

/// Find hyperplane that completes the cone with previously found hyperplanes
/// such that the face allowing concurrent start lies within it.
/// conc_start_faces[i]: concurrent start face for statement $i$
/// evict_pos: position of the hyperplane to be evicted by the one that will
/// enable concurrent start cone_complement_pos: in case of partial concurrent
/// start, the hyperplane that will form the cone with the conc start hyperplane
/// cone_complement_hyps will set to the cone complement hyperplanes found
/// for statements in the band.
static int
find_cone_complement_hyperplane(Band *band, PlutoMatrix *conc_start_faces,
                                unsigned evict_pos, int cone_complement_pos,
                                PlutoConstraints *basecst, PlutoProg *prog,
                                PlutoMatrix **cone_complement_hyps) {
  int nvar = prog->nvar;
  int npar = prog->npar;
  unsigned nstmts = band->loop->nstmts;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  IF_DEBUG(printf("[pluto] find_cone_complement_hyperplane for band\n\t"););
  IF_DEBUG(pluto_band_print(band););

  // lambda_cst is the set of additional constraints added to find the cone
  // complement involving the conic combination multipliers.
  // TODO: improve comment by including info on the constraints that are added.
  PlutoConstraints *lambda_cst = pluto_constraints_alloc(
      2 * nvar * nstmts,
      (npar + 1 + prog->nstmts * (nvar + 1) + 1) + nvar * nstmts, context);

  /* All lambdas >=1 */
  for (unsigned i = 0; i < nstmts; i++) {
    int stmt_offset = npar + 1 + prog->nstmts * (nvar + 1) + i * nvar;
    for (int j = 0; j < nvar; j++) {
      pluto_constraints_add_inequality(lambda_cst);
      lambda_cst->val[lambda_cst->nrows - 1][stmt_offset + j] = 1;
      lambda_cst->val[lambda_cst->nrows - 1][lambda_cst->ncols - 1] = -1;
    }
  }

  /* Now, add the constraints for the new hyperplane to be in the cone
   * of the face and the negatives of the hyperplanes already found
   * (excluding the one being evicted: at `evict_pos'). */
  for (unsigned s = 0; s < nstmts; s++) {
    Stmt *stmt = band->loop->stmts[s];
    int trans_coeff_offset = npar + 1 + stmt->id * (nvar + 1);
    int lambda_offset = npar + 1 + prog->nstmts * (nvar + 1) + s * nvar;
    for (int j = 0; j < nvar; j++) {
      pluto_constraints_add_equality(lambda_cst);
      lambda_cst->val[lambda_cst->nrows - 1][trans_coeff_offset + j] = 1;

      lambda_cst->val[lambda_cst->nrows - 1][lambda_offset] =
          -(conc_start_faces->val[s][j]);

      /* Unless fulldiamondtile is set, enable concurrent start along
       * only one dimension. */
      if (!options->fulldiamondtile) {
        lambda_cst->val[lambda_cst->nrows - 1][lambda_offset + 1] =
            stmt->trans->val[cone_complement_pos][j];
      } else {
        // Full dimensional concurrent start.
        int lambda_k = 0;
        /* Just for the band depth hyperplanes */
        for (unsigned k = band->loop->depth;
             k < band->loop->depth + band->width; k++) {
          if (k != evict_pos && stmt->hyp_types[k] != H_SCALAR) {
            lambda_cst
                ->val[lambda_cst->nrows - 1][lambda_offset + lambda_k + 1] =
                stmt->trans->val[k][j];
            lambda_k++;
          }
        }
      }
      lambda_cst->val[lambda_cst->nrows - 1][lambda_cst->ncols - 1] = 0;
    }
  }

  /*
   * con_start_cst serves the same purpose as Pluto's ILP formulation, but
   * with expanded constraint-width to incorporate lambdas.
   * No need of non-zero solution constraints here.
   */
  PlutoConstraints *con_start_cst = pluto_constraints_dup(basecst);
  PlutoConstraints *boundcst =
      get_coeff_bounding_constraints_for_cone_complement(prog);
  pluto_constraints_add(con_start_cst, boundcst);
  pluto_constraints_free(boundcst);

  for (int i = 0; i < nvar * nstmts; i++) {
    pluto_constraints_add_dim(con_start_cst, basecst->ncols - 1, NULL);
  }

  pluto_constraints_add(con_start_cst, lambda_cst);
  pluto_constraints_free(lambda_cst);

  IF_MORE_DEBUG(printf("Cone complement constraints\n"););
  IF_MORE_DEBUG(pluto_constraints_pretty_print(stdout, con_start_cst););

  int64_t *bestsol =
      pluto_constraints_lexmin(con_start_cst, ALLOW_NEGATIVE_COEFF);
  pluto_constraints_free(con_start_cst);

  if (bestsol == NULL) {
    printf("[pluto] Cone complement hyperplane not found!\n");
    printf("[pluto] No tiled concurrent start possible.\n");
  } else {
    IF_DEBUG(printf("[pluto] Tiled concurrent start possible.\n"););
    for (unsigned j = 0; j < nstmts; j++) {
      Stmt *stmt = band->loop->stmts[j];
      cone_complement_hyps[j] =
          pluto_matrix_alloc(1, stmt->dim + npar + 1, context);
      for (int k = 0; k < nvar; k++) {
        cone_complement_hyps[j]->val[0][k] =
            bestsol[npar + 1 + stmt->id * (nvar + 1) + k];
      }
      /* No parametric shifts. */
      for (int k = nvar; k < nvar + npar; k++) {
        cone_complement_hyps[j]->val[0][k] = 0;
      }
      cone_complement_hyps[j]->val[0][nvar + npar] =
          bestsol[npar + 1 + stmt->id * (nvar + 1) + nvar];

      IF_DEBUG(printf("\tcone_complement(S%d) = ", stmt->id + 1););
      IF_DEBUG(
          pluto_affine_function_print(stdout, cone_complement_hyps[j]->val[0],
                                      nvar, (const char **)stmt->iterators););
      IF_DEBUG(printf("\n"););
    }
    free(bestsol);
  }

  return (cone_complement_hyps[0] == NULL) ? 0 : 1;
}

/*
 * Check if the d'th row for any statement in the band is parallel to the
 * face allowing concurrent start
 */
int is_concurrent_start_face(Band *band, PlutoMatrix *conc_start_faces, int d) {
  for (unsigned s = 0; s < band->loop->nstmts; s++) {
    if (pluto_vector_is_parallel(band->loop->stmts[s]->trans, d,
                                 conc_start_faces, s)) {
      return 1;
    }
  }
  return 0;
}

/*
 * Find the hyperplace parallel to the concurrent start face
 * that will be evicted; if there is none, return
 * the last hyperplane
 */
int find_hyperplane_to_be_evicted(Band *band, PlutoMatrix *conc_start_faces) {
  for (int d = band->loop->depth;
       d < (int)band->loop->depth + (int)band->width - 1; d++) {
    if (is_concurrent_start_face(band, conc_start_faces, d))
      return d;
  }
  /* Return the last one */
  return band->loop->depth + band->width - 1;
}

int get_first_non_scalar_hyperplane(PlutoProg *prog, int start, int end) {
  int i;
  for (i = start; i <= end; i++) {
    if (prog->hProps[i].type == H_LOOP)
      return i;
  }
  return -1;
}

/* Copy h2 into h1 */
static void copy_hyperplane(int64_t *h1, int64_t *h2, int ncols) {
  int j;

  for (j = 0; j < ncols; j++) {
    h1[j] = h2[j];
  }
}

/**
 * Top-level automatic transformation algoritm. Returns 0 on success.
 *
 * All dependences are reset to unsatisfied before starting.
 *
 */
int pluto_auto_transform(PlutoProg *prog) {
  int i, j, s, nsols, depth;
  /* The maximum number of linearly independent solutions needed across all
   * statements */
  int num_ind_sols_req;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  /* The number of linearly independent solutions found (max across all
   * statements) */
  int num_ind_sols_found;
  /* Pluto algo mode -- LAZY or EAGER */
  bool hyp_search_mode;

#if defined GLPK || defined GUROBI
  Graph *fcg;
  int *colour, nVertices;
  int is_skewed = false;
#endif

  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i]->satisfied = false;
  }

  /* Create the data dependence graph */
  prog->ddg = ddg_create(prog);
  ddg_compute_scc(prog);
  for (i = 0; i < prog->ddg->num_sccs; i++) {
    prog->ddg->sccs[i].vertices = NULL;
  }

  Graph *ddg = prog->ddg;
  int nvar = prog->nvar;
  int npar = prog->npar;

  prog->cst_solve_time = 0.0;
  prog->cst_const_time = 0.0;
  prog->scaling_cst_sol_time = 0.0;
  prog->mipTime = 0.0;
  prog->ilpTime = 0.0;
  prog->skew_time = 0.0;
  prog->cst_write_time = 0.0;
  prog->fcg_const_time = 0.0;
  prog->fcg_update_time = 0.0;
  prog->fcg_colour_time = 0.0;
  prog->fcg_dims_scale_time = 0.0;
  prog->fcg_cst_alloc_time = 0.0;
  prog->stencil_check_time = 0.0;

  prog->num_lp_calls = 0;

  if (nstmts == 0)
    return 0;

  PlutoMatrix **orig_trans =
      (PlutoMatrix **)malloc(nstmts * sizeof(PlutoMatrix *));
  PlutoHypType **orig_hyp_types =
      (PlutoHypType **)malloc(nstmts * sizeof(PlutoHypType *));
  int orig_num_hyperplanes = prog->num_hyperplanes;
  HyperplaneProperties *orig_hProps = prog->hProps;

  /* Get rid of any existing transformation */
  for (i = 0; i < nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    /* Save the original transformation */
    orig_trans[i] = stmt->trans;
    orig_hyp_types[i] = stmt->hyp_types;
    /* Pre-allocate a little more to prevent frequent realloc */
    stmt->trans =
        pluto_matrix_alloc(2 * stmt->dim + 1, stmt->dim + npar + 1, context);
    stmt->trans->nrows = 0;
    stmt->hyp_types = NULL;
    stmt->intra_stmt_dep_cst = NULL;
  }

  normalize_domains(prog);

  hyp_search_mode = EAGER;

  prog->num_hyperplanes = 0;
  prog->hProps = NULL;

  /* The number of independent solutions required for the deepest
   * statement */
  num_ind_sols_req = 0;
  for (i = 0; i < nstmts; i++) {
    num_ind_sols_req = PLMAX((unsigned)num_ind_sols_req, stmts[i]->dim);
  }

  depth = 0;

  if (precut(prog, ddg, depth)) {
    /* Distributed based on .fst or .precut file (customized user-supplied
     * fusion structure */
    num_ind_sols_found = pluto_get_max_ind_hyps(prog);
    printf("[pluto] Forced custom fusion structure from .fst/.precut\n");
    IF_DEBUG(
        fprintf(stdout, "%d ind solns in .precut file\n", num_ind_sols_found));
  } else {
    num_ind_sols_found = 0;
    if (options->fuse == kSmartFuse && !options->dfp) {
      cut_scc_dim_based(prog, ddg);
    }
  }

  if (options->fuse == kTypedFuse) {
    /* Mark sccs with stencils to be not distributed. This is accomplished by
     * setting that an scc already has a parallel hyperplane. This enables
     * maxfusion for all the statements in the scc. */
    double tstart = rtclock();
    for (int i = 0; i < ddg->num_sccs; i++) {
      if (is_scc_stencil(i, prog)) {
        ddg->sccs[i].is_scc_stencil = true;
        ddg->sccs[i].has_parallel_hyperplane = true;
        IF_DEBUG(printf("Scc %d has stencil dependence pattern\n", i););
      } else {
        ddg->sccs[i].is_scc_stencil = false;
        ddg->sccs[i].has_parallel_hyperplane = false;
        IF_DEBUG(printf("Scc %d has no stencil dependence patterns\n", i););
      }
    }
    prog->stencil_check_time += rtclock() - tstart;
  }

  /* For diamond tiling */
  bool conc_start_found = false;

  if (options->dfp) {
#if defined GLPK || defined GUROBI
    if (options->fuse == kNoFuse) {
      ddg_compute_scc(prog);
      cut_all_sccs(prog, ddg);
    }
    compute_scc_vertices(prog->ddg);
    IF_DEBUG(printf("[Pluto] Initial DDG\n"););
    IF_DEBUG(pluto_matrix_print(stdout, prog->ddg->adj););
    /* ddg_compute_scc(prog); */
    if (!options->silent) {
      printf("[Pluto] Building fusion conflict graph\n");
    }

    nVertices = 0;
    if (options->scc_cluster) {
      for (unsigned i = 0; i < ddg->num_sccs; i++) {
        if (ddg->sccs[i].max_dim > 0) {
          ddg->sccs[i].fcg_scc_offset = nVertices;
        }
        ddg->sccs[i].is_scc_coloured = false;
        nVertices += ddg->sccs[i].max_dim;
      }
    } else {
      for (unsigned i = 0; i < nstmts; i++) {
        if (stmts[i]->dim > 0) {
          ddg->vertices[i].fcg_stmt_offset = nVertices;
        }
        nVertices += stmts[i]->dim_orig;
      }
    }

    colour = (int *)malloc(nVertices * sizeof(int));
    for (int i = 0; i < nVertices; i++) {
      colour[i] = 0;
    }

    PlutoConstraints *permutecst = get_permutability_constraints(prog);
    IF_DEBUG(pluto_constraints_cplex_print(stdout, permutecst););

    /* The current_colour parameter has to be 1. This is the indicate that
     * colouring has to start from first level. Internally this is also used to
     * introduce must distribute edges in the fusion conflict graph. */
    prog->fcg = build_fusion_conflict_graph(prog, colour, nVertices, 1);

    fcg = prog->fcg;
    fcg->num_coloured_vertices = 0;
    fcg->to_be_rebuilt = false;

    IF_DEBUG(printf("[pluto] Fusion Conflict graph\n"););
    IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););

    prog->total_coloured_stmts = (int *)malloc(nvar * sizeof(int));
    prog->scaled_dims = (int *)malloc(nvar * sizeof(int));
    prog->coloured_dims = 0;
    for (int i = 0; i < nvar; i++) {
      prog->total_coloured_stmts[i] = 0;
      prog->scaled_dims[i] = 0;
    }

    /* This routine frees colour internally */
    find_permutable_dimensions_scc_based(colour, prog);

    if (!options->silent && options->debug) {
      printf("[pluto] Transformations before skewing \n");
      pluto_transformations_pretty_print(prog);
    }

    is_skewed = introduce_skew(prog);

    pluto_dep_satisfaction_reset(prog);
    if (is_skewed && options->diamondtile) {
      conc_start_found = pluto_diamond_tile(prog);
    }
    /* If there are any unsatisfied deps, they have to be
     * distributed at the inner most level. */
    pluto_dep_satisfaction_reset(prog);
    for (int i = 0; i < prog->num_hyperplanes; i++) {
      dep_satisfaction_update(prog, i);
    }
    if (!deps_satisfaction_check(prog)) {
      ddg_update(prog->ddg, prog);
      ddg_compute_scc(prog);
      cut_all_sccs(prog, prog->ddg);
    }

    free(prog->total_coloured_stmts);
    free(prog->scaled_dims);
#endif
  } else {

    do {
      /* Number of linearly independent solutions remaining to be found
       * (maximum across all statements) */
      int num_sols_left;

      if (options->fuse == kNoFuse) {
        ddg_compute_scc(prog);
        cut_all_sccs(prog, ddg);
      }

      num_sols_left = 0;
      for (s = 0; s < nstmts; s++) {
        /* Num linearly independent hyperplanes remaining to be
         * found for a statement; take max across all */
        num_sols_left = PLMAX(num_sols_left,
                              (int)stmts[s]->dim_orig -
                                  (int)pluto_stmt_get_num_ind_hyps(stmts[s]));
      }
      /* Progress in the EAGER mode is made every time a solution is found;
       * thus, the maximum number of linearly independent solutions
       * remaining to be found is the difference between the number required
       * for the deepest statement and the number found so far for the
       * deepest statement (since in EAGER mode, if there was a statement
       * that had fewer than num_ind_sols_found linearly independent
       * hyperplanes,
       * it means it didn't need that many hyperplanes and all of its
       * linearly independent solutions had been found */
      assert(hyp_search_mode == LAZY ||
             num_sols_left == num_ind_sols_req - num_ind_sols_found);

      if (options->per_cc_obj) {
        ddg_compute_cc(prog);
      }
      nsols = find_permutable_hyperplanes(prog, hyp_search_mode, num_sols_left,
                                          depth);

      IF_DEBUG(fprintf(stdout,
                       "[pluto] pluto_auto_transform: band level %d; "
                       "%d hyperplane(s) found\n",
                       depth, nsols));
      IF_DEBUG2(pluto_transformations_pretty_print(prog));

      num_ind_sols_found = pluto_get_max_ind_hyps(prog);

      if (nsols >= 1) {
        /* Diamond tiling: done for the first band of permutable loops */
        if (options->diamondtile && nsols >= 2 && !conc_start_found) {
          conc_start_found = pluto_diamond_tile(prog);
        }

        for (j = 0; j < nsols; j++) {
          /* Mark dependences satisfied by this solution */
          dep_satisfaction_update(prog, stmts[0]->trans->nrows - nsols + j);
          ddg_update(ddg, prog);
        }
      } else {
        /* Satisfy inter-scc dependences via distribution since we have
         * no more fusable loops */

        ddg_compute_scc(prog);

        if (get_num_unsatisfied_inter_scc_deps(prog) >= 1) {
          if (options->fuse == kNoFuse) {
            /* No fuse */
            cut_all_sccs(prog, ddg);
          } else if (options->fuse == kSmartFuse) {
            /* Smart fuse (default) */
            cut_smart(prog, ddg);
          } else {
            /* Max fuse */
            if (depth >= 2 * nvar + 1)
              cut_all_sccs(prog, ddg);
            else
              cut_conservative(prog, ddg);
          }
        } else {
          /* Only one SCC or multiple SCCs with no unsatisfied inter-SCC
           * deps, and no solutions found  */
          if (hyp_search_mode == EAGER) {
            IF_DEBUG(printf("[pluto] Switching to LAZY mode\n"););
            hyp_search_mode = LAZY;
          } else if (!deps_satisfaction_check(prog)) {
            assert(hyp_search_mode == LAZY);
            /* There is a problem; solutions should have been found if
             * there were no inter-scc deps, and some unsatisfied deps
             * existed */
            if (options->debug || options->moredebug) {
              printf("\tNumber of unsatisfied deps: %d\n",
                     get_num_unsatisfied_deps(prog->deps, prog->ndeps));
              printf("\tNumber of unsatisfied inter-scc deps: %d\n",
                     get_num_unsatisfied_inter_scc_deps(prog));
              fprintf(stdout, "[pluto] WARNING: Unfortunately, Pluto cannot "
                              "find any more hyperplanes.\n");
              fprintf(stdout, "\tThis is usually a result of (1) a bug in the "
                              "dependence tester,\n");
              fprintf(stdout,
                      "\tor (2) a bug in Pluto's auto transformation,\n");
              fprintf(stdout, "\tor (3) an inconsistent .fst/.precut in your "
                              "working directory.\n");
              fprintf(stdout, "\tTransformation found so far:\n");
              pluto_transformations_pretty_print(prog);
              pluto_compute_dep_directions(prog);
              pluto_compute_dep_satisfaction(prog);
              pluto_print_dep_directions(prog);
            }
            denormalize_domains(prog);
            printf("[pluto] WARNING: working with original (identity) "
                   "transformation (if they exist)\n");
            /* Restore original ones */
            for (i = 0; i < nstmts; i++) {
              stmts[i]->trans = orig_trans[i];
              stmts[i]->hyp_types = orig_hyp_types[i];
              prog->num_hyperplanes = orig_num_hyperplanes;
              prog->hProps = orig_hProps;
            }
            return 1;
          }
        }
      }
      /* Under LAZY mode, do a precise dep satisfaction check to take
       * care of partial satisfaction (rarely needed) */
      if (hyp_search_mode == LAZY)
        pluto_compute_dep_satisfaction_precise(prog);
      depth++;
    } while (!pluto_transformations_full_ranked(prog) ||
             !deps_satisfaction_check(prog));
  }

  /* Deallocate the fusion conflict graph */
  if (options->dfp) {
#if defined GLPK || defined GUROBI
    ddg = prog->ddg;
    graph_free(prog->fcg);
#endif
  }
  if (options->diamondtile && !conc_start_found) {
    PLUTO_MESSAGE(printf("[pluto] Diamond tiling not possible/useful\n"););
  }

  denormalize_domains(prog);

  for (i = 0; i < nstmts; i++) {
    pluto_matrix_free(orig_trans[i]);
    free(orig_hyp_types[i]);
  }
  free(orig_trans);
  free(orig_hyp_types);
  free(orig_hProps);

  IF_DEBUG(printf("[pluto] pluto_auto_transform: successful, done\n"););

  return 0;
}

int get_num_unsatisfied_deps(Dep **deps, int ndeps) {
  int i, count;

  count = 0;
  for (i = 0; i < ndeps; i++) {
    if (IS_RAR(deps[i]->type))
      continue;
    if (!deps[i]->satisfied) {
      count++;
    }
  }

  return count;
}

int get_num_unsatisfied_inter_stmt_deps(Dep **deps, int ndeps) {
  int i;

  int count = 0;
  for (i = 0; i < ndeps; i++) {
    if (IS_RAR(deps[i]->type))
      continue;
    if (deps[i]->src == deps[i]->dest)
      continue;
    if (!deps[i]->satisfied) {
      count++;
    }
  }

  return count;
}

int get_num_unsatisfied_inter_scc_deps(PlutoProg *prog) {
  int i;

  int count = 0;
  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (IS_RAR(dep->type))
      continue;
    Stmt *src_stmt = prog->stmts[dep->src];
    Stmt *dest_stmt = prog->stmts[dep->dest];
    if (src_stmt->scc_id != dest_stmt->scc_id && !dep->satisfied) {
      count++;
    }
  }

  return count;
}

void ddg_print(Graph *g) { pluto_matrix_print(stdout, g->adj); }

/*
 * Update the DDG - should be called when some dependences
 * are satisfied
 **/
void ddg_update(Graph *g, PlutoProg *prog) {
  int i, j;
  Dep *dep;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] updating DDG\n"););

  for (i = 0; i < g->nVertices; i++)
    for (j = 0; j < g->nVertices; j++)
      g->adj->val[i][j] = 0;

  for (i = 0; i < prog->ndeps; i++) {
    dep = prog->deps[i];
    if (IS_RAR(dep->type))
      continue;
    /* Number of unsatisfied dependences b/w src and dest is stored in the
     * adjacency matrix */
    g->adj->val[dep->src][dep->dest] += !dep_is_satisfied(dep);
  }
}

/*
 * Create the DDG (RAR deps not included) from the unsatisfied deps
 */
Graph *ddg_create(PlutoProg *prog) {
  int i;

  Graph *g = graph_alloc(prog->nstmts, prog->context);

  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    /* no input dep edges in the graph */
    if (IS_RAR(dep->type))
      continue;
    /* remember it's a multi-graph */
    g->adj->val[dep->src][dep->dest] += !dep_is_satisfied(dep);
  }

  return g;
}

/*
 * Get the dimensionality of the stmt with max dimensionality in the SCC
 */
static int get_max_orig_dim_in_scc(PlutoProg *prog, int scc_id) {
  int i;

  int max = -1;
  for (i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    if (stmt->scc_id == scc_id) {
      max = PLMAX(max, (int)stmt->dim_orig);
    }
  }

  return max;
}

/* Number of vertices in a given SCC */
static int get_scc_size(PlutoProg *prog, int scc_id) {
  int i;
  Stmt *stmt;

  int num = 0;
  for (i = 0; i < prog->nstmts; i++) {
    stmt = prog->stmts[i];
    if (stmt->scc_id == scc_id) {
      num++;
    }
  }

  return num;
}

/* Compute the connected components of the graph */
void ddg_compute_cc(PlutoProg *prog) {
  int i;
  int cc_id = -1;
  int num_cc = 0;
  int stmt_id;
  int time = 0;
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] ddg_compute_cc\n"););
  Graph *g = prog->ddg;
  /* Make the graph undirected. */
  Graph *gU = get_undirected_graph(g);
  for (i = 0; i < gU->nVertices; i++) {
    gU->vertices[i].vn = 0;
  }
  for (i = 0; i < gU->nVertices; i++) {
    if (gU->vertices[i].vn == 0) {
      cc_id++;
      num_cc++;
      gU->vertices[i].cc_id = cc_id;
      dfs_vertex(gU, &gU->vertices[i], &time);
      gU->vertices[i].cc_id = cc_id;
    }
    gU->vertices[i].cc_id = cc_id;
    g->vertices[i].cc_id = gU->vertices[i].cc_id;
    stmt_id = g->vertices[i].id;
    assert(stmt_id == i);

    prog->stmts[i]->cc_id = g->vertices[i].cc_id;
  }
  g->num_ccs = num_cc;
  graph_free(gU);
}

/* Compute the SCCs of a graph (using Kosaraju's algorithm) */
void ddg_compute_scc(PlutoProg *prog) {
  int i;

  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] ddg_compute_scc\n"););

  Graph *g = prog->ddg;

  dfs(g);

  Graph *gT = graph_transpose(g);

  dfs_for_scc(gT);

  g->num_sccs = gT->num_sccs;

  for (i = 0; i < g->nVertices; i++) {
    g->vertices[i].scc_id = gT->vertices[i].scc_id;
    int stmt_id = gT->vertices[i].id;
    assert(stmt_id == i);
    prog->stmts[i]->scc_id = g->vertices[i].scc_id;
  }

  for (i = 0; i < g->num_sccs; i++) {
    g->sccs[i].max_dim = get_max_orig_dim_in_scc(prog, i);
    g->sccs[i].size = get_scc_size(prog, i);
    g->sccs[i].id = gT->sccs[i].id;
    g->sccs[i].sol = NULL;
    g->sccs[i].is_parallel = 0;
  }

  graph_free(gT);

  graph_print_sccs(g);
}

/* Get this statement's schedule
 * Schedule format
 * [num sched functions | orig dim iters | params | const ]
 * Number of rows == num sched functions (each row for one hyperplane)
 */
PlutoConstraints *pluto_stmt_get_schedule(const Stmt *stmt) {
  PlutoMatrix *sched, *trans;
  PlutoConstraints *schedcst;

  trans = stmt->trans;
  sched = pluto_matrix_dup(trans);

  for (int i = 0; i < (int)sched->nrows; i++) {
    pluto_matrix_negate_row(sched, sched->nrows - 1 - i);
    pluto_matrix_add_col(sched, 0);
    sched->val[(int)trans->nrows - 1 - i][0] = 1;
  }

  schedcst = pluto_constraints_from_equalities(sched);

  pluto_matrix_free(sched);

  return schedcst;
}

/* Update a dependence with a new constraint added to the statement domain */
void pluto_update_deps(Stmt *stmt, PlutoConstraints *cst, PlutoProg *prog) {
  Stmt **stmts = prog->stmts;

  assert(cst->ncols == stmt->domain->ncols);

  for (int i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (stmts[dep->src] == stmt) {
      PlutoConstraints *cst_l = pluto_constraints_dup(cst);
      Stmt *tstmt = stmts[dep->dest];
      for (unsigned c = 0; c < tstmt->dim; c++) {
        pluto_constraints_add_dim(cst_l, stmt->dim, NULL);
      }
      pluto_constraints_add(dep->dpolytope, cst_l);
      pluto_constraints_free(cst_l);
    }
    if (stmts[dep->dest] == stmt) {
      PlutoConstraints *cst_l = pluto_constraints_dup(cst);
      Stmt *sstmt = stmts[dep->src];
      for (unsigned c = 0; c < sstmt->dim; c++) {
        pluto_constraints_add_dim(cst_l, 0, NULL);
      }
      pluto_constraints_add(dep->dpolytope, cst_l);
      pluto_constraints_free(cst_l);
    }
  }

  for (int i = 0; i < prog->ntransdeps; i++) {
    Dep *dep = prog->transdeps[i];
    if (stmts[dep->src] == stmt) {
      PlutoConstraints *cst_l = pluto_constraints_dup(cst);
      Stmt *tstmt = stmts[dep->dest];
      for (unsigned c = 0; c < tstmt->dim; c++) {
        pluto_constraints_add_dim(cst_l, stmt->dim, NULL);
      }
      pluto_constraints_add(dep->dpolytope, cst_l);
      pluto_constraints_free(cst_l);
    }
    if (stmts[dep->dest] == stmt) {
      PlutoConstraints *cst_l = pluto_constraints_dup(cst);
      Stmt *sstmt = stmts[dep->src];
      for (unsigned c = 0; c < sstmt->dim; c++) {
        pluto_constraints_add_dim(cst_l, 0, NULL);
      }
      pluto_constraints_add(dep->dpolytope, cst_l);
      pluto_constraints_free(cst_l);
    }
  }
}

/// Returns true if these statements completely fused until the innermost level.
/// This assumes one-to-one transformations for statements.
bool pluto_are_stmts_fused(Stmt **stmts, int nstmts, const PlutoProg *prog) {
  /// Find the deepest loop hyperplane for the first statement.
  int d;
  for (d = prog->num_hyperplanes - 1; d >= 0; --d) {
    if (!pluto_is_hyperplane_scalar(stmts[0], d))
      break;
  }
  if (d < 0)
    // No loop hyperplanes at all for at least one statement.
    return false;

  unsigned num;
  Ploop **loops = pluto_get_loops_under(stmts, nstmts, d, prog, &num);
  pluto_loops_free(loops, num);

  return (num == 1);
}

/// Diamond tiling. Returns true if diamond tiling transformation was applied.
bool pluto_diamond_tile(PlutoProg *prog) {
  PlutoContext *context = prog->context;

  IF_DEBUG(printf("[pluto] pluto_diamond_tile\n"));

  bool conc_start_enabled = false;

  /* Get the permutability constraints since a call to
   * detect_transformation_properties with update dep satisfaction levels
   * and we won't get the constraints we want */

  /* Don't free basecst */
  PlutoConstraints *basecst = get_permutability_constraints(prog);

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  unsigned nbands;
  Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

  for (unsigned b = 0; b < nbands; b++) {
    PlutoMatrix **cone_complement_hyps;
    Band *band = bands[b];

    /* Band should not have outer parallelism */
    if (pluto_loop_is_parallel(prog, band->loop))
      continue;

    /* Band should have inner parallelism */
    unsigned ni;
    Ploop **iloops = pluto_get_loops_immediately_inner(band->loop, prog, &ni);
    unsigned i;
    for (i = 0; i < ni; i++) {
      unsigned s;
      for (s = 0; s < band->loop->nstmts; s++) {
        if (!pluto_loop_is_parallel_for_stmt(prog, iloops[i],
                                             band->loop->stmts[s]))
          break;
      }
      if (s < band->loop->nstmts)
        break;
    }
    if (i < ni) {
      pluto_loops_free(iloops, ni);
      continue;
    }

    /* Pure Inner parallelism should be lost via tiling */
    for (i = 0; i < ni; i++) {
      if (!pluto_loop_has_satisfied_dep_with_component(prog, iloops[i]))
        break;
    }
    pluto_loops_free(iloops, ni);
    if (i < ni)
      continue;

    /* Domains should allows point-wise concurrent start */
    PlutoMatrix *conc_start_faces = get_face_with_concurrent_start(prog, band);
    if (!conc_start_faces)
      continue;

    /* face with concurrent start shouldn't be normal to all hyperplanes
     * of all statements in this band */
    unsigned s;
    for (s = 0; s < band->loop->nstmts; s++) {
      unsigned d;
      for (d = band->loop->depth; d < band->loop->depth + band->width; d++) {
        if (!pluto_vector_is_normal(band->loop->stmts[s]->trans, d,
                                    conc_start_faces, s))
          break;
      }
      if (d < band->loop->depth + band->width)
        break;
    }
    if (s == band->loop->nstmts) {
      printf("row normal\n");
      continue;
    }

    cone_complement_hyps =
        (PlutoMatrix **)malloc(band->loop->nstmts * sizeof(PlutoMatrix *));
    for (unsigned i = 0; i < band->loop->nstmts; i++) {
      cone_complement_hyps[i] = NULL;
    }

    int first_loop_hyp = band->loop->depth;
    /*
     * Find hyperplane that will be replaced by the newly found
     * hyperplane
     * Concurrent start pertains to the first band alone
     */
    int evict_pos = find_hyperplane_to_be_evicted(band, conc_start_faces);

    /* If we haven't yet found the cone_complement_pos, just
     * choose the first one as the cone_complement_pos */
    int cone_complement_pos = first_loop_hyp;

    /* If first_loop_hyp hyperplane itself is to be replaced,
     * choose the next one as cone_complement_pos */
    if (evict_pos == first_loop_hyp)
      cone_complement_pos++;
    int conc_start_enabled_band = find_cone_complement_hyperplane(
        band, conc_start_faces, evict_pos, cone_complement_pos, basecst, prog,
        cone_complement_hyps);
    pluto_matrix_free(conc_start_faces);

    /* Re-arrange the transformation matrix if concurrent start
     * was found, store the replaced hyperplane so that it can be
     * put back for the right intra-tile order */
    if (conc_start_enabled_band) {
      IF_DEBUG(
          printf("[pluto] Transformations before concurrent start enable\n"));
      IF_DEBUG(pluto_transformations_pretty_print(prog););
      IF_DEBUG(pluto_print_hyperplane_properties(prog););
      for (i = 0; i < band->loop->nstmts; i++) {
        Stmt *stmt = band->loop->stmts[i];
        /* Since we do concurrent start only once */
        assert(stmt->evicted_hyp == NULL);
        stmt->evicted_hyp = pluto_matrix_alloc(1, stmt->trans->ncols, context);
        copy_hyperplane(stmt->evicted_hyp->val[0], stmt->trans->val[evict_pos],
                        stmt->trans->ncols);
        copy_hyperplane(stmt->trans->val[evict_pos],
                        cone_complement_hyps[i]->val[0], stmt->trans->ncols);
        stmt->evicted_hyp_pos = evict_pos;
      }
    }

    conc_start_enabled |= conc_start_enabled_band;

    for (unsigned i = 0; i < band->loop->nstmts; i++) {
      pluto_matrix_free(cone_complement_hyps[i]);
    }
    free(cone_complement_hyps);
  }

  pluto_bands_free(bands, nbands);

  if (conc_start_enabled) {
    pluto_dep_satisfaction_reset(prog);
    PLUTO_MESSAGE(printf("[pluto] Concurrent start hyperplanes found\n"););
    IF_DEBUG(printf("[pluto] Transformations after concurrent start enable\n"));
    IF_DEBUG(pluto_transformations_pretty_print(prog););
  }

  return conc_start_enabled;
}
