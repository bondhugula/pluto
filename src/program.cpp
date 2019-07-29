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
 * program.c
 *
 * This file contains functions that do the job interfacing the PLUTO
 * core to the frontend and related matters
 *
 */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constraints.h"
#include "isl_support.h"
#include "math_support.h"
#include "pluto.h"
#include "program.h"

#include "osl/body.h"
#include "osl/extensions/arrays.h"
#include "osl/extensions/dependence.h"
#include "osl/extensions/loop.h"
#include "osl/extensions/pluto_unroll.h"
#include "osl/extensions/scatnames.h"
#include "osl/macros.h"
#include "osl/relation_list.h"
#include "osl/scop.h"

#include "cloog/cloog.h"

#include "candl/candl.h"
#include "candl/dependence.h"
#include "candl/options.h"
#include "candl/scop.h"

#include <isl/flow.h>
#include <isl/id.h>
#include <isl/map.h>
#include <isl/mat.h>
#include <isl/set.h>
#include <isl/space.h>
#include <isl/union_map.h>
#include <isl/val.h>

osl_relation_p get_identity_schedule(int dim, int npar);

void pluto_add_dep(PlutoProg *prog, Dep *dep) {
  dep->id = prog->ndeps;
  prog->ndeps++;
  prog->deps = (Dep **)realloc(prog->deps, sizeof(Dep *) * prog->ndeps);
  prog->deps[prog->ndeps - 1] = dep;
}

/*
 * In an [eq -I A c] relation, rows can be ordered any way.
 * Returns the index for the row for the nth output dimension.
 */
int osl_relation_get_row_id_for_nth_dimension(osl_relation_p relation,
                                              int ndim) {
  int nb_ndims_found = 0;
  int row_id = -1;

  if (relation == NULL)
    return OSL_UNDEFINED;

  if ((relation->nb_rows < ndim) || (0 > ndim)) {
    fprintf(stderr, "error: dimension out of bounds");
    exit(1);
  }

  nb_ndims_found = 0;
  for (int i = 0; i < relation->nb_rows; i++) {
    if (!osl_int_zero(relation->precision, relation->m[i][ndim])) {
      nb_ndims_found++;
      row_id = i;
    }
  }
  if (nb_ndims_found == 0) {
    fprintf(stderr, "error: specified dimension not found");
    exit(1);
  }
  if (nb_ndims_found > 1) {
    fprintf(stderr, "error: specified dimension occurs multiple times");
    exit(1);
  }

  return row_id;
}

/*
 * Converts a [eq A c] relation to [A c] Pluto constraints
 */
PlutoConstraints *osl_relation_to_pluto_constraints(osl_relation_p rln) {
  if (rln == NULL)
    return NULL;

  if (rln->nb_local_dims) {
    fprintf(stderr, "[osl_relation_to_pluto_constraints] Cannot handle Local "
                    "Dimensions in a relation.\n");
    exit(1);
  }

  PlutoConstraints *cst =
      pluto_constraints_alloc(rln->nb_rows, rln->nb_columns - 1);
  cst->nrows = rln->nb_rows;

  // copy matrix values
  for (int i = 0; i < rln->nb_rows; i++) {
    cst->is_eq[i] = osl_int_zero(rln->precision, rln->m[i][0]);
    for (unsigned j = 0; j < cst->ncols; j++) {
      cst->val[i][j] = osl_int_get_si(rln->precision, rln->m[i][j + 1]);
    }
  }

  return cst;
}

/*
 * Converts [A c] PLuto constraints to a [eq A c] domain relation
 */
osl_relation_p pluto_constraints_to_osl_domain(PlutoConstraints *cst,
                                               int npar) {
  osl_relation_p rln;

  if (cst == NULL)
    return NULL;

  rln = osl_relation_pmalloc(PLUTO_OSL_PRECISION, cst->nrows, cst->ncols + 1);

  // copy matrix values
  for (int i = 0; i < rln->nb_rows; i++) {
    osl_int_set_si(rln->precision, &rln->m[i][0], cst->is_eq[i] ? 0 : 1);
    for (unsigned j = 0; j < cst->ncols; j++) {
      osl_int_set_si(rln->precision, &rln->m[i][j + 1], cst->val[i][j]);
    }
  }

  rln->type = OSL_TYPE_DOMAIN;
  rln->nb_parameters = npar;
  rln->nb_output_dims = rln->nb_columns - rln->nb_parameters - 2;
  rln->nb_input_dims = 0;
  rln->nb_local_dims = 0;

  return rln;
}

osl_relation_p pluto_constraints_list_to_osl_domain(PlutoConstraints *cst,
                                                    int npar) {
  if (cst == NULL)
    return NULL;

  osl_relation_p list =
      osl_relation_pmalloc(PLUTO_OSL_PRECISION, cst->nrows, cst->ncols + 1);

  list = pluto_constraints_to_osl_domain(cst, npar);

  if (cst->next != NULL)
    list->next = pluto_constraints_list_to_osl_domain(cst->next, npar);

  return list;
}

/*
 * Converts a [eq -I A c] osl access relation to [A c] pluto matrix
 * Note: a[c] and a, having two and one output dimensions respectively
 * in osl, are converted to a one-dimensional pluto matrix.
 */
PlutoMatrix *osl_access_relation_to_pluto_matrix(osl_relation_p smat) {
  int i, j;

  PlutoMatrix *mat;

  if (smat == NULL)
    return NULL;

  if (smat->nb_local_dims) {
    fprintf(stderr, "Cannot handle Local Dimensions in a relation.\n");
    exit(1);
  }

  int nrows =
      smat->nb_rows == 1 ? smat->nb_rows : smat->nb_rows - 1; // skp id line
  int ncols = smat->nb_columns - smat->nb_output_dims - 1;    //-1: skip 1st col
  mat = pluto_matrix_alloc(nrows, ncols);

  // Special case for scalars.
  if (smat->nb_rows == 1) {
    for (j = smat->nb_output_dims + 1; j < smat->nb_columns; j++) {
      mat->val[0][j - (smat->nb_output_dims + 1)] = 0;
    }
  } else {
    // fill in the rest of the information
    for (i = 1; i < smat->nb_rows; i++) {
      int row = osl_relation_get_row_id_for_nth_dimension(smat, i + 1);
      for (j = smat->nb_output_dims + 1; j < smat->nb_columns; j++) {
        mat->val[i - 1][j - (smat->nb_output_dims + 1)] =
            osl_int_get_si(smat->precision, smat->m[row][j]);
      }
    }
  }

  return mat;
}

/* Return the number of lines until the next non-zero element
 * in the first column of "access" or until the end of the matrix.
 */
int access_len(PlutoMatrix *access, int first) {
  unsigned i;

  for (i = first + 1; i < access->nrows; ++i)
    if (access->val[i][0] != 0)
      break;

  return i - first;
}

int is_array(int id, PlutoMatrix *pmat) {
  int is_array = 0;

  for (unsigned i = 0; i < pmat->nrows; i++) {
    if (pmat->val[i][0] == id) {
      if (access_len(pmat, i) > 1)
        is_array = 1;
      else {
        for (unsigned j = 1; j < pmat->ncols; j++)
          if (pmat->val[i][j] != 0)
            is_array = 1;
      }
    }
  }

  return is_array;
}

/*
 * Converts a [A c] pluto matrix to [eq -I A c] osl access relation
 * Note: a[c] and a, having two and one output dimensions respectively
 * in osl, are converted to a one-dimensional pluto matrix.
 */
osl_relation_p pluto_matrix_to_osl_access_relation(PlutoMatrix *pmat) {
  int i = 0;
  int j = 0;

  if (pmat == NULL)
    return NULL;

  // check if it's a scalar
  int nrows = pmat->nrows + is_array(pmat->val[0][0], pmat);
  int ncols = 1 + nrows + pmat->ncols;

  osl_relation_p rl = osl_relation_malloc(nrows, ncols);
  // set the dims outside
  rl->nb_output_dims = nrows;
  // set the type outside

  // first row with array_id
  osl_int_set_si(rl->precision, &rl->m[0][1], -1);

  // rest of the rows
  for (i = 1; i < rl->nb_rows; i++) {
    for (j = 0; j < rl->nb_columns; j++) {
      if (j == i + 1)
        osl_int_set_si(rl->precision, &rl->m[i][j], -1);
      else if (j <= rl->nb_output_dims)
        osl_int_set_si(rl->precision, &rl->m[i][j], 0);
      else if (j < rl->nb_columns)
        osl_int_set_si(rl->precision, &rl->m[i][j],
                       pmat->val[i - 1][j - rl->nb_output_dims - 1]);
    }
  }

  return rl;
}

/*
 * Converts a [eq -I A c] osl scattering to [A c] pluto transformations
 */
PlutoMatrix *osl_scattering_to_pluto_trans(osl_relation_p smat) {
  int i, j;
  PlutoMatrix *mat;

  if (!smat)
    return NULL;

  if (smat->nb_local_dims) {
    fprintf(stderr, "Cannot handle Local Dimensions in a relation.\n");
    exit(1);
  }

  mat = pluto_matrix_alloc(smat->nb_rows,
                           smat->nb_columns - smat->nb_output_dims - 1);
  for (i = 0; i < smat->nb_rows; i++) {
    /* Only equalities in schedule expected */
    assert(osl_int_get_si(smat->precision, smat->m[i][0]) == 0);

    int row = osl_relation_get_row_id_for_nth_dimension(smat, i + 1);
    for (j = smat->nb_output_dims + 1; j < smat->nb_columns; j++) {
      mat->val[i][j - smat->nb_output_dims - 1] =
          osl_int_get_si(smat->precision, smat->m[row][j]);
    }
  }

  return mat;
}

/*
 * Converts a [A c] pluto transformations to a [eq -I A c] osl scattering
 */
osl_relation_p pluto_trans_to_osl_scattering(PlutoMatrix *mat, int npar) {
  int i, j;
  osl_relation_p smat;

  if (!mat)
    return NULL;

  smat = osl_relation_pmalloc(PLUTO_OSL_PRECISION, mat->nrows,
                              mat->nrows + mat->ncols + 1);
  smat->type = OSL_TYPE_SCATTERING;
  smat->nb_parameters = npar;
  smat->nb_output_dims = mat->nrows;
  smat->nb_input_dims = mat->ncols - npar - 1;
  smat->nb_local_dims = 0;

  for (i = 0; i < smat->nb_rows; i++) {
    for (j = 1; j < smat->nb_columns; j++) {

      /* Only equalities in schedule expected */
      if (j == 0) // eq/neq (first) column
        osl_int_set_si(smat->precision, &smat->m[i][j], 0);

      // fill out the output dims
      else if (j == i + 1)
        osl_int_set_si(smat->precision, &smat->m[i][j], -1);
      else if (j <= smat->nb_output_dims) // non diagonal zeros
        osl_int_set_si(smat->precision, &smat->m[i][j], 0);

      // fill out the intput_dims+params+const
      else
        osl_int_set_si(smat->precision, &smat->m[i][j],
                       mat->val[i][j - smat->nb_output_dims - 1]);
    }
  }

  return smat;
}

/*
 * get a list of to-be-vectorized loops from PlutoProg
 */
osl_loop_p pluto_get_vector_loop_list(const PlutoProg *prog) {
  unsigned i, j, nploops;
  osl_loop_p ret_loop = NULL;

  Ploop **ploops = pluto_get_parallel_loops(prog, &nploops);

  for (i = 0; i < nploops; i++) {
    /* Only the innermost ones */
    if (!pluto_loop_is_innermost(ploops[i], prog))
      continue;

    IF_DEBUG(printf("[pluto_get_vector_loop_list] marking loop\n"););
    IF_DEBUG(pluto_loop_print(ploops[i]););

    osl_loop_p newloop = osl_loop_malloc();

    char iter[13];
    snprintf(iter, sizeof(iter), "t%d", ploops[i]->depth + 1);
    newloop->iter = strdup(iter);

    newloop->nb_stmts = ploops[i]->nstmts;
    newloop->stmt_ids = (int *)malloc(ploops[i]->nstmts * sizeof(int));
    for (j = 0; j < ploops[i]->nstmts; j++) {
      newloop->stmt_ids[j] = ploops[i]->stmts[j]->id + 1;
    }

    newloop->directive += CLAST_PARALLEL_VEC;

    // add new loop to looplist
    osl_loop_add(newloop, &ret_loop);
  }

  pluto_loops_free(ploops, nploops);

  return ret_loop;
}

/*
 * Get a list of to-be-parallelized loops frop PlutoProg.
 */
osl_loop_p pluto_get_parallel_loop_list(const PlutoProg *prog,
                                        int vloopsfound) {
  unsigned i, j, nploops;
  osl_loop_p ret_loop = NULL;

  Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);

  IF_DEBUG(printf("[pluto_parallel_loop_list] parallelizable loops\n"););
  IF_DEBUG(pluto_loops_print(ploops, nploops););

  for (i = 0; i < nploops; i++) {
    osl_loop_p newloop = osl_loop_malloc();

    char iter[13];
    snprintf(iter, sizeof(iter), "t%d", ploops[i]->depth + 1);
    newloop->iter = strdup(iter);

    newloop->nb_stmts = ploops[i]->nstmts;
    newloop->stmt_ids = (int *)malloc(ploops[i]->nstmts * sizeof(int));
    unsigned max_depth = 0;
    for (j = 0; j < ploops[i]->nstmts; j++) {
      Stmt *stmt = ploops[i]->stmts[j];
      newloop->stmt_ids[j] = stmt->id + 1;

      if (stmt->trans->nrows > max_depth)
        max_depth = stmt->trans->nrows;
    }

    newloop->directive += CLAST_PARALLEL_OMP;
    char *private_vars = (char *)malloc(128);
    private_vars[0] = '\0';
    if (vloopsfound)
      strcpy(private_vars, "lbv, ubv");
    unsigned depth = ploops[i]->depth + 1;
    for (depth++; depth <= max_depth; depth++) {
      sprintf(private_vars + strlen(private_vars), "t%d,", depth);
    }
    if (strlen(private_vars))
      private_vars[strlen(private_vars) - 1] = '\0'; // remove last comma
    newloop->private_vars = strdup(private_vars);
    free(private_vars);

    // add new loop to looplist
    osl_loop_add(newloop, &ret_loop);
  }

  pluto_loops_free(ploops, nploops);

  return ret_loop;
}

/*
 * Replace the original scop's statements' domains and scatterings
 * by those generated by Pluto
 */
void pluto_populate_scop(osl_scop_p scop, PlutoProg *prog,
                         PlutoOptions *options) {

  int i;
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  int npar = prog->npar;

  osl_statement_p stm = scop->statement;
  /* Fill domains (may have been changed for tiling purposes). */
  for (i = 0; i < nstmts; i++) {

    // replace domain
    if (stm->domain)
      osl_relation_free(stm->domain);
    stm->domain = pluto_constraints_to_osl_domain(stmts[i]->domain, npar);

    // replace scattering
    if (stm->scattering)
      osl_relation_free(stm->scattering);
    stm->scattering = pluto_trans_to_osl_scattering(stmts[i]->trans, npar);

    stm = stm->next;
  }

  // update iterators
  // if domains(iterators) chanaged due to optimizations (tiling, etc.)
  for (stm = scop->statement; stm; stm = stm->next) {
    int niter = stm->domain->nb_columns - scop->context->nb_columns;
    int nb_orig_it = -1;
    osl_body_p stmt_body =
        (osl_body_p)osl_generic_lookup(stm->extension, OSL_URI_BODY);
    if (stmt_body) {
      nb_orig_it = osl_strings_size(stmt_body->iterators);
      if (nb_orig_it != niter) { // update iterators.

        char **iters =
            (char **)malloc(sizeof(char *) * (niter + 1)); //+1 for NULL
        for (i = 0; i < niter - nb_orig_it; ++i) {
          iters[i] = (char *)malloc(sizeof(char) * 16);
          sprintf(iters[i], "fk%d", i);

          // update accesses
          osl_relation_list_p rll = stm->access;
          while (rll) {
            osl_relation_insert_blank_column(rll->elt,
                                             rll->elt->nb_output_dims + 1);
            rll->elt->nb_input_dims++;
            rll = rll->next;
          }
        }
        for (; i < niter; ++i)
          iters[i] = stmt_body->iterators->string[i - niter + nb_orig_it];

        iters[i] = (char *)NULL;

        free(stmt_body->iterators->string);
        stmt_body->iterators->string = iters;
      }
    }
  }

  // update scatnames
  // get max scat dims
  int nb_scatt = 0;
  for (stm = scop->statement; stm; stm = stm->next) {
    int cur_scatt = stm->scattering->nb_output_dims;
    nb_scatt = nb_scatt > cur_scatt ? nb_scatt : cur_scatt;
  }

  // generate scatt names
  osl_strings_p newnames = osl_strings_generate("t", nb_scatt);
  osl_scatnames_p scatt = osl_scatnames_malloc();
  scatt->names = newnames;

  // replace the old scatnames with new one
  osl_generic_remove(&scop->extension, OSL_URI_SCATNAMES);
  osl_generic_p gen = osl_generic_shell(scatt, osl_scatnames_interface());
  osl_generic_add(&scop->extension, gen);

  // update loop information
  // get loops to be marked for parallization and vectorization
  osl_loop_p vll = NULL;
  if (options->prevector) {
    vll = pluto_get_vector_loop_list(prog);
  }
  osl_loop_p pll = NULL;
  if (options->parallel) {
    pll = pluto_get_parallel_loop_list(prog, vll != NULL);
  }
  // concatenate the two lists
  osl_loop_add(vll, &pll);

  if (pll) {
    osl_generic_p loopgen = osl_generic_shell(pll, osl_loop_interface());
    osl_generic_add(&scop->extension, loopgen);
  }

  // Add pluto_unroll extension
  {
    int i;
    HyperplaneProperties *hProps = prog->hProps;
    osl_pluto_unroll_p pluto_unroll = NULL;
    osl_pluto_unroll_p pluto_unroll_base = NULL;

    char buffer[sizeof(i) * CHAR_BIT + 1] = {0};

    for (i = 0; i < prog->num_hyperplanes; i++) {
      if (hProps[i].unroll == UNROLL || hProps[i].unroll == UNROLLJAM) {
        sprintf(buffer, "t%i", i + 1);
        if (pluto_unroll == NULL) {
          pluto_unroll = osl_pluto_unroll_malloc();
          pluto_unroll_base = pluto_unroll;
        } else {
          pluto_unroll->next = osl_pluto_unroll_malloc();
          pluto_unroll = pluto_unroll->next;
        }
        osl_pluto_unroll_fill(pluto_unroll, buffer,
                              hProps[i].unroll == UNROLLJAM, options->ufactor);
      }
    }

    if (pluto_unroll_base) {
      osl_generic_p pluto_unroll_generic =
          osl_generic_shell(pluto_unroll_base, osl_pluto_unroll_interface());
      osl_generic_add(&scop->extension, pluto_unroll_generic);
    }
  }
}

static int get_osl_write_access_position(osl_relation_list_p rl,
                                         osl_relation_p access) {
  int num;

  num = -1;

  osl_relation_list_p tmp = rl;
  for (; tmp; tmp = tmp->next) {

    if ((tmp->elt->type == OSL_TYPE_WRITE) ||
        (tmp->elt->type == OSL_TYPE_MAY_WRITE))
      num++;

    if (tmp->elt == access)
      break;
  }
  assert(num >= 0);
  return num;
}

static int get_osl_read_access_position(osl_relation_list_p rl,
                                        osl_relation_p access) {
  int num;

  num = -1;

  osl_relation_list_p tmp = rl;
  for (; tmp; tmp = tmp->next) {

    if (tmp->elt->type == OSL_TYPE_READ)
      num++;

    if (tmp->elt == access)
      break;
  }
  assert(num >= 0);
  return num;
}

/*
 * Returns a list of write or may_write access relations in a list
 */
osl_relation_list_p osl_access_list_filter_write(osl_relation_list_p list) {

  osl_relation_list_p copy = osl_relation_list_clone(list);
  osl_relation_list_p filtered = NULL;
  osl_relation_list_p previous = NULL;
  osl_relation_list_p trash;
  int first = 1;

  while (copy != NULL) {
    if ((copy->elt != NULL) && ((copy->elt->type == OSL_TYPE_WRITE) ||
                                (copy->elt->type == OSL_TYPE_MAY_WRITE))) {
      if (first) {
        filtered = copy;
        first = 0;
      }

      previous = copy;
      copy = copy->next;
    } else {
      trash = copy;
      if (!first)
        previous->next = copy->next;
      copy = copy->next;
      trash->next = NULL;
      osl_relation_list_free(trash);
    }
  }

  return filtered;
}

/*
 * Returns a list of read access relations in a list
 */
osl_relation_list_p osl_access_list_filter_read(osl_relation_list_p list) {

  osl_relation_list_p copy = osl_relation_list_clone(list);
  osl_relation_list_p filtered = NULL;
  osl_relation_list_p previous = NULL;
  osl_relation_list_p trash;
  int first = 1;

  while (copy != NULL) {
    if ((copy->elt != NULL) && (copy->elt->type == OSL_TYPE_READ)) {
      if (first) {
        filtered = copy;
        first = 0;
      }

      previous = copy;
      copy = copy->next;
    } else {
      trash = copy;
      if (!first)
        previous->next = copy->next;
      copy = copy->next;
      trash->next = NULL;
      osl_relation_list_free(trash);
    }
  }

  return filtered;
}

/*
 * Converts an osl dependence domain to Pluto constraints
 * See osl/extensions/dependence.h for the osl dependence domain matrix format
 */
PlutoConstraints *osl_dep_domain_to_pluto_constraints(osl_dependence_p in_dep) {

  int s_dom_output_dims = in_dep->source_nb_output_dims_domain;
  int t_dom_output_dims = in_dep->target_nb_output_dims_domain;

  int nb_output_dims = in_dep->source_nb_output_dims_domain +
                       in_dep->source_nb_output_dims_access;
  int nb_input_dims = in_dep->target_nb_output_dims_domain +
                      in_dep->target_nb_output_dims_access;

  /* Compute osl domain indexes */
  int osl_ind_source_local_domain = 1 + nb_output_dims + nb_input_dims;
  int osl_ind_source_local_access =
      osl_ind_source_local_domain + in_dep->source_nb_local_dims_domain;
  int osl_ind_target_local_domain =
      osl_ind_source_local_access + in_dep->source_nb_local_dims_access;
  int osl_ind_target_local_access =
      osl_ind_target_local_domain + in_dep->target_nb_local_dims_domain;
  int osl_ind_params =
      osl_ind_target_local_access + in_dep->target_nb_local_dims_access;

  /* Compute pluto constraints domain indexes */
  int pl_ind_target_domain = 1 + in_dep->source_nb_output_dims_domain;
  int pl_ind_params =
      pl_ind_target_domain + in_dep->target_nb_output_dims_domain;

  int rows, cols = 0;

  int nb_pars = in_dep->stmt_source_ptr->domain->nb_parameters;
  int s_dom_rows = in_dep->stmt_source_ptr->domain->nb_rows;
  int t_dom_rows = in_dep->stmt_target_ptr->domain->nb_rows;
  int s_acc_rows = in_dep->ref_source_access_ptr->nb_rows - 1;
  int depth = in_dep->depth;

  //
  rows = s_dom_rows + t_dom_rows +
         (s_acc_rows == 0 ? 1 : s_acc_rows) // special case for 0-dimention
                                            // array(scalar)
         + depth;
  cols = s_dom_output_dims + t_dom_output_dims + nb_pars +
         2; // cols: 2 => eq + const

  PlutoConstraints *cst;

  cst = pluto_constraints_alloc(rows, cols - 1);
  cst->nrows = rows;
  cst->ncols = cols - 1;

  int i = 0;
  int j = 0;
  int osl_constraint = 0;
  int pl_constraint = 0;
  int osl_index = 0;
  int pl_index = 0;

  // copy source domain
  osl_relation_p s_domain = in_dep->stmt_source_ptr->domain;
  for (i = 0; i < s_domain->nb_rows; i++) {

    // copy first column
    if (osl_int_zero(in_dep->domain->precision,
                     in_dep->domain->m[osl_constraint][0])) {
      cst->is_eq[pl_constraint] = 1;
    } else {
      cst->is_eq[pl_constraint] = 0;
    }

    // start of matrix
    osl_index = 1;    // start of src_stmt_domain_output_dims
    pl_index = 1 - 1; // -1 for pluto
    for (j = 0; j < s_dom_output_dims; j++)
      cst->val[pl_constraint][pl_index + j] =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_index + j]);

    // copy localdims - not supprted by converter
    if (s_domain->nb_local_dims) {
      fprintf(stderr, "local dimensions in domain not supported\n");
      exit(1);
    }

    // copy params + constant
    osl_index = osl_ind_params;
    pl_index = pl_ind_params - 1; // -1 for pluto
    for (j = 0; j < nb_pars + 1; j++)
      cst->val[pl_constraint][pl_index + j] =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_index + j]);

    osl_constraint++;
    pl_constraint++;
  }

  // copy target domain
  osl_relation_p t_domain = in_dep->stmt_target_ptr->domain;
  for (i = 0; i < t_domain->nb_rows; i++) {

    // copy first column
    if (osl_int_zero(in_dep->domain->precision,
                     in_dep->domain->m[osl_constraint][0])) {
      cst->is_eq[pl_constraint] = 1;
    } else {
      cst->is_eq[pl_constraint] = 0;
    }

    // start of matrix
    osl_index = 1 + nb_output_dims;
    pl_index = pl_ind_target_domain - 1; // -1 for pluto
    for (j = 0; j < t_dom_output_dims; j++)
      cst->val[pl_constraint][pl_index + j] =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_index + j]);

    // copy local dims - not supported in converter
    if (t_domain->nb_local_dims) {
      fprintf(stderr, "local dimensions in domain not supproted\n");
      exit(1);
    }

    // copy params + constant
    osl_index = osl_ind_params;
    pl_index = pl_ind_params - 1; // -1 for pluto
    for (j = 0; j < nb_pars + 1; j++)
      cst->val[pl_constraint][pl_index + j] =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_index + j]);

    pl_constraint++;
    osl_constraint++;
  }

  // copy source as well as target access
  int osl_s_index = 0;
  int osl_t_index = 0;
  int pl_s_index = 0;
  int pl_t_index = 0;

  osl_relation_p s_access = in_dep->ref_source_access_ptr;
  osl_relation_p t_access = in_dep->ref_target_access_ptr;

  osl_constraint++; // skip the array_id line

  for (i = 0; i < s_acc_rows; i++) {

    // copy first column
    if (osl_int_zero(in_dep->domain->precision,
                     in_dep->domain->m[osl_constraint][0])) {
      cst->is_eq[pl_constraint] = 1;
    } else {
      cst->is_eq[pl_constraint] = 0;
    }

    osl_s_index = 1;
    osl_t_index = 1 + nb_output_dims;
    pl_s_index = 1 - 1;                    // -1 for pluto
    pl_t_index = pl_ind_target_domain - 1; // -1 for pluto

    for (j = 0; j < s_access->nb_input_dims; j++) {
      cst->val[pl_constraint][pl_s_index + j] =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_s_index + j]);
    }

    for (j = 0; j < t_access->nb_input_dims; j++) { // t_acc_dims==s_acc_dims
      cst->val[pl_constraint][pl_t_index + j] = osl_int_get_si(
          in_dep->domain->precision,
          in_dep->domain
              ->m[osl_constraint + s_access->nb_rows][osl_t_index + j]);
    }

    // copy local dimensions - not supported by converter
    if (s_access->nb_local_dims || t_access->nb_local_dims) {
      fprintf(stderr, "local dimensions in Access not supproted\n");
      exit(1);
    }

    // copy params + constant
    osl_index = osl_ind_params;
    pl_index = pl_ind_params - 1; // -1 for pluto
    for (j = 0; j < nb_pars + 1; j++) {
      // get src params
      int src_param =
          osl_int_get_si(in_dep->domain->precision,
                         in_dep->domain->m[osl_constraint][osl_index + j]);
      // get tgt params
      int tgt_param = osl_int_get_si(
          in_dep->domain->precision,
          in_dep->domain->m[osl_constraint + s_access->nb_rows][osl_index + j]);

      tgt_param = -tgt_param; // oppose

      cst->val[pl_constraint][pl_index + j] = src_param - tgt_param;
    }

    pl_constraint++;
    osl_constraint++;
  }

  // copy access equalities
  // skip min_depth
  int min_depth = OSL_min(s_access->nb_output_dims, t_access->nb_output_dims);
  osl_constraint += s_access->nb_rows + min_depth;

  // s_acc_rows calculated by subtracting 1 from acc.nb_rows
  // in case of a scalar this results in 0, still add a constraint for pluto
  if (s_acc_rows == 0)
    pl_constraint++;

  // copy depth
  osl_s_index = 1;
  osl_t_index = 1 + nb_output_dims;
  pl_s_index = 1 - 1;                    // -1 for pluto
  pl_t_index = pl_ind_target_domain - 1; // -1 for pluto
  for (i = 0; i < depth; i++) {
    // copy first column
    if (osl_int_zero(in_dep->domain->precision,
                     in_dep->domain->m[osl_constraint][0])) {
      cst->is_eq[pl_constraint] = 1;
    } else {
      cst->is_eq[pl_constraint] = 0;
    }

    // copy subscript equalities
    cst->val[pl_constraint][pl_s_index + i] =
        osl_int_get_si(in_dep->domain->precision,
                       in_dep->domain->m[osl_constraint][osl_s_index + i]);
    cst->val[pl_constraint][pl_t_index + i] =
        osl_int_get_si(in_dep->domain->precision,
                       in_dep->domain->m[osl_constraint][osl_t_index + i]);

    // copy params -> not applicable here

    // copy const == last column
    cst->val[pl_constraint][cst->ncols - 1] = osl_int_get_si(
        in_dep->domain->precision,
        in_dep->domain->m[osl_constraint][in_dep->domain->nb_columns - 1]);

    osl_constraint++;
    pl_constraint++;
  }

  // return new domain
  return cst;
}

/* Read dependences from candl structures */
static Dep **deps_read(osl_dependence_p candlDeps, PlutoProg *prog) {
  int i, ndeps;
  int spos, tpos;
  Dep **deps;
  int npar = prog->npar;
  Stmt **stmts = prog->stmts;

  ndeps = osl_nb_dependences(candlDeps);

  deps = (Dep **)malloc(ndeps * sizeof(Dep *));

  for (i = 0; i < ndeps; i++) {
    deps[i] = pluto_dep_alloc();
  }

  osl_dependence_p candl_dep = candlDeps;

  candl_dep = candlDeps;

  IF_DEBUG(candl_dependence_pprint(stdout, candl_dep));

  /* Dependence polyhedra information */
  for (i = 0; i < ndeps; i++) {
    Dep *dep = deps[i];
    dep->id = i;
    dep->type = candl_dep->type;
    dep->src = candl_dep->label_source;
    dep->dest = candl_dep->label_target;

    // candl_matrix_print(stdout, candl_dep->domain);
    dep->dpolytope = osl_dep_domain_to_pluto_constraints(candl_dep);
    dep->bounding_poly = pluto_constraints_dup(dep->dpolytope);

    pluto_constraints_set_names_range(
        dep->dpolytope, stmts[dep->src]->iterators, 0, 0, stmts[dep->src]->dim);
    /* suffix the destination iterators with a '*/
    char **dnames = (char **)malloc(stmts[dep->dest]->dim * sizeof(char *));
    for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
      dnames[j] = (char *)malloc(strlen(stmts[dep->dest]->iterators[j]) + 2);
      strcpy(dnames[j], stmts[dep->dest]->iterators[j]);
      strcat(dnames[j], "'");
    }
    pluto_constraints_set_names_range(
        dep->dpolytope, dnames, stmts[dep->src]->dim, 0, stmts[dep->dest]->dim);
    for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
      free(dnames[j]);
    }
    free(dnames);

    pluto_constraints_set_names_range(
        dep->dpolytope, prog->params,
        stmts[dep->src]->dim + stmts[dep->dest]->dim, 0, npar);

    switch (dep->type) {
    case OSL_DEPENDENCE_RAW:
      spos = get_osl_write_access_position(candl_dep->stmt_source_ptr->access,
                                           candl_dep->ref_source_access_ptr);
      dep->src_acc = stmts[dep->src]->writes[spos];
      tpos = get_osl_read_access_position(candl_dep->stmt_target_ptr->access,
                                          candl_dep->ref_target_access_ptr);
      dep->dest_acc = stmts[dep->dest]->reads[tpos];

      break;
    case OSL_DEPENDENCE_WAW:
      spos = get_osl_write_access_position(candl_dep->stmt_source_ptr->access,
                                           candl_dep->ref_source_access_ptr);
      dep->src_acc = stmts[dep->src]->writes[spos];
      tpos = get_osl_write_access_position(candl_dep->stmt_target_ptr->access,
                                           candl_dep->ref_target_access_ptr);
      dep->dest_acc = stmts[dep->dest]->writes[tpos];
      break;
    case OSL_DEPENDENCE_WAR:
      spos = get_osl_read_access_position(candl_dep->stmt_source_ptr->access,
                                          candl_dep->ref_source_access_ptr);
      dep->src_acc = stmts[dep->src]->reads[spos];
      tpos = get_osl_write_access_position(candl_dep->stmt_target_ptr->access,
                                           candl_dep->ref_target_access_ptr);
      dep->dest_acc = stmts[dep->dest]->writes[tpos];
      break;
    case OSL_DEPENDENCE_RAR:
      spos = get_osl_read_access_position(candl_dep->stmt_source_ptr->access,
                                          candl_dep->ref_source_access_ptr);
      dep->src_acc = stmts[dep->src]->reads[spos];
      tpos = get_osl_read_access_position(candl_dep->stmt_target_ptr->access,
                                          candl_dep->ref_target_access_ptr);
      dep->dest_acc = stmts[dep->dest]->reads[tpos];
      break;
    default:
      assert(0);
    }

    /* Get rid of rows that are all zero */
    unsigned r, c;
    bool *remove = (bool *)malloc(sizeof(bool) * dep->dpolytope->nrows);
    for (r = 0; r < dep->dpolytope->nrows; r++) {
      for (c = 0; c < dep->dpolytope->ncols; c++) {
        if (dep->dpolytope->val[r][c] != 0) {
          break;
        }
      }
      if (c == dep->dpolytope->ncols) {
        remove[r] = true;
      } else {
        remove[r] = false;
      }
    }
    unsigned orig_nrows = dep->dpolytope->nrows;
    int del_count = 0;
    for (r = 0; r < orig_nrows; r++) {
      if (remove[r]) {
        pluto_constraints_remove_row(dep->dpolytope, r - del_count);
        del_count++;
      }
    }
    free(remove);

    int src_dim = stmts[dep->src]->dim;
    int target_dim = stmts[dep->dest]->dim;

    assert(candl_dep->source_nb_output_dims_domain +
               candl_dep->target_nb_output_dims_domain +
               candl_dep->stmt_source_ptr->domain->nb_parameters + 1 ==
           src_dim + target_dim + npar + 1);

    candl_dep = candl_dep->next;
  }

  return deps;
}

void pluto_dep_print(FILE *fp, const Dep *dep) {
  fprintf(fp,
          "--- Dep %d from S%d to S%d; satisfied: %d, sat level: %d; Type: ",
          dep->id + 1, dep->src + 1, dep->dest + 1, dep->satisfied,
          dep->satisfaction_level);

  switch (dep->type) {
  case OSL_UNDEFINED:
    fprintf(fp, "UNSET");
    break;
  case OSL_DEPENDENCE_RAW:
    fprintf(fp, "RAW");
    break;
  case OSL_DEPENDENCE_WAR:
    fprintf(fp, "WAR");
    break;
  case OSL_DEPENDENCE_WAW:
    fprintf(fp, "WAW");
    break;
  case OSL_DEPENDENCE_RAR:
    fprintf(fp, "RAR");
    break;
  default:
    fprintf(fp, "unknown");
    break;
  }

  fprintf(fp, "\n");
  if (dep->src_acc) {
    fprintf(fp, "on variable: %s\n", dep->src_acc->name);
  }

  fprintf(fp, "Dependence polyhedron\n");
  pluto_constraints_compact_print(fp, dep->dpolytope);
  fprintf(fp, "\n");
}

void pluto_deps_print(FILE *fp, PlutoProg *prog) {
  int i;
  if (prog->ndeps == 0)
    printf("** No dependences **\n\n");
  for (i = 0; i < prog->ndeps; i++) {
    pluto_dep_print(fp, prog->deps[i]);
  }
}

/* Read statement info from openscop structures (nvar: max domain dim) */
static Stmt **osl_to_pluto_stmts(const osl_scop_p scop) {
  int i, j, k;
  Stmt **stmts;
  int npar, nvar, nstmts, max_sched_rows;
  osl_statement_p scop_stmt;

  npar = scop->context->nb_parameters;
  nstmts = osl_statement_number(scop->statement);

  if (nstmts == 0)
    return NULL;

  /* Max dom dimensionality */
  nvar = -1;
  max_sched_rows = 0;
  scop_stmt = scop->statement;
  for (i = 0; i < nstmts; i++) {
    nvar = PLMAX(nvar, osl_statement_get_nb_iterators(scop_stmt));
    max_sched_rows = PLMAX(max_sched_rows, scop_stmt->scattering->nb_rows);
    scop_stmt = scop_stmt->next;
  }

  stmts = (Stmt **)malloc(nstmts * sizeof(Stmt *));

  scop_stmt = scop->statement;

  for (i = 0; i < nstmts; i++) {
    PlutoConstraints *domain =
        osl_relation_to_pluto_constraints(scop_stmt->domain);
    PlutoMatrix *trans = osl_scattering_to_pluto_trans(scop_stmt->scattering);

    int nb_iter = osl_statement_get_nb_iterators(scop_stmt);

    stmts[i] = pluto_stmt_alloc(nb_iter, domain, trans);

    /* Pad with all zero rows */
    int curr_sched_rows = stmts[i]->trans->nrows;
    for (j = curr_sched_rows; j < max_sched_rows; j++) {
      pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, j);
    }

    pluto_constraints_free(domain);
    pluto_matrix_free(trans);

    Stmt *stmt = stmts[i];

    stmt->id = i;
    stmt->type = ORIG;

    assert(scop_stmt->domain->nb_columns - 1 == (int)stmt->dim + npar + 1);

    for (unsigned j = 0; j < stmt->dim; j++) {
      stmt->is_orig_loop[j] = true;
    }

    /* Tile it if it's tilable unless turned off by .fst/.precut file */
    stmt->tile = 1;

    osl_body_p stmt_body =
        (osl_body_p)osl_generic_lookup(scop_stmt->extension, OSL_URI_BODY);

    for (unsigned j = 0; j < stmt->dim; j++) {
      stmt->iterators[j] = strdup(stmt_body->iterators->string[j]);
    }
    /* Set names for domain dimensions */
    char **names = (char **)malloc((stmt->domain->ncols - 1) * sizeof(char *));
    for (unsigned k = 0; k < stmt->dim; k++) {
      names[k] = stmt->iterators[k];
    }
    osl_strings_p osl_scop_params = NULL;
    if (scop->context->nb_parameters) {
      osl_scop_params = (osl_strings_p)scop->parameters->data;
      for (k = 0; k < npar; k++) {
        names[stmt->dim + k] = osl_scop_params->string[k];
      }
    }
    pluto_constraints_set_names(stmt->domain, names);
    free(names);

    /* Statement text */
    stmt->text = osl_strings_sprint(stmt_body->expression); // appends \n
    stmt->text[strlen(stmt->text) - 1] = '\0'; // remove the \n from end

    /* Read/write accesses */
    osl_relation_list_p wlist = osl_access_list_filter_write(scop_stmt->access);
    osl_relation_list_p rlist = osl_access_list_filter_read(scop_stmt->access);

    osl_relation_list_p rlist_t, wlist_t;
    rlist_t = rlist;
    wlist_t = wlist;

    stmt->nwrites = osl_relation_list_count(wlist);
    stmt->writes =
        (PlutoAccess **)malloc(stmt->nwrites * sizeof(PlutoAccess *));

    stmt->nreads = osl_relation_list_count(rlist);
    stmt->reads = (PlutoAccess **)malloc(stmt->nreads * sizeof(PlutoAccess *));

    osl_arrays_p arrays =
        (osl_arrays_p)osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);

    int count = 0;
    while (wlist != NULL) {
      PlutoMatrix *wmat = osl_access_relation_to_pluto_matrix(wlist->elt);
      stmt->writes[count] = (PlutoAccess *)malloc(sizeof(PlutoAccess));
      stmt->writes[count]->mat = wmat;

      // stmt->writes[count]->symbol = NULL;
      if (arrays) {
        int id = osl_relation_get_array_id(wlist->elt);
        stmt->writes[count]->name = strdup(arrays->names[id - 1]);
      } else {
        stmt->writes[count]->name = NULL;
      }

      count++;
      wlist = wlist->next;
    }

    count = 0;
    while (rlist != NULL) {

      PlutoMatrix *rmat = osl_access_relation_to_pluto_matrix(rlist->elt);
      stmt->reads[count] = (PlutoAccess *)malloc(sizeof(PlutoAccess));
      stmt->reads[count]->mat = rmat;

      // stmt->reads[count]->symbol = NULL;
      if (arrays) {
        int id = osl_relation_get_array_id(rlist->elt);
        stmt->reads[count]->name = strdup(arrays->names[id - 1]);
      } else {
        stmt->reads[count]->name = NULL;
      }

      count++;
      rlist = rlist->next;
    }

    osl_relation_list_free(wlist_t);
    osl_relation_list_free(rlist_t);

    scop_stmt = scop_stmt->next;
  }

  return stmts;
}

void pluto_access_print(FILE *fp, const PlutoAccess *acc, const Stmt *stmt) {
  int j, npar;

  if (!acc) {
    fprintf(fp, "access is NULL\n");
    return;
  }

  npar = acc->mat->ncols - (int)stmt->dim - 1;

  fprintf(fp, "%s", acc->name);
  for (unsigned i = 0; i < acc->mat->nrows; i++) {
    fprintf(fp, "[");
    const char **vars =
        (const char **)malloc((stmt->dim + npar) * sizeof(char *));
    for (unsigned j = 0; j < stmt->dim; j++) {
      vars[j] = stmt->iterators[j];
    }
    for (j = 0; j < npar; j++) {
      if (stmt->domain->names && stmt->domain->names[stmt->dim + j]) {
        vars[stmt->dim + j] = stmt->domain->names[stmt->dim + j];
      } else {
        vars[stmt->dim + j] = "p?";
      }
    }
    pluto_affine_function_print(stdout, acc->mat->val[i], stmt->dim + npar,
                                vars);
    fprintf(fp, "]");
    free(vars);
  }
  fprintf(fp, "\n");
}

void pluto_stmt_print(FILE *fp, const Stmt *stmt) {
  fprintf(fp, "S%d \"%s\"\n", stmt->id + 1, stmt->text);

  fprintf(fp, "ndims: %d; orig_depth: %d\n", stmt->dim, stmt->dim_orig);
  fprintf(fp, "Index set\n");
  pluto_constraints_compact_print(fp, stmt->domain);
  pluto_stmt_transformation_print(stmt);

  if (stmt->nreads == 0) {
    fprintf(fp, "No Read accesses\n");
  } else {
    fprintf(fp, "Read accesses\n");
    for (int i = 0; i < stmt->nreads; i++) {
      pluto_access_print(fp, stmt->reads[i], stmt);
    }
  }

  if (stmt->nwrites == 0) {
    fprintf(fp, "No write access\n");
  } else {
    fprintf(fp, "Write accesses\n");
    for (int i = 0; i < stmt->nwrites; i++) {
      pluto_access_print(fp, stmt->writes[i], stmt);
    }
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    printf("Original loop: %d -> %s\n", i,
           stmt->is_orig_loop[i] ? "yes" : "no");
  }

  fprintf(fp, "\n");
}

/*
 * Checks whether the transformation at 'level' corresponds to a tiled
 * dimension (an affine function of domain dimension divided by a constant),
 * and if yes, returns the function that was tiled (as a function of the
 * domain iterators/param) aka tiling hyperplane, tile_size will point to the
 * tile size / divisor; if no, returns NULL
 *
 * t_i = z is a tile dimension if there exist the following two constraints
 * in the domain:
 *
 * T*z <= f(i_S, p) <= T*z + T - 1, where T is the tile size / divisor
 *
 * Thus, t_i = f(i_S, p)/T
 *
 * pos: position of the supernode in the domain
 *
 */
static int64_t *pluto_check_supernode(const Stmt *stmt, unsigned pos,
                                      int *tile_size) {
  int lb_pos, ub_pos;
  int64_t *tile_hyp;

  PlutoConstraints *dom = stmt->domain;

  lb_pos = -1;
  ub_pos = -1;
  *tile_size = -1;

  for (unsigned r = 0; r < dom->nrows; r++) {
    if (dom->val[r][pos] >= 1 &&
        dom->val[r][dom->ncols - 1] == dom->val[r][pos] - 1) {
      ub_pos = r;
      *tile_size = dom->val[r][pos];
    }

    if (dom->val[r][pos] <= -1 && dom->val[r][dom->ncols - 1] == 0) {
      lb_pos = r;
    }
  }

  if (ub_pos == -1 || lb_pos == -1)
    return NULL;

  int c;
  for (c = 0; c < (int)dom->ncols - 1; c++) {
    if (dom->val[ub_pos][c] != -dom->val[lb_pos][c])
      break;
  }
  if (c < (int)dom->ncols - 1)
    return NULL;

  tile_hyp = (int64_t *)malloc(dom->ncols * sizeof(int64_t));

  for (unsigned c = 0; c < dom->ncols; c++) {
    if (c == pos)
      tile_hyp[c] = 0;
    else
      tile_hyp[c] = dom->val[lb_pos][c];
  }

  return tile_hyp;
}

static int is_skewed(int64_t *func, int len) {
  int count, i;

  count = 0;

  for (i = 0; i < len; i++) {
    if (func[i] != 0)
      count++;
  }

  return count <= 1 ? 0 : 1;
}

/*
 * Prints the 1-d affine function/hyperplane for 'stmt' at depth 'level'
 */
void pluto_stmt_print_hyperplane(FILE *fp, const Stmt *stmt, int level) {
  int npar, j;

  npar = stmt->domain->ncols - (int)stmt->dim - 1;

  char **vars = (char **)malloc((stmt->dim + npar) * sizeof(char *));

  for (unsigned j = 0; j < stmt->dim; j++) {
    vars[j] = strdup(stmt->iterators[j]);
  }
  for (j = 0; j < npar; j++) {
    if (stmt->domain->names && stmt->domain->names[stmt->dim + j]) {
      vars[stmt->dim + j] = stmt->domain->names[stmt->dim + j];
    } else {
      vars[stmt->dim + j] = (char *)"p?";
    }
  }

  for (unsigned j = 0; j < stmt->dim; j++) {
    /* Detect if this dimension is an affine function of other dimensions
     * divided by a constant -- useful to print tiled hyperplanes, the
     * dividing constant being the tile size */
    int div;
    int64_t *super_func;
    super_func = pluto_check_supernode(stmt, j, &div);
    if (super_func) {
      char *tmp;
      tmp = pluto_affine_function_sprint(super_func, stmt->dim + npar,
                                         (const char **)vars);
      free(vars[j]);
      vars[j] = tmp;
      if (is_skewed(super_func, stmt->domain->ncols)) {
        vars[j] = (char *)realloc(vars[j],
                                  1 + strlen(vars[j]) + 2 + log10(div) + 1 + 1);
        sprintf(vars[j], "(%s", tmp = strdup(vars[j]));
        free(tmp);
        sprintf(vars[j] + strlen(vars[j]), ")/%d", div);
      } else {
        vars[j] =
            (char *)realloc(vars[j], strlen(vars[j]) + 1 + log10(div) + 1 + 1);
        sprintf(vars[j] + strlen(vars[j]), "/%d", div);
      }
      free(super_func);
    }
  }

  pluto_affine_function_print(fp, stmt->trans->val[level], stmt->dim + npar,
                              (const char **)vars);

  for (unsigned j = 0; j < stmt->dim; j++) {
    free(vars[j]);
  }

  free(vars);
}

void pluto_stmts_print(FILE *fp, Stmt **stmts, int nstmts) {
  int i;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_print(fp, stmts[i]);
  }
}

void pluto_prog_print(FILE *fp, PlutoProg *prog) {
  int i;

  fprintf(fp, "nvar = %d, npar = %d\n", prog->nvar, prog->npar);
  fprintf(fp, "Parameters: ");

  for (i = 0; i < prog->npar; i++) {
    fprintf(fp, "%s ", prog->params[i]);
  }
  fprintf(fp, "\n");

  pluto_stmts_print(fp, prog->stmts, prog->nstmts);
  pluto_deps_print(fp, prog);
  pluto_transformations_pretty_print(prog);
}

void pluto_dep_free(Dep *dep) {
  pluto_constraints_free(dep->dpolytope);
  pluto_constraints_free(dep->bounding_poly);
  pluto_constraints_free(dep->depsat_poly);
  if (dep->dirvec) {
    free(dep->dirvec);
  }
  if (dep->dirvec) {
    free(dep->satvec);
  }
  pluto_constraints_free(dep->cst);
  pluto_constraints_free(dep->bounding_cst);
  free(dep);
}

/* Set the dimension names of type "type" according to the elements
 * in the array "names".
 */
static __isl_give isl_space *set_names(__isl_take isl_space *space,
                                       enum isl_dim_type type, char **names) {
  int i;

  for (i = 0; i < isl_space_dim(space, type); ++i)
    space = isl_space_set_dim_name(space, type, i, names[i]);

  return space;
}

/* Convert a osl_relation_p containing the constraints of a domain
 * to an isl_set.
 * One shot only; does not take into account the next ptr.
 */
static __isl_give isl_set *osl_relation_to_isl_set(osl_relation_p relation,
                                                   __isl_take isl_space *dim) {
  int i, j;
  int n_eq = 0, n_ineq = 0;
  isl_ctx *ctx;
  isl_mat *eq, *ineq;
  isl_basic_set *bset;

  ctx = isl_space_get_ctx(dim);

  for (i = 0; i < relation->nb_rows; ++i)
    if (osl_int_zero(relation->precision, relation->m[i][0]))
      n_eq++;
    else
      n_ineq++;

  eq = isl_mat_alloc(ctx, n_eq, relation->nb_columns - 1);
  ineq = isl_mat_alloc(ctx, n_ineq, relation->nb_columns - 1);

  n_eq = n_ineq = 0;
  for (i = 0; i < relation->nb_rows; ++i) {
    isl_mat **m;
    int row;

    if (osl_int_zero(relation->precision, relation->m[i][0])) {
      m = &eq;
      row = n_eq++;
    } else {
      m = &ineq;
      row = n_ineq++;
    }

    for (j = 0; j < relation->nb_columns - 1; ++j) {
      int t = osl_int_get_si(relation->precision, relation->m[i][1 + j]);
      *m = isl_mat_set_element_si(*m, row, j, t);
    }
  }

  bset = isl_basic_set_from_constraint_matrices(
      dim, eq, ineq, isl_dim_set, isl_dim_div, isl_dim_param, isl_dim_cst);
  return isl_set_from_basic_set(bset);
}

/* Convert a osl_relation_p describing a union of domains
 * to an isl_set.
 */
static __isl_give isl_set *
osl_relation_list_to_isl_set(osl_relation_p list, __isl_take isl_space *space) {
  isl_set *set;

  set = isl_set_empty(isl_space_copy(space));
  for (; list; list = list->next) {
    isl_set *set_i;
    set_i = osl_relation_to_isl_set(list, isl_space_copy(space));
    set = isl_set_union(set, set_i);
  }

  isl_space_free(space);
  return set;
}

/* Convert an m x ( n + 1) pluto access_matrix_p [d A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *pluto_extract_equalities(isl_ctx *ctx,
                                                    PlutoMatrix *matrix) {
  int i, j;
  int n_col, n;
  isl_mat *eq;

  n_col = matrix->ncols;
  n = matrix->nrows;

  eq = isl_mat_alloc(ctx, n, n + n_col);

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j)
      eq = isl_mat_set_element_si(eq, i, j, 0);
    eq = isl_mat_set_element_si(eq, i, i, -1);
    for (j = 0; j < n_col; ++j) {
      eq = isl_mat_set_element_si(eq, i, n + j, matrix->val[i][j]);
    }
  }

  return eq;
}

/* Convert an m x (1 + m + n + 1) osl_relation_p [d -I A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *extract_equalities_osl(isl_ctx *ctx,
                                                  osl_relation_p relation) {
  int i, j;
  int n_col, n_row;
  isl_mat *eq;

  n_col = relation->nb_columns;
  n_row = relation->nb_rows;

  eq = isl_mat_alloc(ctx, n_row, n_col - 1);

  for (i = 0; i < n_row; ++i) {
    for (j = 0; j < n_col - 1; ++j) {
      int row = osl_relation_get_row_id_for_nth_dimension(relation, i + 1);
      int t = osl_int_get_si(relation->precision, relation->m[row][1 + j]);
      isl_val *v = isl_val_int_from_si(ctx, t);
      eq = isl_mat_set_element_val(eq, i, j, v);
    }
  }

  return eq;
}

/* Convert an m x (1 + m + n + 1) osl_relation_p [d -I A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *
extract_equalities_osl_access(isl_ctx *ctx, osl_relation_p relation) {
  int i, j;
  int n_col, n_row;
  isl_mat *eq;

  n_row = relation->nb_rows == 1 ? 1 : relation->nb_rows - 1;
  n_col = relation->nb_columns - (relation->nb_rows == 1 ? 1 : 2);

  eq = isl_mat_alloc(ctx, n_row, n_col);

  if (relation->nb_rows == 1) {
    isl_val *v = isl_val_negone(ctx);
    eq = isl_mat_set_element_val(eq, 0, 0, v);
    for (j = 1; j < n_col; ++j) {
      v = isl_val_zero(ctx);
      eq = isl_mat_set_element_val(eq, 0, j, v);
    }
  } else {
    for (i = 1; i < relation->nb_rows; ++i) {
      for (j = 2; j < relation->nb_columns; ++j) {
        int row = osl_relation_get_row_id_for_nth_dimension(relation, i + 1);
        int t = osl_int_get_si(relation->precision, relation->m[row][j]);
        isl_val *v = isl_val_int_from_si(ctx, t);
        eq = isl_mat_set_element_val(eq, i - 1, j - 2, v);
      }
    }
  }

  return eq;
}

/* Convert an m x (1 + n + 1) scoplib_matrix_p [d A c]
 * to an m x (m + n + 1) isl_mat [-I A c].
 */
static __isl_give isl_mat *extract_equalities(isl_ctx *ctx, PlutoMatrix *matrix,
                                              int first, int n) {
  int i, j;
  int n_col;
  isl_mat *eq;

  n_col = matrix->ncols;

  eq = isl_mat_alloc(ctx, n, n + n_col);

  for (i = 0; i < n; ++i) {
    isl_val *v = isl_val_zero(ctx);
    for (j = 0; j < n; ++j)
      eq = isl_mat_set_element_val(eq, i, j, v);
    eq = isl_mat_set_element_val(eq, i, i, isl_val_negone(ctx));
    for (j = 0; j < n_col - 1; ++j) {
      int t = matrix->val[first + i][j];
      v = isl_val_int_from_si(ctx, t);
      eq = isl_mat_set_element_val(eq, i, n + j, v);
    }
  }

  return eq;
}

/* Convert a pluto matrix schedule [ A c] to
 * the isl_map { i -> A i + c } in the space prescribed by "dim".
 */
static __isl_give isl_map *
pluto_matrix_schedule_to_isl_map(PlutoMatrix *schedule,
                                 __isl_take isl_space *dim) {
  int n_row, n_col;
  isl_ctx *ctx;
  isl_mat *eq, *ineq;
  isl_basic_map *bmap;

  ctx = isl_space_get_ctx(dim);
  n_row = schedule->nrows;
  n_col = schedule->ncols;

  ineq = isl_mat_alloc(ctx, 0, n_row + n_col);
  eq = extract_equalities(ctx, schedule, 0, n_row);

  bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq, isl_dim_out,
                                                isl_dim_in, isl_dim_div,
                                                isl_dim_param, isl_dim_cst);
  return isl_map_from_basic_map(bmap);
}

/* Convert a osl_relation_p scattering [0 M A c] to
 * the isl_map { i -> A i + c } in the space prescribed by "dim".
 */
static __isl_give isl_map *
osl_scattering_to_isl_map(osl_relation_p scattering,
                          __isl_take isl_space *dim) {
  int n_col;
  isl_ctx *ctx;
  isl_mat *eq, *ineq;
  isl_basic_map *bmap;

  ctx = isl_space_get_ctx(dim);
  n_col = scattering->nb_columns;

  ineq = isl_mat_alloc(ctx, 0, n_col - 1);
  eq = extract_equalities_osl(ctx, scattering);

  bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq, isl_dim_out,
                                                isl_dim_in, isl_dim_div,
                                                isl_dim_param, isl_dim_cst);

  return isl_map_from_basic_map(bmap);
}

/* Convert a osl_relation_list_p describing a series of accesses [eq -I B c]
 * to an isl_union_map with domain "dom" (in space "D").
 * The -I columns identify the output dimensions of the access, the first
 * of them being the identity of the array being accessed.  The remaining
 * output dimensions identiy the array subscripts.
 *
 * Let "A" be array identified by the first entry.
 * The input dimension columns have the form [B c].
 * Each such access is converted to a map { D[i] -> A[B i + c] } * dom.
 *
 */
static __isl_give isl_union_map *
osl_access_list_to_isl_union_map(osl_relation_list_p list,
                                 __isl_take isl_set *dom, char **arrays) {
  int len, n_col;
  isl_ctx *ctx;
  isl_space *space;
  isl_mat *eq, *ineq;
  isl_union_map *res;

  ctx = isl_set_get_ctx(dom);

  space = isl_set_get_space(dom);
  space = isl_space_drop_dims(space, isl_dim_set, 0,
                              isl_space_dim(space, isl_dim_set));
  res = isl_union_map_empty(space);

  for (; list; list = list->next) {

    n_col = list->elt->nb_columns - (list->elt->nb_rows == 1 ? 1 : 2);
    len = list->elt->nb_rows == 1 ? 1 : list->elt->nb_rows - 1;

    isl_basic_map *bmap;
    isl_map *map;
    int arr = osl_relation_get_array_id(list->elt) - 1;

    space = isl_set_get_space(dom);
    space = isl_space_from_domain(space);
    space = isl_space_add_dims(space, isl_dim_out, len);
    space = isl_space_set_tuple_name(space, isl_dim_out, arrays[arr]);

    ineq = isl_mat_alloc(ctx, 0, n_col);
    eq = extract_equalities_osl_access(ctx, list->elt);

    bmap = isl_basic_map_from_constraint_matrices(space, eq, ineq, isl_dim_out,
                                                  isl_dim_in, isl_dim_div,
                                                  isl_dim_param, isl_dim_cst);
    map = isl_map_from_basic_map(bmap);
    map = isl_map_intersect_domain(map, isl_set_copy(dom));
    res = isl_union_map_union(res, isl_union_map_from_map(map));
  }

  isl_set_free(dom);

  return res;
}

/*
 * Like osl_access_list_to_isl_union_map, but just for a single osl access
 * (read or write)
 */
static __isl_give isl_map *
osl_basic_access_to_isl_union_map(osl_relation_p access,
                                  __isl_take isl_set *dom, char **arrays) {
  int len, n_col;
  isl_ctx *ctx;
  isl_space *dim;
  isl_mat *eq, *ineq;

  ctx = isl_set_get_ctx(dom);

  n_col = access->nb_columns - (access->nb_rows == 1 ? 1 : 2);
  len = access->nb_rows == 1 ? 1 : access->nb_rows - 1;

  isl_basic_map *bmap;
  isl_map *map;
  int arr = osl_relation_get_array_id(access) - 1;

  dim = isl_set_get_space(dom);
  dim = isl_space_from_domain(dim);
  dim = isl_space_add_dims(dim, isl_dim_out, len);
  dim = isl_space_set_tuple_name(dim, isl_dim_out, arrays[arr]);

  ineq = isl_mat_alloc(ctx, 0, n_col);
  eq = extract_equalities_osl_access(ctx, access);

  bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq, isl_dim_out,
                                                isl_dim_in, isl_dim_div,
                                                isl_dim_param, isl_dim_cst);
  map = isl_map_from_basic_map(bmap);
  map = isl_map_intersect_domain(map, dom);

  return map;
}

/*
 * Like osl_access_list_to_isl_union_map, but just for a single pluto access
 * (read or write)
 * pos: position (starting row) of the access in 'access'
 */
static __isl_give isl_map *
pluto_basic_access_to_isl_union_map(PlutoMatrix *mat, char *access_name,
                                    __isl_take isl_set *dom) {
  int len, n_col;
  isl_ctx *ctx;
  isl_space *dim;
  isl_mat *eq, *ineq;

  ctx = isl_set_get_ctx(dom);

  dim = isl_set_get_space(dom);
  dim =
      isl_space_drop_dims(dim, isl_dim_set, 0, isl_space_dim(dim, isl_dim_set));

  n_col = mat->ncols;

  isl_basic_map *bmap;
  isl_map *map;
  // int arr = SCOPVAL_get_si(access->p[pos][0]) - 1;

  len = mat->nrows;

  dim = isl_set_get_space(dom);
  dim = isl_space_from_domain(dim);
  dim = isl_space_add_dims(dim, isl_dim_out, len);
  dim = isl_space_set_tuple_name(dim, isl_dim_out, access_name);

  ineq = isl_mat_alloc(ctx, 0, len + n_col);
  eq = pluto_extract_equalities(ctx, mat);

  bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq, isl_dim_out,
                                                isl_dim_in, isl_dim_div,
                                                isl_dim_param, isl_dim_cst);
  map = isl_map_from_basic_map(bmap);
  map = isl_map_intersect_domain(map, dom);

  return map;
}

/* Temporary data structure used inside extract_deps.
 *
 * deps points to the array of Deps being constructed
 * type is the type of the next Dep
 * index is the index of the next Dep in the array.
 */
struct pluto_extra_dep_info {
  Dep **deps;
  Stmt **stmts;
  int type;
  int index;
};

/* Convert an isl_basic_map describing part of a dependence to a Dep.
 * The names of the input and output spaces are of the form S_d or S_d_e
 * with d an integer identifying the statement, e identifying the access
 * (relative to the statement). If it's of the form S_d_e and read/write
 * accesses for the statement are available, source and target accesses
 * are set for the dependence, otherwise not.
 *
 * isl divs are removed; so this is an over-approximation in some cases
 */
static isl_stat basic_map_extract_dep(__isl_take isl_basic_map *bmap,
                                      void *user) {
  Stmt **stmts;
  Dep *dep;
  struct pluto_extra_dep_info *info;

  info = (struct pluto_extra_dep_info *)user;

  stmts = info->stmts;

  bmap = isl_basic_map_remove_divs(bmap);

  dep = info->deps[info->index];

  dep->id = info->index;

  dep->dirvec = NULL;
  dep->type = info->type;
  dep->src = atoi(isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2);

  /* The range space can be wrapped, if we didn't use lastwriter.  For
   * example, [T, N] -> { S_1_w0[t, i] -> [S_0_r2[t', i'] -> a[o2]] : i'
   * = -1 + i and o2 = i and t >= 0 and i >= 3 and i <= -2 + N and t' > t
   * and t' < T } */
  isl_space *space = isl_basic_map_get_space(bmap);
  if (isl_space_range_is_wrapping(space)) {
    bmap = isl_basic_map_range_factor_domain(bmap);
    isl_space_free(space);
    space = isl_basic_map_get_space(bmap);
  }
  dep->dest = atoi(isl_space_get_tuple_name(space, isl_dim_out) + 2);
  dep->dpolytope = isl_basic_map_to_pluto_constraints(bmap);
  dep->bounding_poly = pluto_constraints_dup(dep->dpolytope);
  isl_space_free(space);

  /* Inconsistent dependence if this assertion fails */
  assert(dep->dpolytope->ncols == stmts[dep->src]->dim + stmts[dep->dest]->dim +
                                      stmts[dep->src]->domain->ncols -
                                      stmts[dep->src]->dim);

  pluto_constraints_set_names_range(dep->dpolytope, stmts[dep->src]->iterators,
                                    0, 0, stmts[dep->src]->dim);

  /* Suffix the destination iterators with a '*/
  char **dnames = (char **)malloc(stmts[dep->dest]->dim * sizeof(char *));
  for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
    dnames[j] = (char *)malloc(strlen(stmts[dep->dest]->iterators[j]) + 2);
    strcpy(dnames[j], stmts[dep->dest]->iterators[j]);
    strcat(dnames[j], "'");
  }
  pluto_constraints_set_names_range(
      dep->dpolytope, dnames, stmts[dep->src]->dim, 0, stmts[dep->dest]->dim);
  for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
    free(dnames[j]);
  }
  free(dnames);

  /* parameters */
  pluto_constraints_set_names_range(
      dep->dpolytope, stmts[dep->dest]->domain->names,
      stmts[dep->src]->dim + stmts[dep->dest]->dim, stmts[dep->dest]->dim,
      stmts[dep->dest]->domain->ncols - stmts[dep->dest]->dim - 1);

  if (options->isldepaccesswise &&
      (stmts[dep->src]->reads != NULL && stmts[dep->dest]->reads != NULL)) {
    /* Extract access function information */
    int src_acc_num, dest_acc_num;
    char src_type, dest_type;
    const char *src_name, *dest_name;
    src_name = isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2;
    while (*src_name != '\0' && *(src_name++) != '_')
      ;
    if (*src_name != '\0') {
      src_type = *src_name;
      src_acc_num = atoi(src_name + 1);
    } else
      assert(0); // access function num not encoded in dependence

    dest_name = isl_basic_map_get_tuple_name(bmap, isl_dim_out) + 2;
    while (*dest_name != '\0' && *(dest_name++) != '_')
      ;
    if (*dest_name != '\0') {
      dest_type = *dest_name;
      dest_acc_num = atoi(dest_name + 1);
    } else
      assert(0); // access function num not encoded in dependence

    switch (info->type) {
    case OSL_DEPENDENCE_RAW:
      dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
      break;
    case OSL_DEPENDENCE_WAW:
      dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
      break;
    /*
     * Sometimes dep_war from isl are not WAR deps; there are WAW deps
     * included in the may deps and we can't assume that they are all WAR may
     * deps. Mark them correctly.
     */
    case OSL_DEPENDENCE_WAR:
      dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
      if (src_type == 'w' && dest_type == 'w') {
        /* Fix the type */
        dep->type = OSL_DEPENDENCE_WAW;
        /* This is really a WAW dep */
        dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      } else {
        dep->src_acc = stmts[dep->src]->reads[src_acc_num];
      }
      break;
    case OSL_DEPENDENCE_RAR:
      dep->src_acc = stmts[dep->src]->reads[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
      break;
    default:
      assert(0);
    }
  } else {
    dep->src_acc = NULL;
    dep->dest_acc = NULL;
  }

  info->index++;
  isl_basic_map_free(bmap);
  return isl_stat_ok;
}

/* Extract Pluto dependences from an isl_map */
static isl_stat map_extract_dep(__isl_take isl_map *map, void *user) {
  isl_stat r = isl_map_foreach_basic_map(map, &basic_map_extract_dep, user);
  isl_map_free(map);
  return r;
}

/// Extract a Pluto access function from an access encoded in an isl_basic_map.
static isl_stat
isl_basic_map_extract_access_func(__isl_take isl_basic_map *bmap, void *user) {
  isl_map *map = isl_map_from_basic_map(bmap);

  int dim = isl_map_dim(map, isl_dim_out);
  int ncols =
      isl_map_dim(map, isl_dim_in) + isl_map_dim(map, isl_dim_param) + 1;

  PlutoMatrix *acc_mat = pluto_matrix_alloc(0, ncols);

  for (int i = 0; i < dim; i++) {
    PlutoMatrix *func_onedim = NULL;
    if (isl_map_dim_is_single_valued(map, i)) {
      isl_pw_aff *pw_aff = isl_pw_aff_from_map_dim(map, i);
      /* Best effort: Gets it from the last piece */
      isl_pw_aff_foreach_piece(pw_aff, isl_aff_to_pluto_func, &func_onedim);
      pluto_matrix_add(acc_mat, func_onedim);
      pluto_matrix_free(func_onedim);
      isl_pw_aff_free(pw_aff);
    } else {
      // Multi-valued map for access function? Setting this to zero.
      pluto_matrix_add_row(acc_mat, 0);
      pluto_matrix_zero_row(acc_mat, 0);
    }
  }
  struct pluto_access_meta_info *info = (struct pluto_access_meta_info *)user;
  (*info->accs)[info->index] = (PlutoAccess *)malloc(sizeof(PlutoAccess));
  PlutoAccess *acc = (*info->accs)[info->index];
  acc->name = strdup(isl_basic_map_get_tuple_name(bmap, isl_dim_out));
  acc->mat = acc_mat;
  info->index++;
  isl_map_free(map);
  return isl_stat_ok;
}

/// Extract a Pluto access function from an isl_map.
isl_stat isl_map_extract_access_func(__isl_take isl_map *map, void *user) {
  /* Extract a PlutoAccess from every isl_basic_map */
  isl_stat r =
      isl_map_foreach_basic_map(map, &isl_basic_map_extract_access_func, user);
  isl_map_free(map);
  return r;
}

/// Extract read and writes accesses encoded in isl_union_map structures and
/// populate Pluto stmt's accesses. The access relation maps statement domain
/// instances to elements in a data space. The name on the output tuple provides
/// array name.
void extract_accesses_for_pluto_stmt(Stmt *stmt,
                                     __isl_keep isl_union_map *reads,
                                     __isl_keep isl_union_map *writes) {
  isl_union_map_foreach_map(reads, &isl_map_count, &stmt->nreads);
  isl_union_map_foreach_map(writes, &isl_map_count, &stmt->nwrites);

  int npar = stmt->domain->ncols - stmt->dim - 1;

  struct pluto_access_meta_info e_reads = {&stmt->reads, 0, stmt->dim, npar};
  struct pluto_access_meta_info e_writes = {&stmt->writes, 0, stmt->dim, npar};

  if (stmt->nreads > 0) {
    stmt->reads = (PlutoAccess **)malloc(stmt->nreads * sizeof(PlutoAccess *));
  }
  if (stmt->nwrites > 0) {
    stmt->writes =
        (PlutoAccess **)malloc(stmt->nwrites * sizeof(PlutoAccess *));
  }
  for (int j = 0; j < stmt->nreads; j++) {
    stmt->reads[j] = NULL;
  }
  for (int j = 0; j < stmt->nwrites; j++) {
    stmt->writes[j] = NULL;
  }

  isl_union_map_foreach_map(reads, &isl_map_extract_access_func, &e_reads);
  isl_union_map_foreach_map(writes, &isl_map_extract_access_func, &e_writes);

  // Verify correct extraction.
  for (int j = 0; j < stmt->nreads; j++) {
    assert(stmt->reads[j]->mat->ncols == stmt->dim + npar + 1 &&
           "Inconsistent access matrix");
  }
  for (int j = 0; j < stmt->nwrites; j++) {
    // Verify correct extraction.
    assert(stmt->writes[j]->mat->ncols == stmt->dim + npar + 1 &&
           "Inconsistent access matrix");
  }
}

/* Extract deps from isl union maps into Pluto Deps */
int extract_deps(Dep **deps, int first, Stmt **stmts,
                 __isl_keep isl_union_map *umap, int type) {
  struct pluto_extra_dep_info info = {deps, stmts, type, first};

  isl_union_map_foreach_map(umap, &map_extract_dep, &info);

  return info.index - first;
}

osl_names_p get_scop_names(osl_scop_p scop) {

  // generate temp names
  osl_names_p names = osl_scop_names(scop);

  // if scop has names substitute them for temp names
  if (scop->context->nb_parameters) {
    osl_strings_free(names->parameters);
    names->parameters =
        osl_strings_clone((osl_strings_p)scop->parameters->data);
  }

  osl_arrays_p arrays =
      (osl_arrays_p)osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
  if (arrays) {
    osl_strings_free(names->arrays);
    names->arrays = osl_arrays_to_strings(arrays);
  }

  return names;
}

// Compute dependences using ISL.
// If options->lastwriter is false, then
//       RAW deps are those from any earlier write to a read
//       WAW deps are those from any earlier write to a write
//       WAR deps are those from any earlier read to a write
//       RAR deps are those from any earlier read to a read
//  If options->lastwriter is true, then
//       RAW deps are those from the last write to a read
//       WAW deps are those from the last write to a write
//       WAR deps are those from any earlier read not masked by an intermediate
//       write to a write
//       RAR deps are those from the last read to a read
//
//  The RAR deps are only computed if options->rar is set.
void compute_deps_isl(__isl_keep isl_union_map *reads,
                      __isl_keep isl_union_map *writes,
                      __isl_keep isl_union_map *schedule,
                      __isl_keep isl_union_map *empty, isl_union_map **dep_raw,
                      isl_union_map **dep_war, isl_union_map **dep_waw,
                      isl_union_map **dep_rar, isl_union_map **trans_dep_war,
                      isl_union_map **trans_dep_waw) {
  assert(options && "options not set");

  if (options->lastwriter) {
    // Compute RAW dependences with last writer (no transitive dependences).
    isl_union_map_compute_flow(
        isl_union_map_copy(reads), isl_union_map_copy(writes),
        isl_union_map_copy(empty), isl_union_map_copy(schedule), dep_raw, NULL,
        NULL, NULL);
    // Compute WAW and WAR dependences without transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(writes),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), dep_waw,
        dep_war, NULL, NULL);
    // Compute WAR dependences with transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL,
        trans_dep_war, NULL, NULL);
    // Compute WAW dependences with transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL,
        trans_dep_waw, NULL, NULL);
    if (options->rar) {
      // Compute RAR dependences without transitive dependences.
      isl_union_map_compute_flow(
          isl_union_map_copy(reads), isl_union_map_copy(reads),
          isl_union_map_copy(empty), isl_union_map_copy(schedule), dep_rar,
          NULL, NULL, NULL);
    }
  } else {
    // Without lastwriter, compute transitive dependences.
    // RAW dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(reads), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL, dep_raw,
        NULL, NULL);
    // WAR dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL, dep_war,
        NULL, NULL);
    // WAW dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL, dep_waw,
        NULL, NULL);
    if (options->rar) {
      // RAR dependences.
      isl_union_map_compute_flow(
          isl_union_map_copy(reads), isl_union_map_copy(empty),
          isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL,
          dep_rar, NULL, NULL);
    }
  }

  if (options->isldepcoalesce) {
    *dep_raw = isl_union_map_coalesce(*dep_raw);
    *dep_war = isl_union_map_coalesce(*dep_war);
    *dep_waw = isl_union_map_coalesce(*dep_waw);

    if (options->lastwriter) {
      *trans_dep_war = isl_union_map_coalesce(*trans_dep_war);
      *trans_dep_waw = isl_union_map_coalesce(*trans_dep_waw);
    }
  }
}

/* Compute dependences based on the iteration domain and access
 * information in "scop" and put the result in "prog".
 */
static void compute_deps(osl_scop_p scop, PlutoProg *prog,
                         PlutoOptions *options) {
  int i, racc_num, wacc_num;
  int nstmts = osl_statement_number(scop->statement);
  isl_space *space;
  isl_space *param_space;
  isl_set *context;
  isl_union_map *dep_raw, *dep_war, *dep_waw, *dep_rar, *trans_dep_war;
  isl_union_map *trans_dep_waw;
  osl_statement_p stmt;
  osl_strings_p scop_params = NULL;

  if (!options->silent) {
    printf("[pluto] compute_deps (isl%s)\n",
           options->lastwriter ? " with lastwriter" : "");
  }

  isl_ctx *ctx = isl_ctx_alloc();
  assert(ctx);

  osl_names_p names = get_scop_names(scop);

  space = isl_space_set_alloc(ctx, scop->context->nb_parameters, 0);
  if (scop->context->nb_parameters) {
    scop_params = (osl_strings_p)scop->parameters->data;
    space = set_names(space, isl_dim_param, scop_params->string);
  }
  param_space = isl_space_params(isl_space_copy(space));
  context = osl_relation_to_isl_set(scop->context, param_space);

  if (!options->rar)
    dep_rar = isl_union_map_empty(isl_space_copy(space));
  isl_union_map *writes = isl_union_map_empty(isl_space_copy(space));
  isl_union_map *reads = isl_union_map_empty(isl_space_copy(space));
  isl_union_map *schedule = isl_union_map_empty(isl_space_copy(space));
  isl_union_map *empty = isl_union_map_empty(space);

  if (!options->isldepaccesswise) {
    /* Leads to fewer dependences. Each dependence may not have a unique
     * source/target access relating to it, since a union is taken
     * across all reads for a statement (and writes) for a particualr
     * array. Relationship between a dependence and associated dependent
     * data / array elements is lost, and some analyses may not work with
     * such a representation
     */
    for (i = 0, stmt = scop->statement; i < nstmts; ++i, stmt = stmt->next) {
      isl_set *dom;
      isl_map *schedule_i;
      isl_union_map *read_i;
      isl_union_map *write_i;
      char name[20];

      snprintf(name, sizeof(name), "S_%d", i);

      int niter = osl_statement_get_nb_iterators(stmt);
      space = isl_space_set_alloc(ctx, scop->context->nb_parameters, niter);
      if (scop->context->nb_parameters) {
        scop_params = (osl_strings_p)scop->parameters->data;
        space = set_names(space, isl_dim_param, scop_params->string);
      }
      if (niter) {
        osl_body_p stmt_body =
            (osl_body_p)osl_generic_lookup(stmt->extension, OSL_URI_BODY);
        space = set_names(space, isl_dim_set, stmt_body->iterators->string);
      }
      space = isl_space_set_tuple_name(space, isl_dim_set, name);
      dom = osl_relation_list_to_isl_set(stmt->domain, space);
      dom = isl_set_intersect_params(dom, isl_set_copy(context));

      space = isl_space_alloc(ctx, scop->context->nb_parameters, niter,
                              2 * niter + 1);
      if (scop->context->nb_parameters) {
        scop_params = (osl_strings_p)scop->parameters->data;
        space = set_names(space, isl_dim_param, scop_params->string);
      }
      if (niter) {
        osl_body_p stmt_body =
            (osl_body_p)osl_generic_lookup(stmt->extension, OSL_URI_BODY);
        space = set_names(space, isl_dim_in, stmt_body->iterators->string);
      }
      space = isl_space_set_tuple_name(space, isl_dim_in, name);
      schedule_i = osl_scattering_to_isl_map(stmt->scattering, space);

      osl_relation_list_p rlist = osl_access_list_filter_read(stmt->access);
      osl_relation_list_p wlist = osl_access_list_filter_write(stmt->access);

      osl_arrays_p arrays =
          (osl_arrays_p)osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
      if (arrays) {
        osl_strings_free(names->arrays);
        names->arrays = osl_arrays_to_strings(arrays);
      }

      read_i = osl_access_list_to_isl_union_map(rlist, isl_set_copy(dom),
                                                names->arrays->string);
      write_i = osl_access_list_to_isl_union_map(wlist, isl_set_copy(dom),
                                                 names->arrays->string);

      reads = isl_union_map_union(reads, read_i);
      writes = isl_union_map_union(writes, write_i);
      schedule =
          isl_union_map_union(schedule, isl_union_map_from_map(schedule_i));

      osl_relation_list_free(rlist);
      osl_relation_list_free(wlist);
    }
  } else {
    /* Each dependence is for a particular source and target access. Use
     * <stmt, access> pair while relating to accessed data so each
     * dependence can be associated to a unique source and target access
     */

    for (i = 0, stmt = scop->statement; i < nstmts; ++i, stmt = stmt->next) {
      isl_set *dom;

      racc_num = 0;
      wacc_num = 0;

      osl_relation_list_p access = stmt->access;
      for (; access; access = access->next) {
        isl_map *read_pos;
        isl_map *write_pos;
        isl_map *schedule_i;

        char name[25];

        if (access->elt->type == OSL_TYPE_READ) {
          snprintf(name, sizeof(name), "S_%d_r%d", i, racc_num);
        } else {
          snprintf(name, sizeof(name), "S_%d_w%d", i, wacc_num);
        }

        int niter = osl_statement_get_nb_iterators(stmt);
        space = isl_space_set_alloc(ctx, scop->context->nb_parameters, niter);
        if (scop->context->nb_parameters) {
          scop_params = (osl_strings_p)scop->parameters->data;
          space = set_names(space, isl_dim_param, scop_params->string);

          osl_strings_free(names->parameters);
          names->parameters = osl_strings_clone(scop_params);
        }
        if (niter) {
          osl_body_p stmt_body =
              (osl_body_p)osl_generic_lookup(stmt->extension, OSL_URI_BODY);
          space = set_names(space, isl_dim_set, stmt_body->iterators->string);

          osl_strings_free(names->iterators);
          names->iterators = osl_strings_clone(stmt_body->iterators);
        }
        space = isl_space_set_tuple_name(space, isl_dim_set, name);
        dom = osl_relation_list_to_isl_set(stmt->domain, space);
        dom = isl_set_intersect_params(dom, isl_set_copy(context));

        space = isl_space_alloc(ctx, scop->context->nb_parameters, niter,
                                2 * niter + 1);
        if (scop->context->nb_parameters) {
          scop_params = (osl_strings_p)scop->parameters->data;
          space = set_names(space, isl_dim_param, scop_params->string);
        }
        if (niter) {
          osl_body_p stmt_body =
              (osl_body_p)osl_generic_lookup(stmt->extension, OSL_URI_BODY);
          space = set_names(space, isl_dim_in, stmt_body->iterators->string);
        }
        space = isl_space_set_tuple_name(space, isl_dim_in, name);

        schedule_i = osl_scattering_to_isl_map(stmt->scattering, space);

        osl_arrays_p arrays =
            (osl_arrays_p)osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
        if (arrays) {
          osl_strings_free(names->arrays);
          names->arrays = osl_arrays_to_strings(arrays);
        }

        if (access->elt->type == OSL_TYPE_READ) {
          read_pos = osl_basic_access_to_isl_union_map(access->elt, dom,
                                                       names->arrays->string);
          reads = isl_union_map_union(reads, isl_union_map_from_map(read_pos));
        } else {
          write_pos = osl_basic_access_to_isl_union_map(access->elt, dom,
                                                        names->arrays->string);
          writes =
              isl_union_map_union(writes, isl_union_map_from_map(write_pos));
        }

        schedule =
            isl_union_map_union(schedule, isl_union_map_from_map(schedule_i));

        if (access->elt->type == OSL_TYPE_READ) {
          racc_num++;
        } else {
          wacc_num++;
        }
      }
    }
  }

  compute_deps_isl(reads, writes, schedule, empty, &dep_raw, &dep_war, &dep_waw,
                   &dep_rar, &trans_dep_war, &trans_dep_waw);

  prog->ndeps = 0;
  isl_union_map_foreach_map(dep_raw, &isl_map_count, &prog->ndeps);
  isl_union_map_foreach_map(dep_war, &isl_map_count, &prog->ndeps);
  isl_union_map_foreach_map(dep_waw, &isl_map_count, &prog->ndeps);
  isl_union_map_foreach_map(dep_rar, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  prog->ndeps = 0;
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_raw,
                              OSL_DEPENDENCE_RAW);
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_war,
                              OSL_DEPENDENCE_WAR);
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_waw,
                              OSL_DEPENDENCE_WAW);
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_rar,
                              OSL_DEPENDENCE_RAR);

  if (options->lastwriter) {
    prog->ntransdeps = 0;
    isl_union_map_foreach_map(dep_raw, &isl_map_count, &prog->ntransdeps);
    isl_union_map_foreach_map(trans_dep_war, &isl_map_count, &prog->ntransdeps);
    isl_union_map_foreach_map(trans_dep_waw, &isl_map_count, &prog->ntransdeps);
    isl_union_map_foreach_map(dep_rar, &isl_map_count, &prog->ntransdeps);

    if (prog->ntransdeps >= 1) {
      prog->transdeps = (Dep **)malloc(prog->ntransdeps * sizeof(Dep *));
      for (i = 0; i < prog->ntransdeps; i++) {
        prog->transdeps[i] = pluto_dep_alloc();
      }
      prog->ntransdeps = 0;
      prog->ntransdeps +=
          extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, dep_raw,
                       OSL_DEPENDENCE_RAW);
      prog->ntransdeps +=
          extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts,
                       trans_dep_war, OSL_DEPENDENCE_WAR);
      prog->ntransdeps +=
          extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts,
                       trans_dep_waw, OSL_DEPENDENCE_WAW);
      prog->ntransdeps +=
          extract_deps(prog->transdeps, prog->ntransdeps, prog->stmts, dep_rar,
                       OSL_DEPENDENCE_RAR);
    }

    isl_union_map_free(trans_dep_war);
    isl_union_map_free(trans_dep_waw);
  }

  isl_union_map_free(dep_raw);
  isl_union_map_free(dep_war);
  isl_union_map_free(dep_waw);
  isl_union_map_free(dep_rar);

  isl_union_map_free(writes);
  isl_union_map_free(reads);
  isl_union_map_free(schedule);
  isl_union_map_free(empty);
  isl_set_free(context);

  if (names)
    osl_names_free(names);

  isl_ctx_free(ctx);
}

PlutoMatrix *get_identity_schedule_new(int dim, int npar) {
  PlutoMatrix *smat = pluto_matrix_alloc(2 * dim + 1, dim + npar + 1);

  int i, j;
  for (i = 0; i < 2 * dim + 1; i++)
    for (j = 0; j < dim + 1 + npar; j++)
      smat->val[i][j] = 0;

  for (i = 0; i < dim; i++)
    smat->val[i][i] = 1;

  return smat;
}

int read_codegen_context_from_file(PlutoConstraints *codegen_context) {
  FILE *fp = fopen("codegen.context", "r");

  if (fp) {
    IF_DEBUG(printf("[Pluto] Reading from codegen.context\n"););
    PlutoConstraints *cc = pluto_constraints_read(fp);
    if (cc && cc->ncols == codegen_context->ncols) {
      pluto_constraints_add(codegen_context, cc);
      return 0;
    }
    IF_DEBUG(printf("[WARNING] Failed to read from codegen.context\n"););
  }

  return 1;
}

/*
 * Extract necessary information from clan_scop to create PlutoProg - a
 * representation of the program sufficient to be used throughout Pluto.
 * PlutoProg also includes dependences; so candl is run here.
 */
PlutoProg *scop_to_pluto_prog(osl_scop_p scop, PlutoOptions *options) {
  int i, max_sched_rows, npar;

  PlutoProg *prog = pluto_prog_alloc();

  /* Program parameters */
  npar = scop->context->nb_parameters;

  osl_strings_p osl_scop_params = NULL;
  if (npar >= 1)
    osl_scop_params = (osl_strings_p)scop->parameters->data;

  for (i = 0; i < npar; i++) {
    pluto_prog_add_param(prog, osl_scop_params->string[i], prog->npar);
  }

  pluto_constraints_free(prog->context);
  prog->context = osl_relation_to_pluto_constraints(scop->context);

  if (options->codegen_context != -1) {
    for (i = 0; i < prog->npar; i++) {
      pluto_constraints_add_inequality(prog->codegen_context);
      prog->codegen_context->val[i][i] = 1;
      prog->codegen_context->val[i][prog->codegen_context->ncols - 1] =
          -options->codegen_context;
    }
  }
  read_codegen_context_from_file(prog->codegen_context);

  prog->nstmts = osl_statement_number(scop->statement);
  prog->options = options;

  /* Data variables in the program */
  osl_arrays_p arrays =
      (osl_arrays_p)osl_generic_lookup(scop->extension, OSL_URI_ARRAYS);
  if (arrays == NULL) {
    prog->num_data = 0;
    fprintf(stderr, "warning: arrays extension not found\n");
  } else {
    prog->num_data = arrays->nb_names;
    prog->data_names = (char **)malloc(prog->num_data * sizeof(char *));
    for (i = 0; i < prog->num_data; i++) {
      prog->data_names[i] = strdup(arrays->names[i]);
    }
  }

  osl_statement_p scop_stmt = scop->statement;

  prog->nvar = osl_statement_get_nb_iterators(scop_stmt);
  max_sched_rows = 0;
  for (i = 0; i < prog->nstmts; i++) {
    int stmt_num_iter = osl_statement_get_nb_iterators(scop_stmt);
    prog->nvar = PLMAX(prog->nvar, stmt_num_iter);
    max_sched_rows = PLMAX(max_sched_rows, scop_stmt->scattering->nb_rows);
    scop_stmt = scop_stmt->next;
  }

  prog->stmts = osl_to_pluto_stmts(scop);
  prog->scop = scop;

  /* Compute dependences */
  if (options->isldep) {
    compute_deps(scop, prog, options);
  } else {
    /*  Using Candl */
    candl_options_p candlOptions = candl_options_malloc();
    if (options->rar) {
      candlOptions->rar = 1;
    }
    /* No longer supported */
    candlOptions->lastwriter = options->lastwriter;
    candlOptions->scalar_privatization = options->scalpriv;
    // candlOptions->verbose = 1;

    /* Add more infos (depth, label, ...) */
    /* Needed by Candl */
    candl_scop_usr_init(scop);

    osl_dependence_p candl_deps = candl_dependence(scop, candlOptions);
    prog->deps = deps_read(candl_deps, prog);
    prog->ndeps = osl_nb_dependences(candl_deps);
    candl_options_free(candlOptions);
    osl_dependence_free(candl_deps);

    candl_scop_usr_cleanup(scop); // undo candl_scop_user_init

    prog->transdeps = NULL;
    prog->ntransdeps = 0;
  }

  /* Add hyperplanes */
  if (prog->nstmts >= 1) {
    for (i = 0; i < max_sched_rows; i++) {
      pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_UNKNOWN);
      prog->hProps[prog->num_hyperplanes - 1].type =
          (i % 2) ? H_LOOP : H_SCALAR;
    }
  }

  /* Hack for linearized accesses */
  FILE *lfp = fopen(".linearized", "r");
  FILE *nlfp = fopen(".nonlinearized", "r");
  char tmpstr[256];
  char linearized[256];
  if (lfp && nlfp) {
    for (i = 0; i < prog->nstmts; i++) {
      rewind(lfp);
      rewind(nlfp);
      while (!feof(lfp) && !feof(nlfp)) {
        fgets(tmpstr, 256, nlfp);
        fgets(linearized, 256, lfp);
        if (strstr(tmpstr, prog->stmts[i]->text)) {
          prog->stmts[i]->text = (char *)realloc(
              prog->stmts[i]->text, sizeof(char) * (strlen(linearized) + 1));
          strcpy(prog->stmts[i]->text, linearized);
        }
      }
    }
    fclose(lfp);
    fclose(nlfp);
  }

  return prog;
}

/* Get an upper bound for transformation coefficients to prevent spurious
 * transformations that represent shifts or skews proportional to trip counts:
 * this happens when loop bounds are constants
 */
int pluto_prog_get_largest_const_in_domains(const PlutoProg *prog) {
  int max = 0;
  for (int i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    for (unsigned r = 0; r < stmt->domain->nrows; r++) {
      max = PLMAX(max, stmt->domain->val[r][stmt->domain->ncols - 1]);
    }
  }

  return max - 1;
}

PlutoProg *pluto_prog_alloc() {
  PlutoProg *prog = (PlutoProg *)malloc(sizeof(PlutoProg));

  prog->nstmts = 0;
  prog->stmts = NULL;
  prog->npar = 0;
  prog->nvar = 0;
  prog->params = NULL;
  prog->context = pluto_constraints_alloc(1, prog->npar + 1);
  prog->codegen_context = pluto_constraints_alloc(1, prog->npar + 1);
  prog->deps = NULL;
  prog->ndeps = 0;
  prog->transdeps = NULL;
  prog->ntransdeps = 0;
  prog->ddg = NULL;
  prog->fcg = NULL;
  prog->hProps = NULL;
  prog->num_hyperplanes = 0;
  prog->decls = (char *)malloc(16384 * 9);
  prog->data_names = NULL;
  prog->num_data = 0;

  strcpy(prog->decls, "");

  prog->globcst = NULL;

  prog->num_parameterized_loops = -1;

  return prog;
}

void pluto_prog_free(PlutoProg *prog) {
  int i;

  /* Free dependences */
  for (i = 0; i < prog->ndeps; i++) {
    pluto_dep_free(prog->deps[i]);
  }
  free(prog->deps);

  for (i = 0; i < prog->ntransdeps; i++) {
    pluto_dep_free(prog->transdeps[i]);
  }
  free(prog->transdeps);

  /* Free DDG */
  if (prog->ddg != NULL) {
    graph_free(prog->ddg);
  }

  if (prog->hProps != NULL) {
    free(prog->hProps);
  }

  for (i = 0; i < prog->npar; i++) {
    free(prog->params[i]);
  }
  if (prog->npar >= 1) {
    free(prog->params);
  }

  /* Statements */
  for (i = 0; i < prog->nstmts; i++) {
    pluto_stmt_free(prog->stmts[i]);
  }
  if (prog->nstmts >= 1) {
    free(prog->stmts);
  }

  pluto_constraints_free(prog->context);
  pluto_constraints_free(prog->codegen_context);

  pluto_constraints_free(prog->globcst);

  free(prog->decls);

  for (i = 0; i < prog->num_data; i++) {
    free(prog->data_names[i]);
  }
  free(prog->data_names);

  free(prog);
}

PlutoOptions *pluto_options_alloc() {
  PlutoOptions *options;

  options = (PlutoOptions *)malloc(sizeof(PlutoOptions));

  /* Initialize to default */
  options->flic = 0;
  options->tile = 1;
  options->intratileopt = 1;
  options->dynschedule = 0;
  options->dynschedule_graph = 0;
  options->dynschedule_graph_old = 0;
  options->dyn_trans_deps_tasks = 0;
  options->debug = 0;
  options->moredebug = 0;
  options->scancount = 0;
  options->parallel = 1;
  options->innerpar = 0;
  options->identity = 0;

  options->pet = 0;

  /* Enable one dimension of concurrent startup by default */
  options->diamondtile = 1;
  options->fulldiamondtile = 0;

  options->per_cc_obj = 0;

  options->iss = 0;
  options->unrolljam = 0;

  /* Unroll/jam factor */
  options->ufactor = 8;

  /* Ignore input deps */
  options->rar = 0;

  /* Override for first and last levels to tile */
  options->ft = -1;
  options->lt = -1;

  /* Override for first and last cloog options */
  options->cloogf = -1;
  options->cloogl = -1;

  options->cloogsh = 0;

  options->cloogbacktrack = 1;

  options->multipar = 0;
  options->l2tile = 0;
  options->prevector = 1;
  options->fuse = kSmartFuse;

  /* Experimental */
  options->delayed_cut = 0;
  options->hybridcut = 0;

  /* Default context is no context */
  options->codegen_context = -1;

  options->coeff_bound = -1;

  options->forceparallel = 0;

  options->bee = 0;

  options->isldep = 0;
  options->isldepaccesswise = 1;
  /* Disabled due to a potential bug in coalescing. Reproduce with
   * examples/heat-2d/heat-2d.c - coalescing dep_raw leads to no hyperplanes
   * being found. */
  options->isldepcoalesce = 0;

  options->candldep = 0;

  options->pipsolve = 0;
  options->islsolve = 1;
  options->glpk = 0;
  options->gurobi = 0;

  options->lp = 0;
  options->dfp = 0;
  options->ilp = 0;

  options->lpcolour = 0;
  options->scc_cluster = 0;

  options->readscop = 0;

  options->lastwriter = 0;

  options->nodepbound = 0;

  options->scalpriv = 0;

  options->silent = 0;

  options->out_file = NULL;

  options->time = 1;

  return options;
}

/* Add global/program parameter at position 'pos' */
void pluto_prog_add_param(PlutoProg *prog, const char *param, int pos) {
  int i, j;

  for (i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    pluto_constraints_add_dim(
        stmt->domain, stmt->domain->ncols - 1 - prog->npar + pos, param);
    pluto_matrix_add_col(stmt->trans,
                         stmt->trans->ncols - 1 - prog->npar + pos);

    for (j = 0; j < stmt->nwrites; j++) {
      pluto_matrix_add_col(stmt->writes[j]->mat, stmt->dim + pos);
    }
    for (j = 0; j < stmt->nreads; j++) {
      pluto_matrix_add_col(stmt->reads[j]->mat, stmt->dim + pos);
    }
  }
  for (i = 0; i < prog->ndeps; i++) {
    pluto_constraints_add_dim(
        prog->deps[i]->dpolytope,
        prog->deps[i]->dpolytope->ncols - 1 - prog->npar + pos, NULL);
  }
  pluto_constraints_add_dim(prog->context,
                            prog->context->ncols - 1 - prog->npar + pos, param);
  pluto_constraints_add_dim(prog->codegen_context,
                            prog->codegen_context->ncols - 1 - prog->npar + pos,
                            param);

  prog->params =
      (char **)realloc(prog->params, sizeof(char *) * (prog->npar + 1));

  for (i = prog->npar - 1; i >= pos; i--) {
    prog->params[i + 1] = prog->params[i];
  }

  prog->params[pos] = strdup(param);
  prog->npar++;
}

void pluto_options_free(PlutoOptions *options) {
  if (options->out_file != NULL) {
    free(options->out_file);
  }
  free(options);
}

/* pos: position of domain iterator
 * time_pos: position of time iterator; iter: domain iterator; supply -1
 * if you don't want a scattering function row added for it */
void pluto_stmt_add_dim(Stmt *stmt, unsigned pos, int time_pos,
                        const char *iter, PlutoHypType hyp_type,
                        PlutoProg *prog) {
  int i, npar;

  npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(pos <= stmt->dim);
  assert(time_pos <= (int)stmt->trans->nrows);
  assert(stmt->dim + npar + 1 == stmt->domain->ncols);

  pluto_constraints_add_dim(stmt->domain, pos, NULL);
  stmt->dim++;
  stmt->iterators =
      (char **)realloc(stmt->iterators, stmt->dim * sizeof(char *));
  for (i = stmt->dim - 2; i >= (int)pos; i--) {
    stmt->iterators[i + 1] = stmt->iterators[i];
  }
  stmt->iterators[pos] = strdup(iter);

  /* Stmt should always have a transformation */
  assert(stmt->trans != NULL);
  pluto_matrix_add_col(stmt->trans, pos);

  if (time_pos != -1) {
    pluto_matrix_add_row(stmt->trans, time_pos);
    stmt->trans->val[time_pos][pos] = 1;

    stmt->hyp_types = (PlutoHypType *)realloc(
        stmt->hyp_types, sizeof(PlutoHypType) * stmt->trans->nrows);
    for (i = stmt->trans->nrows - 2; i >= time_pos; i--) {
      stmt->hyp_types[i + 1] = stmt->hyp_types[i];
    }
    stmt->hyp_types[time_pos] = hyp_type;
  }

  /* Update is_orig_loop */
  stmt->is_orig_loop =
      (bool *)realloc(stmt->is_orig_loop, sizeof(bool) * stmt->dim);
  for (i = stmt->dim - 2; i >= (int)pos; i--) {
    stmt->is_orig_loop[i + 1] = stmt->is_orig_loop[i];
  }
  stmt->is_orig_loop[pos] = true;

  for (i = 0; i < stmt->nwrites; i++) {
    pluto_matrix_add_col(stmt->writes[i]->mat, pos);
  }
  for (i = 0; i < stmt->nreads; i++) {
    pluto_matrix_add_col(stmt->reads[i]->mat, pos);
  }

  /* Update dependences */
  for (i = 0; i < prog->ndeps; i++) {
    if (prog->deps[i]->src == stmt->id) {
      pluto_constraints_add_dim(prog->deps[i]->dpolytope, pos, NULL);
      pluto_constraints_add_dim(prog->deps[i]->bounding_poly, pos, NULL);
    }
    if (prog->deps[i]->dest == stmt->id) {
      pluto_constraints_add_dim(prog->deps[i]->dpolytope,
                                prog->stmts[prog->deps[i]->src]->dim + pos,
                                NULL);
      pluto_constraints_add_dim(prog->deps[i]->bounding_poly,
                                prog->stmts[prog->deps[i]->src]->dim + pos,
                                NULL);
    }
  }

  for (i = 0; i < prog->ntransdeps; i++) {
    assert(prog->transdeps[i] != NULL);
    if (prog->transdeps[i]->src == stmt->id) {
      pluto_constraints_add_dim(prog->transdeps[i]->dpolytope, pos, NULL);
      pluto_constraints_add_dim(prog->transdeps[i]->bounding_poly, pos, NULL);
    }
    if (prog->transdeps[i]->dest == stmt->id) {
      pluto_constraints_add_dim(prog->transdeps[i]->dpolytope,
                                prog->stmts[prog->transdeps[i]->src]->dim + pos,
                                NULL);
      pluto_constraints_add_dim(prog->transdeps[i]->bounding_poly,
                                prog->stmts[prog->transdeps[i]->src]->dim + pos,
                                NULL);
    }
  }
}

/* Warning: use it only to knock off a dummy dimension (unrelated to
 * anything else */
void pluto_stmt_remove_dim(Stmt *stmt, unsigned pos, PlutoProg *prog) {
  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(pos <= stmt->dim);
  assert(stmt->dim + npar + 1 == stmt->domain->ncols);

  pluto_constraints_remove_dim(stmt->domain, pos);
  stmt->dim--;

  if (stmt->iterators != NULL) {
    free(stmt->iterators[pos]);
    for (int i = pos; i <= (int)stmt->dim - 1; i++) {
      stmt->iterators[i] = stmt->iterators[i + 1];
    }
    stmt->iterators =
        (char **)realloc(stmt->iterators, stmt->dim * sizeof(char *));
  }

  pluto_matrix_remove_col(stmt->trans, pos);

  /* Update is_orig_loop */
  for (int i = pos; i <= (int)stmt->dim - 1; i++) {
    stmt->is_orig_loop[i] = stmt->is_orig_loop[i + 1];
  }
  stmt->is_orig_loop =
      (bool *)realloc(stmt->is_orig_loop, sizeof(bool) * stmt->dim);

  for (int i = 0; i < stmt->nwrites; i++) {
    pluto_matrix_remove_col(stmt->writes[i]->mat, pos);
  }

  for (int i = 0; i < stmt->nreads; i++) {
    pluto_matrix_remove_col(stmt->reads[i]->mat, pos);
  }

  /* Update deps */
  for (int i = 0; i < prog->ndeps; i++) {
    if (prog->deps[i]->src == stmt->id) {
      pluto_constraints_remove_dim(prog->deps[i]->dpolytope, pos);
    }
    if (prog->deps[i]->dest == stmt->id) {
      pluto_constraints_remove_dim(prog->deps[i]->dpolytope,
                                   prog->stmts[prog->deps[i]->src]->dim + pos);
    }
  }

  for (int i = 0; i < prog->ntransdeps; i++) {
    assert(prog->transdeps[i] != NULL);
    if (prog->transdeps[i]->src == stmt->id) {
      pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope, pos);
    }
    if (prog->transdeps[i]->dest == stmt->id) {
      pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope,
                                   prog->stmts[prog->transdeps[i]->src]->dim +
                                       pos);
    }
  }
}

void pluto_stmt_add_hyperplane(Stmt *stmt, PlutoHypType type, unsigned pos) {
  assert(pos <= stmt->trans->nrows);

  pluto_matrix_add_row(stmt->trans, pos);

  stmt->hyp_types = (PlutoHypType *)realloc(
      stmt->hyp_types, sizeof(PlutoHypType) * stmt->trans->nrows);
  for (int i = stmt->trans->nrows - 2; i >= (int)pos; i--) {
    stmt->hyp_types[i + 1] = stmt->hyp_types[i];
  }
  stmt->hyp_types[pos] = type;

  if (stmt->first_tile_dim >= (int)pos)
    stmt->first_tile_dim++;
  if (stmt->last_tile_dim >= (int)pos)
    stmt->last_tile_dim++;
}

void pluto_prog_add_hyperplane(PlutoProg *prog, int pos,
                               PlutoHypType hyp_type) {
  int i;

  prog->num_hyperplanes++;
  prog->hProps = (HyperplaneProperties *)realloc(
      prog->hProps, prog->num_hyperplanes * sizeof(HyperplaneProperties));

  for (i = prog->num_hyperplanes - 2; i >= pos; i--) {
    prog->hProps[i + 1] = prog->hProps[i];
  }
  /* Initialize some */
  prog->hProps[pos].unroll = NO_UNROLL;
  prog->hProps[pos].prevec = 0;
  prog->hProps[pos].band_num = -1;
  prog->hProps[pos].dep_prop = UNKNOWN;
  prog->hProps[pos].type = hyp_type;
}

/*
 * Create a new statement (see also pluto_stmt_dup)
 */
Stmt *pluto_create_stmt(int dim, const PlutoConstraints *domain,
                        const PlutoMatrix *trans, char **iterators,
                        const char *text, PlutoStmtType type) {
  Stmt *stmt = pluto_stmt_alloc(dim, domain, trans);

  stmt->type = type;

  stmt->text = strdup(text);

  for (unsigned i = 0; i < stmt->dim; i++) {
    stmt->iterators[i] = strdup(iterators[i]);
  }

  pluto_constraints_set_names_range(stmt->domain, stmt->iterators, 0, 0,
                                    stmt->dim);

  /* TODO: Set names for parameters */

  return stmt;
}

/* Pad statement transformations so that they all equal number
 * of rows */
void pluto_pad_stmt_transformations(PlutoProg *prog) {
  int i, nstmts;

  nstmts = prog->nstmts;
  Stmt **stmts = prog->stmts;

  /* Pad all trans if necessary with zeros */
  unsigned max_nrows = 0;
  for (i = 0; i < nstmts; i++) {
    if (stmts[i]->trans != NULL) {
      max_nrows = PLMAX(max_nrows, stmts[i]->trans->nrows);
    }
  }

  if (max_nrows >= 1) {
    for (i = 0; i < nstmts; i++) {
      if (stmts[i]->trans == NULL) {
        stmts[i]->trans =
            pluto_matrix_alloc(max_nrows, stmts[i]->dim + prog->npar + 1);
        stmts[i]->trans->nrows = 0;
      }

      unsigned curr_rows = stmts[i]->trans->nrows;

      /* Add all zero rows */
      for (unsigned j = curr_rows; j < max_nrows; j++) {
        pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
      }
    }

    unsigned old_hyp_num = prog->num_hyperplanes;
    for (unsigned i = old_hyp_num; i < max_nrows; i++) {
      /* This is not really H_SCALAR, but this is the best we can do */
      pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);
    }
  }
}

/* Add statement to program; can't reuse arg stmt pointer any more */
void pluto_add_given_stmt(PlutoProg *prog, Stmt *stmt) {
  prog->stmts =
      (Stmt **)realloc(prog->stmts, ((prog->nstmts + 1) * sizeof(Stmt *)));

  stmt->id = prog->nstmts;

  prog->nvar = PLMAX(prog->nvar, (int)stmt->dim);
  prog->stmts[prog->nstmts] = stmt;
  prog->nstmts++;

  pluto_pad_stmt_transformations(prog);
}

/* Create a statement and add it to the program
 * iterators: domain iterators
 * trans: schedule/transformation
 * domain: domain
 * text: statement text
 */
void pluto_add_stmt(PlutoProg *prog, const PlutoConstraints *domain,
                    const PlutoMatrix *trans, char **iterators,
                    const char *text, PlutoStmtType type) {
  int nstmts;

  assert(trans != NULL);
  assert(trans->ncols == domain->ncols);

  nstmts = prog->nstmts;

  prog->stmts = (Stmt **)realloc(prog->stmts, ((nstmts + 1) * sizeof(Stmt *)));

  Stmt *stmt = pluto_create_stmt(domain->ncols - prog->npar - 1, domain, trans,
                                 iterators, text, type);
  stmt->id = nstmts;

  /* Initialize intra statement deps to Null. Will be updated when fcg is built
   */
  stmt->intra_stmt_dep_cst = NULL;

  prog->nvar = PLMAX(prog->nvar, (int)stmt->dim);

  prog->stmts[nstmts] = stmt;
  prog->nstmts++;

  pluto_pad_stmt_transformations(prog);
}

Dep *pluto_dep_alloc() {
  Dep *dep = (Dep *)malloc(sizeof(Dep));

  dep->id = -1;
  dep->satvec = NULL;
  dep->dpolytope = NULL;
  dep->bounding_poly = NULL;
  dep->depsat_poly = NULL;
  dep->satisfied = false;
  dep->satisfaction_level = -1;
  dep->dirvec = NULL;
  dep->src_acc = NULL;
  dep->dest_acc = NULL;
  dep->cst = NULL;
  dep->bounding_cst = NULL;
  dep->src_unique_dpolytope = NULL;

  return dep;
}

Dep *pluto_dep_dup(Dep *d) {
  Dep *dep = (Dep *)malloc(sizeof(Dep));

  dep->id = d->id;
  dep->src = d->src;
  dep->dest = d->dest;
  dep->src_acc = d->src_acc;
  dep->dest_acc = d->dest_acc;
  dep->dpolytope = pluto_constraints_dup(d->dpolytope);
  dep->bounding_poly = pluto_constraints_dup(d->bounding_poly);

  dep->src_unique_dpolytope =
      d->src_unique_dpolytope ? pluto_constraints_dup(d->src_unique_dpolytope)
                              : NULL;

  dep->depsat_poly =
      d->depsat_poly ? pluto_constraints_dup(d->depsat_poly) : NULL;
  dep->satvec = NULL; // TODO
  dep->type = d->type;
  dep->satisfied = d->satisfied;
  dep->satisfaction_level = d->satisfaction_level;
  dep->dirvec = NULL; // TODO
  dep->cst = d->cst ? pluto_constraints_dup(d->cst) : NULL;
  dep->bounding_cst =
      d->bounding_cst ? pluto_constraints_dup(d->bounding_cst) : NULL;

  return dep;
}

/*
 * Only very essential information is needed to allocate; rest can be
 * populated as needed
 */
Stmt *pluto_stmt_alloc(unsigned dim, const PlutoConstraints *domain,
                       const PlutoMatrix *trans) {
  /* Have to provide a transformation */
  assert(trans != NULL);

  Stmt *stmt = (Stmt *)malloc(sizeof(Stmt));

  /* id will be assigned when added to PlutoProg */
  stmt->id = -1;
  stmt->dim = dim;
  stmt->dim_orig = dim;
  if (domain != NULL) {
    stmt->domain = pluto_constraints_dup(domain);
  } else {
    stmt->domain = NULL;
  }

  stmt->trans = pluto_matrix_dup(trans);

  stmt->hyp_types =
      (PlutoHypType *)malloc(stmt->trans->nrows * sizeof(PlutoHypType));
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    stmt->hyp_types[i] = H_LOOP;
  }

  stmt->text = NULL;
  stmt->tile = 1;
  stmt->num_tiled_loops = 0;
  stmt->reads = NULL;
  stmt->writes = NULL;
  stmt->nreads = 0;
  stmt->nwrites = 0;

  /* For diamond tiling */
  stmt->evicted_hyp = NULL;
  stmt->evicted_hyp_pos = -1;

  stmt->first_tile_dim = 0;
  stmt->last_tile_dim = -1;

  stmt->type = STMT_UNKNOWN;
  stmt->ploop_id = -1;

  if (dim >= 1) {
    stmt->is_orig_loop = (bool *)malloc(dim * sizeof(bool));
    stmt->iterators = (char **)malloc(sizeof(char *) * dim);
    for (unsigned i = 0; i < stmt->dim; i++) {
      stmt->iterators[i] = NULL;
    }
  } else {
    stmt->is_orig_loop = NULL;
    stmt->iterators = NULL;
  }

  return stmt;
}

PlutoAccess *pluto_access_dup(const PlutoAccess *acc) {
  assert(acc);

  PlutoAccess *nacc = (PlutoAccess *)malloc(sizeof(PlutoAccess));
  nacc->mat = pluto_matrix_dup(acc->mat);
  nacc->name = strdup(acc->name);
  nacc->sym_id = acc->sym_id;

  return nacc;
}

void pluto_access_free(PlutoAccess *acc) {
  if (acc) {
    pluto_matrix_free(acc->mat);
    free(acc->name);
    free(acc);
  }
}

void pluto_stmt_free(Stmt *stmt) {
  pluto_constraints_free(stmt->domain);

  pluto_matrix_free(stmt->trans);

  free(stmt->hyp_types);
  free(stmt->text);

  for (unsigned j = 0; j < stmt->dim; j++) {
    if (stmt->iterators[j] != NULL) {
      free(stmt->iterators[j]);
    }
  }

  /* If dim is zero, iterators, is_orig_loop are NULL */
  if (stmt->iterators != NULL) {
    free(stmt->iterators);
    free(stmt->is_orig_loop);
  }

  PlutoAccess **writes = stmt->writes;
  PlutoAccess **reads = stmt->reads;

  if (writes != NULL) {
    for (int i = 0; i < stmt->nwrites; i++) {
      pluto_access_free(writes[i]);
    }
    free(writes);
  }
  if (reads != NULL) {
    for (int i = 0; i < stmt->nreads; i++) {
      pluto_access_free(reads[i]);
    }
    free(reads);
  }

  pluto_matrix_free(stmt->evicted_hyp);

  free(stmt);
}

/* Get transformed domain */
PlutoConstraints *pluto_get_new_domain(const Stmt *stmt) {
  PlutoConstraints *sched;

  PlutoConstraints *newdom = pluto_constraints_dup(stmt->domain);
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    pluto_constraints_add_dim(newdom, 0, NULL);
  }

  sched = pluto_stmt_get_schedule(stmt);

  pluto_constraints_intersect(newdom, sched);

  pluto_constraints_project_out(newdom, stmt->trans->nrows, stmt->dim);

  pluto_constraints_free(sched);

  return newdom;
}

/*
 * Checks if the range of the variable at depth 'depth' can be bound by a
 * constant; returns the constant of -1 if it can't be
 *
 * WARNING: If cnst is a list, looks at just the first element
 *
 * TODO: Not general now: difference being constant can be implied through
 * other inequalities
 *
 * */
int get_const_bound_difference(const PlutoConstraints *cnst, int depth) {
  int constdiff, c, _lcm;

  assert(cnst != NULL);
  PlutoConstraints *cst = pluto_constraints_dup(cnst);

  pluto_constraints_project_out(cst, depth + 1, cst->ncols - 1 - depth - 1);
  assert(depth >= 0 && depth <= (int)cst->ncols - 2);

  constdiff = INT_MAX;

  unsigned r;
  for (r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0)
      break;
  }
  /* Variable doesn't appear */
  if (r == cst->nrows)
    return -1;

  /* Scale rows so that the coefficient of depth var is the same */
  _lcm = 1;
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0)
      _lcm = lcm(_lcm, llabs(cst->val[r][depth]));
  }
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0) {
      for (unsigned c = 0; c < cst->ncols; c++) {
        cst->val[r][c] = cst->val[r][c] * (_lcm / llabs(cst->val[r][depth]));
      }
    }
  }

  /* Equality to a function of parameters/constant implies single point */
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r] && cst->val[r][depth] != 0) {
      for (c = depth + 1; c < (int)cst->ncols - 1; c++) {
        if (cst->val[r][c] != 0) {
          break;
        }
      }
      if (c == (int)cst->ncols - 1) {
        constdiff = 1;
        // printf("constdiff is 1\n");
      }
    }
  }

  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r])
      continue;
    if (cst->val[r][depth] <= -1) {
      /* Find a lower bound with constant difference */
      for (unsigned r1 = 0; r1 < cst->nrows; r1++) {
        if (cst->is_eq[r1])
          continue;
        if (cst->val[r1][depth] >= 1) {
          for (c = 0; c < (int)cst->ncols - 1; c++) {
            if (cst->val[r1][c] + cst->val[r][c] != 0) {
              break;
            }
          }
          if (c == (int)cst->ncols - 1) {
            constdiff = PLMIN(
                constdiff,
                floorf(cst->val[r][c] / (float)-cst->val[r][depth]) +
                    ceilf(cst->val[r1][c] / (float)cst->val[r1][depth]) + 1);
          }
        }
      }
    }
  }
  pluto_constraints_free(cst);

  if (constdiff == INT_MAX) {
    return -1;
  }

  /* Sometimes empty sets imply negative difference */
  /* It basically means zero points */
  if (constdiff <= -1)
    constdiff = 0;

  return constdiff;
}

#define MINF 0
#define MAXF 1

/* Get expression for pos^{th} constraint in cst;
 * Returned string should be freed with 'free' */
char *get_expr(PlutoConstraints *cst, int pos, const char **params,
               int bound_type) {
  int c, sum;

  char *expr = (char *)malloc(512);
  strcpy(expr, "");

  // printf("Get expr\n");
  // pluto_constraints_print(stdout, cst);

  if (bound_type == MINF)
    assert(cst->val[pos][0] <= -1);
  else
    assert(cst->val[pos][0] >= 1);

  sum = 0;
  for (c = 1; c < (int)cst->ncols - 1; c++) {
    sum += llabs(cst->val[pos][c]);
  }

  if (sum == 0) {
    /* constant */
    if (bound_type == MINF) {
      sprintf(expr + strlen(expr), "%d",
              (int)floorf(cst->val[pos][cst->ncols - 1] /
                          -(float)cst->val[pos][0]));
    } else {
      sprintf(
          expr + strlen(expr), "%d",
          (int)ceilf(-cst->val[pos][cst->ncols - 1] / (float)cst->val[pos][0]));
    }
  } else {
    /* if it's being divided by 1, make it better by not putting
     * floor/ceil */
    if (llabs(cst->val[pos][0]) != 1) {
      if (bound_type == MINF) {
        sprintf(expr + strlen(expr), "floorf((");
      } else {
        sprintf(expr + strlen(expr), "ceilf((");
      }
    }

    for (c = 1; c < (int)cst->ncols - 1; c++) {
      if (cst->val[pos][c] != 0) {
        if (bound_type == MINF) {
          sprintf(expr + strlen(expr),
                  (cst->val[pos][c] >= 1) ? "+%ld*%s" : "%ld*%s",
                  cst->val[pos][c], params[c - 1]);
        } else {
          sprintf(expr + strlen(expr),
                  (cst->val[pos][c] <= -1) ? "+%ld*%s" : "%ld*%s",
                  -cst->val[pos][c], params[c - 1]);
        }
      }
    }

    if (cst->val[pos][c] != 0) {
      if (bound_type == MINF) {
        sprintf(expr + strlen(expr), (cst->val[pos][c] >= 1) ? "+%ld" : "%ld",
                cst->val[pos][c]);
      } else {
        sprintf(expr + strlen(expr), (cst->val[pos][c] <= -1) ? "+%ld" : "%ld",
                -cst->val[pos][c]);
      }
    }

    /* If it's being divided by 1, make it better by not putting
     * floor/ceil. */
    if (llabs(cst->val[pos][0]) != 1) {
      sprintf(expr + strlen(expr), ")/(float)%ld)",
              (bound_type == MINF) ? -cst->val[pos][0] : cst->val[pos][0]);
    }
  }

  return expr;
}

/*
 * Get min or max of all upper or lower bounds (resp).
 * Returned string should be freed with free
 */
char *get_func_of_expr(PlutoConstraints *cst, int offset, int bound_type,
                       const char **params) {
  char *expr, *expr1;

  char *fexpr = (char *)malloc(512);

  strcpy(fexpr, "");

  char func[5];
  if (bound_type == MINF)
    strcpy(func, "min(");
  else
    strcpy(func, "max(");

  if (cst->nrows - offset == 1) {
    expr = get_expr(cst, offset, params, bound_type);
    strcat(fexpr, expr);
  } else {
    /* cst->nrows >= 2 */
    expr = get_expr(cst, offset, params, bound_type);
    strcat(fexpr, func);
    strcat(fexpr, expr);
    expr1 = get_func_of_expr(cst, offset + 1, bound_type, params);
    strcat(fexpr, ",");
    strcat(fexpr, expr1);
    strcat(fexpr, ")");
    free(expr1);
  }
  free(expr);

  return fexpr;
}

/* Return the size of the parametric bounding box for a (contiguous)
 * block of dimensions
 * start: position of start of block
 * num: number of dimensions in block
 * npar: number of parameters in terms of which expression will be computed;
 * these are assumed to be the last 'npar' variables of cst
 * parmas: strings for 'npar' parameters
 * Return: expression describing the maximum number of points 'block'
 * vars traverse for any value of '0..start-1' variables
 *
 * This function is constant-aware, i.e., if possible, it will exploit the
 * fact that the range of a variable is bounded by a constant. The underlying
 * call to get_parametric_extent_const for each of the 'num' dimensions
 * achieves this.
 */
char *get_parametric_bounding_box(const PlutoConstraints *cst, int start,
                                  int num, int npar, const char **params) {
  int k;
  char *buf_size;

  buf_size = (char *)malloc(2048 * 8);
  strcpy(buf_size, "(");

  const PlutoConstraints *cst_tmp = cst;
  while (cst_tmp != NULL) {
    sprintf(buf_size + strlen(buf_size), "+1");
    for (k = 0; k < num; k++) {
      char *extent;
      get_parametric_extent_const(cst_tmp, start + k, npar, params, &extent,
                                  NULL);
      sprintf(buf_size + strlen(buf_size), "*(%s)", extent);
      free(extent);
    }
    cst_tmp = cst_tmp->next;
  }
  sprintf(buf_size + strlen(buf_size), ")");

  return buf_size;
}

/*  Parametric extent of the pos^th variable in cst
 *  Extent computation is constant-aware, i.e., look when it can be
 *  bounded by a constant; if not, just a difference of max and min
 *  expressions of parameters is returned;  last 'npar'  ones are
 *  treated as parameters; *extent should be freed by 'free'
 */
void get_parametric_extent_const(const PlutoConstraints *cst, int pos, int npar,
                                 const char **params, char **extent,
                                 char **p_lbexpr) {
  int constdiff = get_const_bound_difference(cst, pos);

  if ((p_lbexpr == NULL) && (constdiff != -1)) {
    *extent = (char *)malloc(sizeof(int) * 8);
    sprintf(*extent, "%d", constdiff);
  } else {
    get_parametric_extent(cst, pos, npar, params, extent, p_lbexpr);
  }
}

/* Get lower and upper bound expression as a function of parameters for pos^th
 * variable; last npar in cst are treated as parameters
 * lbexpr and ubexpr should be freed with free
 * */
void get_lb_ub_expr(const PlutoConstraints *cst, int pos, int npar,
                    const char **params, char **lbexpr, char **ubexpr) {
  PlutoConstraints *lb, *ub, *lbs, *ubs;
  char *lbe, *ube;

  PlutoConstraints *dup = pluto_constraints_dup(cst);

  pluto_constraints_project_out(dup, 0, pos);
  pluto_constraints_project_out(dup, 1, dup->ncols - npar - 1 - 1);

  lbs = pluto_constraints_alloc(dup->nrows, dup->ncols);
  ubs = pluto_constraints_alloc(dup->nrows, dup->ncols);

  for (unsigned i = 0; i < dup->nrows; i++) {
    if (dup->is_eq[i] && dup->val[i][0] != 0) {
      lb = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(lbs, lb);
      pluto_constraints_free(lb);

      ub = pluto_constraints_select_row(dup, i);
      pluto_constraints_negate_row(ub, 0);
      pluto_constraints_add(ubs, ub);
      pluto_constraints_free(ub);
    }
    if (dup->val[i][0] >= 1) {
      /* Lower bound */
      lb = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(lbs, lb);
      pluto_constraints_free(lb);
    } else if (dup->val[i][0] <= -1) {
      /* Upper bound */
      ub = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(ubs, ub);
      pluto_constraints_free(ub);
    }
  }

  assert(lbs->nrows >= 1);
  assert(ubs->nrows >= 1);
  pluto_constraints_free(dup);

  lbe = get_func_of_expr(lbs, 0, MAXF, params);
  ube = get_func_of_expr(ubs, 0, MINF, params);

  *lbexpr = lbe;
  *ubexpr = ube;

  pluto_constraints_free(lbs);
  pluto_constraints_free(ubs);
}

/*
 * Get expression for difference of upper and lower bound of pos^th variable
 * in cst in terms of parameters;  last 'npar' dimensions of cst are treated
 * as parameters; *extent should be freed by 'free'
 */
void get_parametric_extent(const PlutoConstraints *cst, int pos, int npar,
                           const char **params, char **extent,
                           char **p_lbexpr) {
  char *lbexpr, *ubexpr;

  get_lb_ub_expr(cst, pos, npar, params, &lbexpr, &ubexpr);

  if (!strcmp(lbexpr, ubexpr)) {
    *extent = strdup("1");
  } else {
    *extent =
        (char *)malloc(strlen(lbexpr) + strlen(ubexpr) + strlen(" -  + 1") + 1);
    sprintf(*extent, "%s - %s + 1", ubexpr, lbexpr);
  }
  if (p_lbexpr != NULL) {
    *p_lbexpr = (char *)malloc(strlen(lbexpr) + 1);
    strcpy(*p_lbexpr, lbexpr);
  }
  free(lbexpr);
  free(ubexpr);
}

/* Get Alpha matrix (A matrix - INRIA transformation representation */
PlutoMatrix *get_alpha(const Stmt *stmt, const PlutoProg *prog) {
  PlutoMatrix *a;
  a = pluto_matrix_alloc(stmt->dim, stmt->dim);

  unsigned r = 0;
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    if (stmt->hyp_types[i] == H_LOOP ||
        stmt->hyp_types[i] == H_TILE_SPACE_LOOP) {
      for (unsigned c = 0; c < stmt->dim; c++) {
        a->val[r][c] = stmt->trans->val[i][c];
      }
      r++;
      if (r == stmt->dim)
        break;
    }
  }

  assert(r == stmt->dim);

  return a;
}

bool pluto_is_hyperplane_scalar(const Stmt *stmt, int level) {
  assert(level <= (int)stmt->trans->nrows - 1);

  for (unsigned j = 0; j < stmt->dim; j++) {
    if (stmt->trans->val[level][j] != 0)
      return false;
  }

  return true;
}

int pluto_is_hyperplane_loop(const Stmt *stmt, int level) {
  return !pluto_is_hyperplane_scalar(stmt, level);
}

/* Get the remapping matrix: maps time iterators back to the domain
 * iterators; divs: divisors for the rows */
PlutoMatrix *pluto_stmt_get_remapping(const Stmt *stmt, int **divs) {
  PlutoMatrix *trans = stmt->trans;
  PlutoMatrix *remap = pluto_matrix_dup(trans);

  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  *divs = (int *)malloc(sizeof(int) * (stmt->dim + npar + 1));

  for (unsigned i = 0; i < remap->nrows; i++) {
    pluto_matrix_negate_row(remap, remap->nrows - 1 - i);
    pluto_matrix_add_col(remap, 0);
    remap->val[trans->nrows - 1 - i][0] = 1;
  }

  /* Bring the stmt iterators to the left */
  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_move_col(remap, remap->nrows + i, i);
  }

  assert(stmt->dim <= remap->nrows);

  for (unsigned i = 0; i < stmt->dim; i++) {
    if (remap->val[i][i] == 0) {
      unsigned k;
      for (k = i + 1; k < remap->nrows; k++) {
        if (remap->val[k][i] != 0)
          break;
      }
      if (k < remap->nrows) {
        pluto_matrix_interchange_rows(remap, i, k);
      } else {
        /* Can't associate domain iterator with time iterator */
        /* Shouldn't happen with a full-ranked transformation */
        printf("Can't associate domain iterator #%d with time iterators\n",
               i + 1);
        pluto_matrix_print(stdout, remap);
        assert(0);
      }
    }
    assert(remap->val[i][i] != 0);
    for (unsigned k = i + 1; k < remap->nrows; k++) {
      if (remap->val[k][i] == 0)
        continue;
      int64_t _lcm = lcm(remap->val[k][i], remap->val[i][i]);
      int64_t factor1 = _lcm / remap->val[k][i];
      for (unsigned j = i; j < remap->ncols; j++) {
        remap->val[k][j] = remap->val[k][j] * factor1 -
                           remap->val[i][j] * (_lcm / remap->val[i][i]);
      }
    }
  }

  /* Solve upper triangular system now */
  for (int i = stmt->dim - 1; i >= 0; i--) {
    assert(remap->val[i][i] != 0);
    for (int k = i - 1; k >= 0; k--) {
      if (remap->val[k][i] == 0)
        continue;
      int64_t _lcm = lcm(remap->val[k][i], remap->val[i][i]);
      int64_t factor1 = _lcm / remap->val[k][i];
      for (unsigned j = 0; j < remap->ncols; j++) {
        remap->val[k][j] = remap->val[k][j] * (factor1)-remap->val[i][j] *
                           (_lcm / remap->val[i][i]);
      }
    }
  }

  assert(remap->nrows >= stmt->dim);
  for (int i = remap->nrows - 1; i >= (int)stmt->dim; i--) {
    pluto_matrix_remove_row(remap, remap->nrows - 1);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    assert(remap->val[i][i] != 0);
    if (remap->val[i][i] <= -1) {
      pluto_matrix_negate_row(remap, i);
    }
    (*divs)[i] = llabs(remap->val[i][i]);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_remove_col(remap, 0);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_negate_row(remap, i);
  }

  /* Identity for the parameter and constant part */
  for (int i = 0; i < npar + 1; i++) {
    pluto_matrix_add_row(remap, remap->nrows);
    remap->val[remap->nrows - 1][remap->ncols - npar - 1 + i] = 1;
    (*divs)[stmt->dim + i] = 1;
  }

  assert(remap->nrows == stmt->dim + npar + 1 && "Invalid remapping matrix");

  return remap;
}

void pluto_prog_params_print(const PlutoProg *prog) {
  for (int i = 0; i < prog->npar; i++) {
    printf("%s\n", prog->params[i]);
  }
}

/// Constructs the new access function for `acc' resulting from the
/// transformation attached to 'stmt'.
PlutoMatrix *pluto_get_new_access_func(const PlutoMatrix *acc, const Stmt *stmt,
                                       int **divs) {
  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(acc->ncols == stmt->dim + npar + 1 &&
         "Invalid access function matrix");

  *divs = (int *)malloc(sizeof(int) * acc->nrows);

  int *remap_divs;
  PlutoMatrix *remap = pluto_stmt_get_remapping(stmt, &remap_divs);
  assert(remap->nrows == stmt->dim + npar + 1);

  int _lcm = 1;
  for (unsigned r = 0; r < remap->nrows; r++) {
    assert(remap_divs[r] != 0);
    _lcm = lcm(_lcm, remap_divs[r]);
  }
  for (unsigned r = 0; r < remap->nrows; r++) {
    for (unsigned c = 0; c < remap->ncols; c++) {
      remap->val[r][c] = (remap->val[r][c] * _lcm) / remap_divs[r];
    }
  }

  PlutoMatrix *newacc = pluto_matrix_product(acc, remap);

  for (unsigned r = 0; r < newacc->nrows; r++) {
    (*divs)[r] = _lcm;
  }

  assert(newacc->ncols == stmt->trans->nrows + npar + 1);

  pluto_matrix_free(remap);
  free(remap_divs);

  return newacc;
}

/* Separates a list of statements at level 'level' */
void pluto_separate_stmts(PlutoProg *prog, Stmt **stmts, int num, int level,
                          int offset) {
  int i, nstmts, k;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
  }
  for (k = 0; k < num; k++) {
    stmts[k]->trans->val[level][stmts[k]->trans->ncols - 1] = offset + 1 + k;
  }

  pluto_prog_add_hyperplane(prog, level, H_SCALAR);
  prog->hProps[level].dep_prop = SEQ;
}

/* Separates a statement from the rest (places it later) at level 'level';
 * this is done by inserting a scalar dimension separating them */
void pluto_separate_stmt(PlutoProg *prog, const Stmt *stmt, int level) {
  int i, nstmts;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
  }
  stmt->trans->val[level][stmt->trans->ncols - 1] = 1;

  pluto_prog_add_hyperplane(prog, level, H_SCALAR);
  prog->hProps[level].dep_prop = SEQ;
}

int pluto_stmt_is_member_of(int stmt_id, Stmt **slist, int len) {
  int i;
  for (i = 0; i < len; i++) {
    if (stmt_id == slist[i]->id)
      return 1;
  }
  return 0;
}

int pluto_stmt_is_subset_of(Stmt **s1, int n1, Stmt **s2, int n2) {
  int i;

  for (i = 0; i < n1; i++) {
    if (!pluto_stmt_is_member_of(s1[i]->id, s2, n2))
      return 0;
  }

  return 1;
}

/* Add new to accs if it's an access to a variable not already contained in
 * accs */
void add_if_new_var(PlutoAccess ***accs, int *num, PlutoAccess *new_acc) {
  int i;

  for (i = 0; i < *num; i++) {
    if (!strcmp((*accs)[i]->name, new_acc->name)) {
      break;
    }
  }

  if (i == *num) {
    *accs = (PlutoAccess **)realloc(*accs, (*num + 1) * sizeof(PlutoAccess *));
    (*accs)[*num] = new_acc;
    (*num)++;
  }
}

/* Get all write accesses in the program */
PlutoAccess **pluto_get_all_waccs(const PlutoProg *prog, int *num) {
  int i;

  PlutoAccess **accs = NULL;
  *num = 0;

  for (i = 0; i < prog->nstmts; i++) {
    assert(prog->stmts[i]->nwrites == 1);
    add_if_new_var(&accs, num, prog->stmts[i]->writes[0]);
  }
  return accs;
}

int pluto_get_max_ind_hyps_non_scalar(const PlutoProg *prog) {
  int max, i;

  max = 0;

  for (i = 0; i < prog->nstmts; i++) {
    max = PLMAX(max, pluto_stmt_get_num_ind_hyps_non_scalar(prog->stmts[i]));
  }

  return max;
}

/*
 * The maximum number of linearly independent hyperplanes across all
 * statements
 */
int pluto_get_max_ind_hyps(const PlutoProg *prog) {
  unsigned max = 0;

  for (int i = 0; i < prog->nstmts; i++) {
    max = PLMAX(max, pluto_stmt_get_num_ind_hyps(prog->stmts[i]));
  }

  return max;
}

int pluto_stmt_get_num_ind_hyps_non_scalar(const Stmt *stmt) {
  PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

  /* Ignore padding dimensions, params, and constant part */
  for (unsigned i = stmt->dim_orig; i < stmt->trans->ncols; i++) {
    pluto_matrix_remove_col(tprime, stmt->dim_orig);
  }
  unsigned j = 0;
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    if (stmt->hyp_types[i] == H_SCALAR) {
      pluto_matrix_remove_row(tprime, i - j);
      j++;
    }
  }

  unsigned isols = pluto_matrix_get_rank(tprime);
  pluto_matrix_free(tprime);

  return isols;
}

unsigned pluto_stmt_get_num_ind_hyps(const Stmt *stmt) {
  PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

  /* Ignore padding dimensions, params, and constant part */
  for (unsigned i = stmt->dim_orig; i < stmt->trans->ncols; i++) {
    pluto_matrix_remove_col(tprime, stmt->dim_orig);
  }

  unsigned isols = pluto_matrix_get_rank(tprime);
  pluto_matrix_free(tprime);

  return isols;
}

/*
 * Are all transformations full column-ranked?
 */
int pluto_transformations_full_ranked(PlutoProg *prog) {
  for (int i = 0; i < prog->nstmts; i++) {
    if (pluto_stmt_get_num_ind_hyps(prog->stmts[i]) <
        prog->stmts[i]->dim_orig) {
      return 0;
    }
  }

  return 1;
}

struct acc_info {
  char *prefix;
  int acc_num;
  isl_union_map **new_maps;
  isl_union_map **schedule;
  isl_map *base_schedule;
};

osl_relation_p get_identity_schedule(int dim, int npar) {

  int i, j;
  osl_relation_p rln = osl_relation_pmalloc(PLUTO_OSL_PRECISION, 2 * dim + 1,
                                            dim + npar + 1 + 1);

  // copy matrix values
  for (i = 0; i < rln->nb_rows; i++) {
    osl_int_set_si(rln->precision, &rln->m[i][0], 0);
    for (j = 0; j < rln->nb_columns; j++) {
      osl_int_set_si(rln->precision, &rln->m[i][j], 0);
    }
  }

  for (i = 1; i < dim; i++) {
    osl_int_set_si(rln->precision, &rln->m[2 * i - 1][i], 1);
  }

  rln->type = OSL_TYPE_SCATTERING;
  rln->nb_parameters = npar;
  rln->nb_output_dims = dim;
  rln->nb_input_dims = 0;
  rln->nb_local_dims = 0;

  return rln;
}

/*
 * Return clone of a statement
 */
Stmt *pluto_stmt_dup(const Stmt *stmt) {
  Stmt *nstmt = pluto_stmt_alloc(stmt->dim, stmt->domain, stmt->trans);

  nstmt->dim_orig = stmt->dim_orig;
  nstmt->type = stmt->type;

  for (unsigned i = 0; i < stmt->dim; i++) {
    nstmt->iterators[i] = strdup(stmt->iterators[i]);
    nstmt->is_orig_loop[i] = stmt->is_orig_loop[i];
  }
  if (stmt->text)
    nstmt->text = strdup(stmt->text);

  nstmt->nreads = stmt->nreads;
  nstmt->nwrites = stmt->nwrites;

  nstmt->reads = (PlutoAccess **)malloc(nstmt->nreads * sizeof(PlutoAccess *));
  nstmt->writes =
      (PlutoAccess **)malloc(nstmt->nwrites * sizeof(PlutoAccess *));

  for (int i = 0; i < stmt->nreads; i++) {
    nstmt->reads[i] = pluto_access_dup(stmt->reads[i]);
  }

  for (int i = 0; i < stmt->nwrites; i++) {
    nstmt->writes[i] = pluto_access_dup(stmt->writes[i]);
  }

  return nstmt;
}

static void decrement_stmt_id(PlutoProg *prog, int id) {
  int i;

  prog->stmts[id]->id--;

  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (dep->src == id) {
      dep->src--;
    }
    if (dep->dest == id) {
      dep->dest--;
    }
  }
}

/* Add statement to program; can't reuse arg stmt pointer any more */
void pluto_remove_stmt(PlutoProg *prog, int stmt_id) {
  int i;

  pluto_stmt_free(prog->stmts[stmt_id]);

  for (i = stmt_id; i < prog->nstmts - 1; i++) {
    prog->stmts[i] = prog->stmts[i + 1];
    decrement_stmt_id(prog, prog->stmts[i]->id);
  }

  prog->nstmts--;

  prog->stmts =
      (Stmt **)realloc(prog->stmts, ((prog->nstmts) * sizeof(Stmt *)));

  for (i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, (int)prog->stmts[i]->dim);
  }
}

void pluto_transformations_pretty_print(const PlutoProg *prog) {
  int nstmts, i;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_transformation_print(prog->stmts[i]);
  }
}

void pluto_transformation_print_level(const PlutoProg *prog, int level) {
  int nstmts, i;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    fprintf(stdout, "h(S%d) = ", i + 1);
    pluto_stmt_print_hyperplane(stdout, prog->stmts[i], level);
    if (i < nstmts - 1)
      fprintf(stdout, ", ");
  }
  printf("\n");
}

/* List properties of newly found hyperplanes */
void pluto_print_hyperplane_properties(const PlutoProg *prog) {
  int j, numH;
  HyperplaneProperties *hProps;

  hProps = prog->hProps;
  numH = prog->num_hyperplanes;

  if (numH == 0) {
    fprintf(stdout, "No hyperplanes\n");
  }

  /* Note that loop properties are calculated for each dimension in the
   * transformed space (common for all statements) */
  for (j = 0; j < numH; j++) {
    fprintf(stdout, "t%d --> ", j + 1);
    switch (hProps[j].dep_prop) {
    case PARALLEL:
      fprintf(stdout, "parallel ");
      break;
    case SEQ:
      fprintf(stdout, "serial   ");
      break;
    case PIPE_PARALLEL:
      fprintf(stdout, "fwd_dep  ");
      break;
    default:
      fprintf(stdout, "unknown  ");
      break;
    }
    switch (hProps[j].type) {
    case H_LOOP:
      fprintf(stdout, "loop  ");
      break;
    case H_SCALAR:
      fprintf(stdout, "scalar");
      break;
    case H_TILE_SPACE_LOOP:
      fprintf(stdout, "tLoop ");
      break;
    default:
      fprintf(stdout, "unknown  ");
      break;
    }
    fprintf(stdout, " (band %d)", hProps[j].band_num);
    fprintf(stdout, hProps[j].unroll ? "ujam" : "no-ujam");
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
}

void pluto_transformations_print(const PlutoProg *prog) {
  int i;

  for (i = 0; i < prog->nstmts; i++) {
    printf("T_(S%d) \n", prog->stmts[i]->id + 1);
    pluto_matrix_print(stdout, prog->stmts[i]->trans);
  }
}

void pluto_stmt_transformation_print(const Stmt *stmt) {
  fprintf(stdout, "T(S%d): ", stmt->id + 1);
  printf("(");
  for (unsigned level = 0; level < stmt->trans->nrows; level++) {
    pluto_stmt_print_hyperplane(stdout, stmt, level);
    if ((int)level <= (int)stmt->trans->nrows - 2)
      printf(", ");
  }
  printf(")\n");

  printf("loop types (");
  for (unsigned level = 0; level < stmt->trans->nrows; level++) {
    if (level > 0)
      printf(", ");
    if (stmt->hyp_types[level] == H_SCALAR)
      printf("scalar");
    else if (stmt->hyp_types[level] == H_LOOP)
      printf("loop");
    else if (stmt->hyp_types[level] == H_TILE_SPACE_LOOP)
      printf("tloop");
    else
      printf("unknown");
  }
  printf(")\n\n");
}

/*
 * Computes the dependence polyhedron between the source iterators of dep1 and
 * dep2
 * domain1:  source iterators of dep1
 * domain2:  source iterators of dep2
 * dep1: first dependence
 * dep2: second dependence
 * access_matrix: access function matrix for source iterators of dep2. pass NULL
 * to use access function in dep2
 * returns dependence polyhedron
 */
PlutoConstraints *pluto_find_dependence(PlutoConstraints *domain1,
                                        PlutoConstraints *domain2, Dep *dep1,
                                        Dep *dep2, PlutoProg *prog,
                                        PlutoMatrix *access_matrix) {
  isl_ctx *ctx = isl_ctx_alloc();
  assert(ctx);

  isl_space *dim = isl_space_set_alloc(ctx, prog->npar, 0);
  dim = set_names(dim, isl_dim_param, prog->params);
  isl_space *param_space = isl_space_params(isl_space_copy(dim));
  isl_set *context = osl_relation_list_to_isl_set(
      pluto_constraints_to_osl_domain(prog->context, prog->npar), param_space);

  isl_union_map *empty = isl_union_map_empty(isl_space_copy(dim));
  isl_union_map *write = isl_union_map_empty(isl_space_copy(dim));
  isl_union_map *read = isl_union_map_empty(isl_space_copy(dim));
  isl_union_map *schedule = isl_union_map_empty(dim);

  // Add the source iterators of dep1 and corresponding access function to isl
  PlutoConstraints *source_iterators = domain2;
  PlutoAccess *access = dep2->src_acc;
  Stmt *s = prog->stmts[dep2->src];
  int domain_dim = source_iterators->ncols - prog->npar - 1;
  char **iter = (char **)malloc(domain_dim * sizeof(char *));

  for (int i = 0; i < domain_dim; i++) {
    iter[i] = (char *)malloc(12 * sizeof(char));
    sprintf(iter[i], "d%d", i + 1);
  }

  // assert(domain_dim <= s->dim);
  isl_map *read_pos;
  isl_map *write_pos;
  isl_map *schedule_i;

  char name[20];

  snprintf(name, sizeof(name), "S_%d_r%d", 0, 0);

  dim = isl_space_set_alloc(ctx, prog->npar, domain_dim);
  dim = set_names(dim, isl_dim_param, prog->params);
  dim = set_names(dim, isl_dim_set, iter);
  dim = isl_space_set_tuple_name(dim, isl_dim_set, name);

  isl_set *dom = osl_relation_list_to_isl_set(
      pluto_constraints_list_to_osl_domain(source_iterators, prog->npar), dim);

  dom = isl_set_intersect_params(dom, isl_set_copy(context));

  dim = isl_space_alloc(ctx, prog->npar, domain_dim, 2 * domain_dim + 1);
  dim = set_names(dim, isl_dim_param, prog->params);
  dim = set_names(dim, isl_dim_in, iter);
  dim = isl_space_set_tuple_name(dim, isl_dim_in, name);

  PlutoMatrix *i_schedule = get_identity_schedule_new(domain_dim, prog->npar);
  schedule_i = pluto_matrix_schedule_to_isl_map(i_schedule, dim);

  int *divs = NULL;
  if (access_matrix == NULL)
    read_pos = pluto_basic_access_to_isl_union_map(
        pluto_get_new_access_func(access->mat, s, &divs), access->name, dom);
  else
    read_pos =
        pluto_basic_access_to_isl_union_map(access_matrix, access->name, dom);
  read = isl_union_map_union(read, isl_union_map_from_map(read_pos));
  free(divs);

  schedule = isl_union_map_union(schedule, isl_union_map_from_map(schedule_i));

  for (int i = 0; i < domain_dim; i++) {
    free(iter[i]);
  }

  free(iter);

  // Add the source iterators of dep2 and corresponding access function to isl
  source_iterators = domain1;
  access = dep1->src_acc;
  s = prog->stmts[dep1->src];
  domain_dim = source_iterators->ncols - prog->npar - 1;

  iter = (char **)malloc(domain_dim * sizeof(char *));

  for (int i = 0; i < domain_dim; i++) {
    iter[i] = (char *)malloc(12 * sizeof(char));
    sprintf(iter[i], "d%d", i + 1);
  }

  snprintf(name, sizeof(name), "S_%d_w%d", 0, 0);

  dim = isl_space_set_alloc(ctx, prog->npar, domain_dim);
  dim = set_names(dim, isl_dim_param, prog->params);
  dim = set_names(dim, isl_dim_set, iter);
  dim = isl_space_set_tuple_name(dim, isl_dim_set, name);

  dom = osl_relation_list_to_isl_set(
      pluto_constraints_list_to_osl_domain(source_iterators, prog->npar), dim);

  dom = isl_set_intersect_params(dom, isl_set_copy(context));

  dim = isl_space_alloc(ctx, prog->npar, domain_dim, 2 * domain_dim + 1);
  dim = set_names(dim, isl_dim_param, prog->params);
  dim = set_names(dim, isl_dim_in, iter);
  dim = isl_space_set_tuple_name(dim, isl_dim_in, name);

  // osl_relation_free(smat);
  i_schedule = get_identity_schedule_new(domain_dim, prog->npar);
  schedule_i = pluto_matrix_schedule_to_isl_map(i_schedule, dim);

  write_pos = pluto_basic_access_to_isl_union_map(
      pluto_get_new_access_func(access->mat, s, &divs), access->name, dom);
  write = isl_union_map_union(write, isl_union_map_from_map(write_pos));

  schedule = isl_union_map_union(schedule, isl_union_map_from_map(schedule_i));

  isl_union_map *dep_raw;
  isl_union_map_compute_flow(
      isl_union_map_copy(read), isl_union_map_copy(empty),
      isl_union_map_copy(write), isl_union_map_copy(schedule), NULL, &dep_raw,
      NULL, NULL);

  dep_raw = isl_union_map_coalesce(dep_raw);

  int ndeps = 0;
  isl_union_map_foreach_map(dep_raw, &isl_map_count, &ndeps);

  if (ndeps == 0) {
    return NULL;
  }

  Dep **deps = (Dep **)malloc(ndeps * sizeof(Dep *));
  for (int i = 0; i < ndeps; i++) {
    deps[i] = pluto_dep_alloc();
  }
  ndeps = 0;
  ndeps += extract_deps(deps, ndeps, prog->stmts, dep_raw, OSL_DEPENDENCE_RAW);

  PlutoConstraints *tdpoly = NULL;
  for (int i = 0; i < ndeps; i++) {
    if (tdpoly == NULL)
      tdpoly = pluto_constraints_dup(deps[i]->dpolytope);
    else
      pluto_constraints_unionize(tdpoly, deps[i]->dpolytope);
  }

  isl_union_map_free(dep_raw);

  isl_union_map_free(empty);
  isl_union_map_free(write);
  isl_union_map_free(read);
  isl_union_map_free(schedule);
  isl_set_free(context);

  isl_ctx_free(ctx);

  for (int i = 0; i < domain_dim; i++) {
    free(iter[i]);
  }

  free(iter);

  return tdpoly;
}
