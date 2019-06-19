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
#include "transforms.h"
#include "constraints.h"
#include "pluto.h"
#include "program.h"

#include "assert.h"

/// Sink statement (domain); depth: 0-indexed.
void pluto_sink_statement(Stmt *stmt, int depth, int val, PlutoProg *prog) {
  assert(stmt->dim == stmt->domain->ncols - prog->npar - 1);

  char iter[13];
  snprintf(iter, sizeof(iter), "d%d", stmt->dim);

  pluto_stmt_add_dim(stmt, depth, -1, iter, H_SCALAR, prog);

  pluto_constraints_set_var(stmt->domain, depth, val);
  stmt->is_orig_loop[depth] = false;
}

/// Stripmine 'dim'th time dimension of stmt by stripmine factor; use
/// 'supernode' as the name of the supernode in the domain.
void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode,
                     PlutoProg *prog) {
  pluto_stmt_add_dim(stmt, 0, dim, supernode, H_TILE_SPACE_LOOP, prog);

  PlutoConstraints *domain = stmt->domain;

  pluto_constraints_add_inequality(domain);
  domain->val[domain->nrows - 1][0] = -factor;
  assert(stmt->trans->ncols == domain->ncols);
  for (int i = 1; i < (int)stmt->trans->ncols - 1; i++) {
    domain->val[domain->nrows - 1][i] = stmt->trans->val[dim + 1][i];
  }

  pluto_constraints_add_inequality(domain);
  domain->val[domain->nrows - 1][0] = factor;
  assert(stmt->trans->ncols == domain->ncols);
  for (int i = 1; i < (int)stmt->trans->ncols - 1; i++) {
    domain->val[domain->nrows - 1][i] = -stmt->trans->val[dim + 1][i];
  }
  domain->val[domain->nrows - 1][stmt->trans->ncols] += factor;
}

/// Interchange loops for a stmt.
void pluto_stmt_loop_interchange(Stmt *stmt, int level1, int level2) {
  for (unsigned j = 0; j < stmt->trans->ncols; j++) {
    int64_t tmp = stmt->trans->val[level1][j];
    stmt->trans->val[level1][j] = stmt->trans->val[level2][j];
    stmt->trans->val[level2][j] = tmp;
  }
}

void pluto_interchange(PlutoProg *prog, int level1, int level2) {
  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  for (int k = 0; k < nstmts; k++) {
    pluto_stmt_loop_interchange(stmts[k], level1, level2);
  }

  HyperplaneProperties hTmp = prog->hProps[level1];
  prog->hProps[level1] = prog->hProps[level2];
  prog->hProps[level2] = hTmp;
}

void pluto_sink_transformation(Stmt *stmt, unsigned pos) {
  int npar = stmt->domain->ncols - stmt->dim - 1;

  assert(pos <= stmt->trans->nrows);
  assert(stmt->dim + npar + 1 == stmt->domain->ncols);

  /* Stmt should always have a transformation */
  assert(stmt->trans != NULL);

  pluto_matrix_add_row(stmt->trans, pos);

  stmt->hyp_types = (PlutoHypType *)realloc(
      stmt->hyp_types, sizeof(PlutoHypType) * stmt->trans->nrows);
  for (int i = stmt->trans->nrows - 2; i >= (int)pos; i--) {
    stmt->hyp_types[i + 1] = stmt->hyp_types[i];
  }
  stmt->hyp_types[pos] = H_SCALAR;
}

/// Make loop the innermost loop; all loops below move up by one. If the loop
/// nest is tiled, then the loop can be  moved acorss scalar hyperplanes
/// in the innermost permutable band by setting the boolean variable
/// move_across_scalar_hyperplanes in the caller. The maximum depth to which a
/// loop can be moved is given by last_level, which is set to the last level in
/// the band containing the loop, in case of a tiled loop nest.
void pluto_make_innermost_loop(Ploop *loop, unsigned last_level,
                               bool move_across_scalar_hyperplanes,
                               PlutoProg *prog) {
  assert(prog->num_hyperplanes >= 1);
  assert(last_level <= prog->num_hyperplanes);

  if (!move_across_scalar_hyperplanes) {
    unsigned last_depth = prog->num_hyperplanes - 1;
    for (unsigned i = 0; i < loop->nstmts; i++) {
      Stmt *stmt = loop->stmts[i];
      unsigned d;
      for (d = loop->depth; d < stmt->trans->nrows; d++) {
        if (pluto_is_hyperplane_scalar(stmt, d)) {
          break;
        }
      }
      last_depth = PLMIN(last_depth, d - 1);
    }

    for (unsigned i = 0; i < loop->nstmts; i++) {
      Stmt *stmt = loop->stmts[i];
      for (unsigned d = loop->depth; d < last_depth; d++) {
        pluto_stmt_loop_interchange(stmt, d, d + 1);
      }
    }
    return;
  }

  unsigned last_depth = last_level;
  for (unsigned i = 0; i < loop->nstmts; i++) {
    Stmt *stmt = loop->stmts[i];
    /* Current level that has to be made the innermost. */
    unsigned current_level = loop->depth;
    for (unsigned d = loop->depth + 1; d < last_depth; d++) {
      if (pluto_is_hyperplane_scalar(stmt, d))
        continue;
      pluto_stmt_loop_interchange(stmt, current_level, d);
      current_level = d;
    }
  }
}
