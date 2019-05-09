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
#include <stdlib.h>
#include <stdio.h>

#include "pluto.h"
#include "program.h"
#include "ast_transform.h"

#include "cloog/cloog.h"

/*
 * Clast-based parallel loop marking */
void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
                         CloogOptions *cloogOptions) {
  unsigned i, j, nloops, nstmts, nploops;
  struct clast_for **loops;
  int *stmts;
  assert(root != NULL);

  Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);

  IF_DEBUG(printf("[pluto_mark_parallel] parallel loops\n"););
  IF_DEBUG(pluto_loops_print(ploops, nploops););

  for (i = 0; i < nploops; i++) {
    char iter[13];
    sprintf(iter, "t%d", ploops[i]->depth + 1);
    int *stmtids = (int *)malloc(ploops[i]->nstmts * sizeof(int));
    int max_depth = 0;
    for (j = 0; j < ploops[i]->nstmts; j++) {
      Stmt *stmt = ploops[i]->stmts[j];
      if (stmt->trans->nrows > max_depth)
        max_depth = stmt->trans->nrows;
      stmtids[j] = stmt->id + 1;
    }

    IF_DEBUG(printf("\tLooking for loop\n"););
    IF_DEBUG(pluto_loop_print(ploops[i]););
    // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

    ClastFilter filter = { iter, stmtids, (int)ploops[i]->nstmts, subset };
    clast_filter(root, filter, &loops, (int *)&nloops, &stmts, (int *)&nstmts);

    /* There should be at least one */
    if (nloops == 0) {
      /* Sometimes loops may disappear (1) tile size larger than trip count
       * 2) it's a scalar dimension but can't be determined from the
       * trans matrix */
      printf("Warning: parallel poly loop not found in AST\n");
      continue;
    } else {
      for (j = 0; j < nloops; j++) {
        loops[j]->parallel = CLAST_PARALLEL_NOT;
        char *private_vars = (char *)malloc(512);
        strcpy(private_vars, "lbv,ubv");
        if (options->parallel) {
          IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
          loops[j]->parallel = CLAST_PARALLEL_OMP;
          int depth = ploops[i]->depth + 1;
          for (depth++; depth <= max_depth; depth++) {
            sprintf(private_vars + strlen(private_vars), ",t%d", depth);
          }
        }
        loops[j]->private_vars = strdup(private_vars);
        free(private_vars);
      }
    }
    free(stmtids);
    free(loops);
    free(stmts);
  }

  pluto_loops_free(ploops, nploops);
}

/*
 * Clast-based vector loop marking */
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
                       CloogOptions *cloogOptions) {
  unsigned i, j, nloops, nstmts, nploops;
  struct clast_for **loops;
  int *stmts;
  assert(root != NULL);

  Ploop **ploops = pluto_get_parallel_loops(prog, &nploops);

  IF_DEBUG(printf("[pluto_mark_vector] parallel loops\n"););
  IF_DEBUG(pluto_loops_print(ploops, nploops););

  for (i = 0; i < nploops; i++) {
    /* Only the innermost ones */
    if (!pluto_loop_is_innermost(ploops[i], prog))
      continue;

    IF_DEBUG(printf("[pluto_mark_vector] marking loop vectorizable\n"););
    IF_DEBUG(pluto_loop_print(ploops[i]););
    char iter[13];
    sprintf(iter, "t%d", ploops[i]->depth + 1);
    int *stmtids = (int *)malloc(ploops[i]->nstmts * sizeof(int));
    for (j = 0; j < ploops[i]->nstmts; j++) {
      stmtids[j] = ploops[i]->stmts[j]->id + 1;
    }

    ClastFilter filter = { iter, stmtids, (int)ploops[i]->nstmts, subset };
    clast_filter(root, filter, &loops, (int *)&nloops, &stmts, (int *)&nstmts);

    /* There should be at least one */
    if (nloops == 0) {
      /* Sometimes loops may disappear (1) tile size larger than trip count
       * 2) it's a scalar dimension but can't be determined from the
       * trans matrix */
      printf("[pluto] pluto_mark_vector: WARNING: vectorizable poly loop not "
             "found in AST\n");
      continue;
    }
    for (j = 0; j < nloops; j++) {
      loops[j]->parallel += CLAST_PARALLEL_VEC;
    }
    free(stmtids);
    free(loops);
    free(stmts);
  }

  pluto_loops_free(ploops, nploops);
}
