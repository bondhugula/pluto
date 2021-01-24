/*
 * Pluto: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE in the
 * top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "ast_transform.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "program.h"

#include "cloog/cloog.h"

/// Clast-based parallel loop marking.
void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
                         CloogOptions *cloogOptions) {
  assert(root != NULL);
  unsigned nploops;
  Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);
  Band **pbands = (Band **)malloc(nploops * sizeof(Band *));
  PlutoContext *context = prog->context;
  PlutoOptions *options = prog->context->options;

  IF_DEBUG(printf("[pluto_mark_parallel] parallel loops\n"););
  IF_DEBUG(pluto_loops_print(ploops, nploops););

  for (unsigned i = 0; i < nploops; i++) {
    int innermost_split_level;
    pbands[i] =
        pluto_get_parallel_band(ploops[i], prog, &innermost_split_level);
    assert(pbands[i]->width > 0);
  }
  unsigned npbands = nploops;

  /* Iterate over bands and for each loop in the band mark it as parallel. If
   * the loop is removed by the clast, then recursively mark the inner loop in
   * the given band until atleast one loop is marked or no parallel loops are
   * found. Note that the number of bands is equal to the number of dominating
   * parallel loops. */
  for (unsigned k = 0; k < npbands; k++) {
    /* Iterate over loops in a band */
    Band *band = pbands[k];
    Ploop *curr_loop = pbands[k]->loop;
    do {
      char iter[13];
      sprintf(iter, "t%d", curr_loop->depth + 1);
      int *stmtids = (int *)malloc(curr_loop->nstmts * sizeof(int));
      unsigned max_depth = 0;
      for (unsigned j = 0; j < curr_loop->nstmts; j++) {
        Stmt *stmt = curr_loop->stmts[j];
        if (stmt->trans->nrows > max_depth)
          max_depth = stmt->trans->nrows;
        stmtids[j] = stmt->id + 1;
      }

      IF_DEBUG(printf("\tLooking for loop\n"););
      IF_DEBUG(pluto_loop_print(curr_loop););

      unsigned nloops, nstmts;
      int *stmts;
      struct clast_for **loops;
      ClastFilter filter = {iter, stmtids, (int)curr_loop->nstmts, subset};
      clast_filter(root, filter, &loops, (int *)&nloops, &stmts,
                   (int *)&nstmts);

      /* There should be at least one loop.*/
      if (nloops == 0) {
        /* Sometimes loops may disappear (1) tile size larger than trip count
         * 2) it's a scalar dimension but can't be determined from the
         * trans matrix. */

        /* Reiterate over inner loops in the band. */
        unsigned ninloops;
        Ploop **inloops =
            pluto_get_loops_immediately_inner(curr_loop, prog, &ninloops);

        /* TODO: Replace this. Ideally a cloog clast has to be created for each
         * loop till it succeeds. Hence a worklist based algrithm is necessary.
         */
        if (ninloops == 0 ||
            inloops[0]->depth > band->loop->depth + band->width) {
          /* No parallel loops in this band. */
          printf("Warning: parallel poly loop not found in AST\n");
          curr_loop = NULL;
        } else {
          curr_loop = inloops[0];
        }
      } else {
        for (unsigned j = 0; j < nloops; j++) {
          loops[j]->parallel = CLAST_PARALLEL_NOT;
          char *private_vars = (char *)malloc(512);
          strcpy(private_vars, "lbv,ubv");
          if (options->parallel) {
            IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
            loops[j]->parallel = CLAST_PARALLEL_OMP;
            unsigned depth = curr_loop->depth + 1;
            depth++;
            for (; depth <= max_depth; depth++) {
              sprintf(private_vars + strlen(private_vars), ",t%u", depth);
            }
          }
          loops[j]->private_vars = strdup(private_vars);
          free(private_vars);
        }
        /* A loop has been marked parallel. Move to the next band. */
        free(stmtids);
        free(loops);
        free(stmts);
        break;
      }
      free(stmtids);
      free(loops);
      free(stmts);
    } while (curr_loop != NULL);
  }

  pluto_loops_free(ploops, nploops);
  pluto_bands_free(pbands, npbands);
}

/// Marks loops in the cloog clast for unroll jam. Currently it marks all
/// possible candidates. A way to restrict it must be implemented.
void pluto_mark_unroll_jam(struct clast_stmt *root, const PlutoProg *prog,
                           CloogOptions *cloogOptions, unsigned ufactor) {

  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;
  assert(root != NULL);
  unsigned num_ujloops;
  Ploop **ujloops = pluto_get_unroll_jam_loops(prog, &num_ujloops);
  if (num_ujloops == 0) {
    if (!options->silent) {
      printf("[pluto-unroll-jam] No unroll jam loop candidates found\n");
    }
    return;
  }
  IF_DEBUG(printf("Unroll jam candidates \n"));
  IF_DEBUG(pluto_loops_print(ujloops, num_ujloops));

  for (unsigned i = 0; i < num_ujloops; i++) {
    char iter[13];
    sprintf(iter, "t%d", ujloops[i]->depth + 1);
    unsigned *stmtids =
        (unsigned *)malloc(ujloops[i]->nstmts * sizeof(unsigned));
    for (unsigned j = 0; j < ujloops[i]->nstmts; j++) {
      stmtids[j] = ujloops[i]->stmts[j]->id + 1;
    }

    struct clast_for **loops;
    unsigned nloops, nstmts;
    /* TODO: clast_filter takes an arugument of type int**. However the usage of
     * this needs to be checked. */
    int *stmts;
    ClastFilter filter = {iter, (int *)stmtids, (int)ujloops[i]->nstmts,
                          subset};
    clast_filter(root, filter, &loops, (int *)&nloops, &stmts, (int *)&nstmts);

    /* There should be at least one */
    if (nloops == 0) {
      /* Sometimes loops may disappear (1) tile size larger than trip count
       * 2) it's a scalar dimension but can't be determined from the
       * trans matrix */
      printf("[pluto] pluto_unroll_jam: WARNING: Unroll-jam poly loop not "
             "found in AST\n");
      free(stmtids);
      free(loops);
      free(stmts);
      continue;
    }
    for (unsigned j = 0; j < nloops; j++) {
      loops[j]->unroll_type += clast_unroll_and_jam;
      loops[j]->ufactor = ufactor;
    }
    free(stmtids);
    free(loops);
    free(stmts);
  }

  pluto_loops_free(ujloops, num_ujloops);
}

/*
 * Clast-based vector loop marking */
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
                       CloogOptions *cloogOptions) {
  unsigned i, j, nloops, nstmts, nploops;
  struct clast_for **loops;
  int *stmts;
  assert(root != NULL);
  PlutoContext *context = prog->context;

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

    ClastFilter filter = {iter, stmtids, (int)ploops[i]->nstmts, subset};
    clast_filter(root, filter, &loops, (int *)&nloops, &stmts, (int *)&nstmts);

    /* There should be at least one */
    if (nloops == 0) {
      /* Sometimes loops may disappear (1) tile size larger than trip count
       * 2) it's a scalar dimension but can't be determined from the
       * trans matrix */
      printf("[pluto] pluto_mark_vector: WARNING: vectorizable poly loop not "
             "found in AST\n");
      free(stmtids);
      free(loops);
      free(stmts);
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
