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
#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"

int is_invariant(Stmt *stmt, PlutoAccess *acc, int depth) {
  int *divs;
  PlutoMatrix *newacc = pluto_get_new_access_func(stmt, acc->mat, &divs);
  assert(depth <= (int)newacc->ncols - 1);
  unsigned i;
  for (i = 0; i < newacc->nrows; i++) {
    if (newacc->val[i][depth] != 0)
      break;
  }
  int is_invariant = (i == newacc->nrows);
  pluto_matrix_free(newacc);
  free(divs);
  return is_invariant;
}

#define SHORT_STRIDE 4
int has_spatial_reuse(Stmt *stmt, PlutoAccess *acc, int depth) {
  int *divs;
  PlutoMatrix *newacc = pluto_get_new_access_func(stmt, acc->mat, &divs);
  assert(depth <= (int)newacc->ncols - 1);

  /* Scalars */
  if (newacc->nrows == 0)
    return 0;

  for (int i = 0; i < (int)newacc->nrows - 1; i++) {
    /* No spatial reuse when the func is varying at a non-innermost dim */
    if (newacc->val[i][depth] != 0) {
      pluto_matrix_free(newacc);
      free(divs);
      return 0;
    }
  }

  if (newacc->val[newacc->nrows - 1][depth] >= 1 &&
      newacc->val[newacc->nrows - 1][depth] <= SHORT_STRIDE) {
    pluto_matrix_free(newacc);
    free(divs);
    return 1;
  }

  pluto_matrix_free(newacc);
  free(divs);

  return 0;
}

unsigned get_num_invariant_accesses(Ploop *loop) {
  /* All statements under the loop, all accesses for the statement */
  unsigned ni = 0;
  for (unsigned i = 0; i < loop->nstmts; i++) {
    Stmt *stmt = loop->stmts[i];
    for (int j = 0; j < stmt->nreads; j++) {
      ni += is_invariant(stmt, stmt->reads[j], loop->depth);
    }
    for (int j = 0; j < stmt->nwrites; j++) {
      ni += is_invariant(stmt, stmt->writes[j], loop->depth);
    }
  }
  return ni;
}

unsigned get_num_spatial_accesses(Ploop *loop) {
  /* All statements under the loop, all accesses for the statement */
  unsigned ns = 0;
  for (unsigned i = 0; i < loop->nstmts; i++) {
    Stmt *stmt = loop->stmts[i];
    for (int j = 0; j < stmt->nreads; j++) {
      ns += has_spatial_reuse(stmt, stmt->reads[j], loop->depth);
    }
    for (int j = 0; j < stmt->nwrites; j++) {
      ns += has_spatial_reuse(stmt, stmt->writes[j], loop->depth);
    }
  }
  return ns;
}

unsigned get_num_accesses(Ploop *loop) {
  /* All statements under the loop, all accesses for the statement */
  unsigned ns = 0;
  for (unsigned i = 0; i < loop->nstmts; i++) {
    ns += loop->stmts[i]->nreads + loop->stmts[i]->nwrites;
  }

  return ns;
}

int getDeepestNonScalarLoop(PlutoProg *prog) {
  int loop;

  for (loop = prog->num_hyperplanes - 1; loop >= 0; loop--) {
    if (prog->hProps[loop].type != H_SCALAR) {
      break;
    }
  }

  return loop;
}

/* Check if loop is amenable to straightforward vectorization */
int pluto_loop_is_vectorizable(Ploop *loop, PlutoProg *prog) {
  /* LIMITATION: it is possible (rarely) that a loop is not parallel at this
   * position, but, once made innermost, is parallel. We aren't checking
   * if it would be parallel at its new position
   */
  if (!pluto_loop_is_parallel(prog, loop))
    return 0;
  unsigned a = get_num_accesses(loop);
  unsigned s = get_num_spatial_accesses(loop);
  unsigned t = get_num_invariant_accesses(loop);
  /* Vectorize only if each access has either spatial or temporal
   * reuse */
  /* if accesses haven't been provided, a would be 0 */
  if (a >= 1 && a == s + t)
    return 1;

  return 0;
}

/* Detect upto two loops to register tile (unroll-jam) */
int pluto_detect_mark_unrollable_loops(PlutoProg *prog) {
  int bandStart, bandEnd;
  int numUnrollableLoops;
  int loop, i;

  HyperplaneProperties *hProps = prog->hProps;

  /* Loops to be unroll-jammed come from the innermost tilable band; there
   * is trivially always one hyperplane in this band; discount the last one
   * in this band if it's vectorizable. If the innermost tilable doesn't
   * give two loops to unroll-jam, look for parallel loops from inner to
   * outer to fill up the quota of two */

  getInnermostTilableBand(prog, &bandStart, &bandEnd);

  numUnrollableLoops = 0;

  int lastloop = getDeepestNonScalarLoop(prog);

  IF_DEBUG(fprintf(stdout,
                   "[Pluto post transform] Innermost tilable band: t%d--t%d\n",
                   bandStart + 1, bandEnd + 1));

  for (i = 0; i < prog->num_hyperplanes; i++) {
    prog->hProps[i].unroll = NO_UNROLL;
  }

  /* NOTE: CLooG iterators are t0 to t<num>-1 */

  if (bandEnd == lastloop && bandStart < bandEnd) {
    /* Leave alone the vectorizable loop */
    if (hProps[bandEnd].dep_prop == PARALLEL && options->prevector == 1) {
      for (i = PLMAX(bandEnd - 2, bandStart); i <= bandEnd - 1; i++) {
        if (hProps[i].type == H_TILE_SPACE_LOOP)
          continue;
        prog->hProps[i].unroll = UNROLLJAM;
        numUnrollableLoops++;
      }
    } else {
      if (hProps[bandEnd - 1].type != H_TILE_SPACE_LOOP) {
        hProps[bandEnd - 1].unroll = UNROLLJAM;
        numUnrollableLoops++;
      }
      if (hProps[bandEnd].type != H_TILE_SPACE_LOOP) {
        hProps[bandEnd].unroll = UNROLL;
        numUnrollableLoops++;
      }
    }
  } else {
    /* Can unroll only the last loop of course - leave alone if it's
     * vectorizable  */
    if (hProps[lastloop].dep_prop != PARALLEL || options->prevector == 0) {
      hProps[lastloop].unroll = UNROLL;
      numUnrollableLoops = 1;
    }
  }

  if (numUnrollableLoops < 2) {
    /* Any parallel loop at any level can be unrolled */
    for (loop = bandStart - 1; loop >= 0; loop--) {
      if (hProps[loop].dep_prop == PARALLEL &&
          hProps[loop].type != H_TILE_SPACE_LOOP) {
        hProps[loop].unroll = UNROLLJAM;
        numUnrollableLoops++;
        if (numUnrollableLoops == UNROLLJAM)
          break;
      }
    }
  }

  IF_DEBUG(fprintf(
      stdout, "[Pluto post transform] Detected %d unroll/jammable loops\n\n",
      numUnrollableLoops));

  return numUnrollableLoops;
}

/* Create a .unroll - empty .unroll if no unroll-jammable loops */
int gen_unroll_file(PlutoProg *prog) {
  int i;

  HyperplaneProperties *hProps = prog->hProps;
  FILE *unrollfp = fopen(".unroll", "w");

  if (!unrollfp) {
    printf("Error opening .unroll file for writing\n");
    return -1;
  }

  for (i = 0; i < prog->num_hyperplanes; i++) {
    if (hProps[i].unroll == UNROLL) {
      fprintf(unrollfp, "t%d Unroll %d\n", i + 1, options->ufactor);
    } else if (hProps[i].unroll == UNROLLJAM) {
      fprintf(unrollfp, "t%d UnrollJam %d\n", i + 1, options->ufactor);
    }
  }

  fclose(unrollfp);
  return 0;
}

/// Optimize the intra-tile loop order for locality and vectorization.
int pluto_intra_tile_optimize_band(Band *band, int num_tiled_levels,
                                   PlutoProg *prog) {
  /* Band has to be the innermost band as well */
  if (!pluto_is_band_innermost(band, num_tiled_levels)) {
    return 0;
  }

  unsigned nloops;
  Ploop **loops = pluto_get_loops_under(
      band->loop->stmts, band->loop->nstmts,
      band->loop->depth + num_tiled_levels * band->width, prog, &nloops);

  int max_score = INT_MIN;
  Ploop *best_loop = NULL;
  for (unsigned l = 0; l < nloops; l++) {
    int a, s, t, v, score;
    a = get_num_accesses(loops[l]);
    s = get_num_spatial_accesses(loops[l]);
    t = get_num_invariant_accesses(loops[l]);
    v = pluto_loop_is_vectorizable(loops[l], prog);
    /*
     * Penalize accesses which will have neither spatial nor temporal
     * reuse (i.e., non-contiguous ones); high priority for vectorization
     * (if it's vectorizable it also means everything in it has
     * either spatial or temporal reuse, so there is no big tradeoff)
     * TODO: tune this further
     */
    score = (2 * s + 4 * t + 8 * v - 16 * (a - s - t)) * loops[l]->nstmts;
    /* Using >= since we'll take the last one if all else is the same */
    if (score >= max_score) {
      max_score = score;
      best_loop = loops[l];
    }
    IF_DEBUG(
        printf("[pluto-intra-tile-opt] Score for loop %d: %d\n", l, score));
    IF_DEBUG(pluto_loop_print(loops[l]));
  }

  printf("Best loop\n");
  pluto_loop_print(best_loop);
  if (best_loop && !pluto_loop_is_innermost(best_loop, prog)) {
    IF_DEBUG(printf("[pluto] intra_tile_opt: loop to be made innermost: "););
    IF_DEBUG(pluto_loop_print(best_loop););

    /* The last level in the innermost permutable band. This is true only if the
     * outermost permutable band and innermost permutable band are the same. */
    unsigned last_level =
        band->loop->depth + num_tiled_levels * band->width + band->width;
    pluto_make_innermost_loop(best_loop, last_level, prog);
    pluto_loops_free(loops, nloops);
    return 1;
  }

  pluto_loops_free(loops, nloops);
  return 0;
}

/* is_tiled: is the band tiled */
int pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled) {
  unsigned nbands;
  int retval;
  Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

  retval = 0;
  for (unsigned i = 0; i < nbands; i++) {
    retval |= pluto_intra_tile_optimize_band(bands[i], is_tiled, prog);
  }
  pluto_bands_free(bands, nbands);

  if (retval) {
    /* Detect properties again */
    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);
    if (!options->silent) {
      printf("[pluto] After intra-tile optimize\n");
      pluto_transformations_pretty_print(prog);
    }
  }

  return retval;
}

int get_outermost_parallel_loop(const PlutoProg *prog) {
  int parallel_loop, loop;
  HyperplaneProperties *hProps = prog->hProps;

  parallel_loop = -1;
  for (loop = 0; loop < prog->num_hyperplanes; loop++) {
    if (hProps[loop].dep_prop == PARALLEL && hProps[loop].type != H_SCALAR) {
      parallel_loop = loop;

      // Just the outermost parallel one
      break;
    }
  }

  return parallel_loop;
}
