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

#include "post_transform.h"

#include "program.h"
#include "transforms.h"

int is_invariant(Stmt *stmt, PlutoAccess *acc, int depth) {
  int *divs;
  PlutoMatrix *newacc = pluto_get_new_access_func(acc->mat, stmt, &divs);
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
  PlutoMatrix *newacc = pluto_get_new_access_func(acc->mat, stmt, &divs);
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

/// Detects up to two loops to register tile (unroll-jam). Returns the number
/// of loops marked for register tiling.
int pluto_detect_mark_register_tile_loops(PlutoProg *prog) {
  int bandStart, bandEnd;
  int numRegTileLoops;
  int loop, i;

  HyperplaneProperties *hProps = prog->hProps;

  /* Loops to be unroll-jammed come from the innermost tilable band; there
   * is trivially always one hyperplane in this band; discount the last one
   * in this band if it's vectorizable. If the innermost tilable doesn't
   * give two loops to unroll-jam, look for parallel loops from inner to
   * outer to fill up the quota of two */

  getInnermostTilableBand(prog, &bandStart, &bandEnd);

  numRegTileLoops = 0;

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
        numRegTileLoops++;
      }
    } else {
      if (hProps[bandEnd - 1].type != H_TILE_SPACE_LOOP) {
        hProps[bandEnd - 1].unroll = UNROLLJAM;
        numRegTileLoops++;
      }
      if (hProps[bandEnd].type != H_TILE_SPACE_LOOP) {
        hProps[bandEnd].unroll = UNROLL;
        numRegTileLoops++;
      }
    }
  } else {
    /* Can unroll only the last loop of course - leave alone if it's
     * vectorizable  */
    if (hProps[lastloop].dep_prop != PARALLEL || options->prevector == 0) {
      hProps[lastloop].unroll = UNROLL;
      numRegTileLoops = 1;
    }
  }

  if (numRegTileLoops < 2) {
    /* Any parallel loop at any level can be unrolled */
    for (loop = bandStart - 1; loop >= 0; loop--) {
      if (hProps[loop].dep_prop == PARALLEL &&
          hProps[loop].type != H_TILE_SPACE_LOOP) {
        hProps[loop].unroll = UNROLLJAM;
        numRegTileLoops++;
        if (numRegTileLoops == UNROLLJAM)
          break;
      }
    }
  }

  IF_DEBUG(fprintf(
      stdout, "[Pluto post transform] Detected %d unroll/jammable loops\n\n",
      numRegTileLoops));

  return numRegTileLoops;
}

/* Create a .regtile - empty .regtile if no unroll-jammable loops. */
int gen_reg_tile_file(PlutoProg *prog) {
  int i;

  HyperplaneProperties *hProps = prog->hProps;
  FILE *unrollfp = fopen(".regtile", "w");

  if (!unrollfp) {
    printf("Error opening .regtile file for writing\n");
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
  unsigned num_new_levels =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;
  /* Band has to be the innermost band as well */
  if (!pluto_is_band_innermost(band, num_tiled_levels, num_new_levels)) {
    return 0;
  }

  unsigned nloops;
  Ploop **loops = pluto_get_loops_under(
      band->loop->stmts, band->loop->nstmts,
      band->loop->depth + num_tiled_levels * band->width + num_new_levels, prog,
      &nloops);

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
    if (score > max_score) {
      max_score = score;
      best_loop = loops[l];
    }
    IF_DEBUG(
        printf("[pluto-intra-tile-opt] Score for loop %d: %d\n", l, score));
    IF_DEBUG(pluto_loop_print(loops[l]));
  }

  if (best_loop && !pluto_loop_is_innermost(best_loop, prog)) {
    IF_DEBUG(printf("[pluto] intra_tile_opt: loop to be made innermost: "););
    IF_DEBUG(pluto_loop_print(best_loop););

    /* The last level in the innermost permutable band. This is true only if the
     * outermost permutable band and innermost permutable band are the same. */
    unsigned last_level = band->loop->depth + num_tiled_levels * band->width +
                          band->width + num_new_levels;
    bool move_across_scalar_hyperplanes = false;
    /* Move loop across scalar hyperplanes only if the loop nest is tiled.*/
    if (num_tiled_levels > 0) {
      move_across_scalar_hyperplanes = true;
    }
    pluto_make_innermost_loop(best_loop, last_level,
                              move_across_scalar_hyperplanes, prog);
    pluto_loops_free(loops, nloops);
    return 1;
  }

  pluto_loops_free(loops, nloops);
  return 0;
}

/// This routine is called only when the band is not tiled ?
// is_tiled: is the band tiled
int pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled) {
  unsigned nbands;
  int retval;
  Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

  retval = 0;

  /* This assertion is required because the num_levels_introduced is set to
   * zero. If this routine is called on a tiled loopnest then
   * num_levels_introduced should be set to the number of scalar hyperplanes
   * introduced by post tile distribution */
  assert(is_tiled == false);
  /* unsigned num_levels_introduced = 0; */

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

/// Any reuse between s1 and s2 at hyperplane depth 'depth'
static int has_reuse(Stmt *s1, Stmt *s2, int depth, PlutoProg *prog) {
  for (int i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (((dep->src == s1->id && dep->dest == s2->id) ||
         (dep->src == s2->id && dep->dest == s1->id)) &&
        dep->satvec[depth]) {
      return 1;
    }
  }

  return 0;
}

/// Returns the DDG after removing the deps satisfied by the outermost scalar
/// dimensions.
Graph *get_ddg_for_outermost_non_scalar_level(PlutoProg *prog, Band *band) {
  pluto_dep_satisfaction_reset(prog);
  Graph *new_ddg = ddg_create(prog);
  for (int i = 0; i < prog->num_hyperplanes; i++) {
    if (!pluto_is_hyperplane_scalar(band->loop->stmts[0], i)) {
      break;
    }
    dep_satisfaction_update(prog, i);
  }
  ddg_update(new_ddg, prog);
  return new_ddg;
}

/// Given a band, the routine returns the first non-scalar depth (among the
/// intra tile iterators) at which some statement has a scalar hyperplane.
unsigned get_first_scalar_hyperplane(const Band *band, const PlutoProg *prog) {
  unsigned last_loop_depth = band->loop->depth + 2 * band->width - 1;
  Ploop *loop = band->loop;
  for (unsigned i = loop->depth + band->width; i <= last_loop_depth; i++) {
    if (pluto_is_depth_scalar(loop, i)) {
      continue;
    }
    for (int j = 0; j < loop->nstmts; j++) {
      if (pluto_is_hyperplane_scalar(loop->stmts[j], i)) {
        return i;
      }
    }
  }
  return last_loop_depth + 1;
}

/// Returns SCCs whose statements lie in the current band. This check is
/// incorrect if the outer and inner permutable bands are not same.
int *get_sccs_in_band(Graph *ddg, Band *band, PlutoProg *prog, int *nsccs) {

  int num = 0;
  int *scc_list = (int *)malloc(ddg->num_sccs * sizeof(int));
  for (int i = 0; i < band->loop->nstmts; i++) {
    scc_list[num++] = band->loop->stmts[i]->scc_id;
  }
  *nsccs = num;
  return scc_list;
}

/// Returns true if the list contains SCCs of different dimensionalities
bool contains_sccs_with_diff_dims(int *scc_list, int nsccs, Graph *ddg) {
  for (int i = 0; i < nsccs - 1; i++) {
    /* we many also want to check if the DDGs are connected */
    if (ddg->sccs[i].max_dim != ddg->sccs[i + 1].max_dim)
      return true;
  }
  return false;
}

/// Returns the first depth in the intra tile iterators at which SCCs have
/// different dimensionalities
unsigned get_depth_with_different_scc_dims(Graph *new_ddg, Band *band,
                                           PlutoProg *prog) {

  /* Assumes that there is only 1 level of tiling */
  unsigned new_levels_introduced =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;

  int last_loop_depth =
      band->loop->depth + 2 * band->width - 1 + new_levels_introduced;

  for (int i = 0; i < (int)band->loop->depth - 1; i++) {
    dep_satisfaction_update(prog, i);
  }

  Graph *ddg = prog->ddg;
  prog->ddg = new_ddg;

  unsigned band_begin = band->loop->depth;
  unsigned band_end =
      band->loop->depth + band->width + band->post_tile_dist_hyp_in_band;
  unsigned inner_band_level = band_end;

  /* Iterate over the inter tile dimensions and keep updating the ddg based on
   * this level. Correspondingly move the intra tile dimension */
  for (int i = band_begin; i < band_end; i++) {
    if (pluto_is_depth_scalar(band->loop, i)) {
      dep_satisfaction_update(prog, i);
      continue;
    }
    ddg_update(new_ddg, prog);
    ddg_compute_scc(prog);
    int nsccs = 0;
    int *scc_list = get_sccs_in_band(new_ddg, band, prog, &nsccs);
    bool diff_dim = contains_sccs_with_diff_dims(scc_list, nsccs, new_ddg);
    if (diff_dim) {
      prog->ddg = ddg;
      return inner_band_level;
    }

    /* Move the intra tile hyperplane to the next hyperplane of type loop */
    for (int j = inner_band_level + 1; j < last_loop_depth; j++) {
      if (!pluto_is_depth_scalar(band->loop, j)) {
        break;
      }
      inner_band_level++;
    }
    dep_satisfaction_update(prog, i);
  }
  prog->ddg = ddg;
  return last_loop_depth + 1;
}

int distribute_band_dim_based(Band *band, PlutoProg *prog, int num_tiled_levels,
                              Band **bands, int nbands) {
  unsigned new_levels_introduced =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;

  int last_loop_depth = band->loop->depth +
                        (num_tiled_levels + 1) * band->width - 1 +
                        new_levels_introduced;

  Graph *new_ddg = get_ddg_for_outermost_non_scalar_level(prog, band);
  IF_DEBUG(pluto_matrix_print(stdout, new_ddg->adj););
  unsigned depth = get_depth_with_different_scc_dims(new_ddg, band, prog);

  if (depth == last_loop_depth + 1) {
    return 0;
  }
  printf("Cutting band at depth %d\n", depth);
  return 0;
}

/// See comments for pluto_post_tile_distribute. num_levels_introduced indicates
/// the number of scalar dimensions introduced by the distributing the intratile
/// iterators of the previous bands.
int pluto_post_tile_distribute_band(Band *band, PlutoProg *prog,
                                    int num_tiled_levels, Band **bands,
                                    int nbands) {
  if (band->loop->nstmts == 1)
    return 0;

  /* Distribute based on dimensionalities */
  int retval =
      distribute_band_dim_based(band, prog, num_tiled_levels, bands, nbands);
  if (retval) {
    return retval;
  }

  /* Distribute if there is no reuse*/
  unsigned new_levels_introduced =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;

  int last_loop_depth = band->loop->depth +
                        (num_tiled_levels + 1) * band->width - 1 +
                        new_levels_introduced;

  // printf("last loop depth %d\n", last_loop_depth);

  /* Find depth to distribute statements */
  int depth = last_loop_depth;

  /* This doesn't strictly check for validity of distribution, but only finds a
   * set
   * of loops from innermost which do not satisfy any inter-statement
   * dependences -- this is sufficient to distribute on the outermost among
   * those. Strictly speaking, should have checked whether there is a cycle
   * of dependences between statements at a given depth and all deeper
   * depths (both loops and scalar dims)
   */

  /* Assumption that the innermost and outerost bands are the same ?  */
  for (depth = band->loop->depth + band->width + new_levels_introduced;
       depth <= last_loop_depth; depth++) {
    if (!pluto_satisfies_inter_stmt_dep(prog, band->loop, depth))
      break;
  }

  IF_DEBUG(pluto_band_print(band););
  IF_DEBUG(printf("band_loop_depth %d, depth %d, last_loop_depth %d\n",
                  band->loop->depth, depth, last_loop_depth););

  /* There are no dimensions such that there inter statment deps are not
   * satisfied */
  if (depth == last_loop_depth + 1) {
    return 0;
  }

  /* Look for the first loop with reuse score = 0 */
  for (; depth <= last_loop_depth; depth++) {
    int rscore = 0;
    for (int i = 0; i < band->loop->nstmts; i++) {
      for (int j = i + 1; j < band->loop->nstmts; j++) {
        rscore +=
            has_reuse(band->loop->stmts[i], band->loop->stmts[j], depth, prog);
      }
    }
    if (rscore == 0)
      break;
  }

  if (depth > last_loop_depth)
    return 0;

  IF_DEBUG(printf("[pluto] post_tile_distribute on band\n\t"););
  IF_DEBUG(pluto_band_print(band););
  IF_DEBUG(printf("distributing at depth %d\n", depth););

  /* Distribute statements */
  pluto_separate_stmts(prog, band->loop->stmts, band->loop->nstmts, depth, 0);

  /* band->post_tile_dist_hyp_out_band++; */
  for (int i = 0; i < nbands; i++) {
    if (depth >= bands[i]->loop->depth + bands[i]->width) {
      bands[i]->post_tile_dist_hyp_out_band++;
    } else {
      bands[i]->post_tile_dist_hyp_in_band++;
    }
    if (options->debug) {
      printf("Band %d, in hyp: %d, out hyp: %d\n", i,
             bands[i]->post_tile_dist_hyp_in_band,
             bands[i]->post_tile_dist_hyp_out_band);
    }
  }

  return 1;
}

/// Distribute statements that are fused in the innermost level of a tile if
/// there is no reuse between them (when valid); this is mainly to take care of
/// cache capacity misses / pollution after index set splitting has been
/// performed using mid-point cutting.
int pluto_post_tile_distribute(PlutoProg *prog, Band **bands, int nbands,
                               int num_tiled_levels) {
  if (num_tiled_levels == 0) {
    return 0;
  }
  int retval = 0;
  for (int i = 0; i < nbands; i++) {
    printf("Distributing band of width %d\n", bands[i]->width);
    retval |= pluto_post_tile_distribute_band(bands[i], prog, num_tiled_levels,
                                              bands, nbands);
  }

  if (retval) {
    /* Detect properties again */
    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);
    if (!options->silent) {
      printf("[pluto] After post-tile distribution\n");
      pluto_transformations_pretty_print(prog);
    }
  }

  return retval;
}
