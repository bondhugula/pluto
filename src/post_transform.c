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
#include <string.h>

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
int pluto_detect_mark_register_tile_loops(PlutoProg *prog,
                                          unsigned num_tiled_levels) {
  int bandStart, bandEnd;
  int numRegTileLoops;
  int loop, i;

  HyperplaneProperties *hProps = prog->hProps;

  /* Loops to be unroll-jammed come from the innermost tilable band; there
   * is trivially always one hyperplane in this band; discount the last one
   * in this band if it's vectorizable. If the innermost tilable doesn't
   * give two loops to unroll-jam, look for parallel loops from inner to
   * outer to fill up the quota of two */

  /* getInnermostTilableBand(prog, &bandStart, &bandEnd); */
  Band **ibands;
  unsigned nbands;
  ibands =
      pluto_get_innermost_permutable_bands(prog, num_tiled_levels, &nbands);
  IF_DEBUG(printf("Innermost permutatble bands \n"));
  pluto_bands_print(ibands, nbands);

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

/// Returns the pointer to the best loop that can be moved to the innermost
/// level. The loop is treated as vectorizable if 1) it is parallel 2) has
/// stride 0 or stride 1 accesses (exhibits spatial locality). A weighted sum of
/// these factors is taken into account. Any access that hinders spatial
/// locality is peanalized.
Ploop *get_best_vectorizable_loop(Ploop **loops, int nloops, PlutoProg *prog) {
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
  return best_loop;
}

/// For each statement in the input outermost band, the routine returns a per
/// statement band. The returned band->loop will have a single statement. This
/// is used to find statements in the input band that are distributed in the
/// tile space.
Band **get_per_stmt_band(Band *band, unsigned *nstmt_bands) {
  unsigned nstmts = band->loop->nstmts;
  unsigned nbands = 0;
  Band **per_stmt_bands = (Band **)malloc(nstmts * sizeof(Band *));
  for (int i = 0; i < nstmts; i++) {
    /* Create a Ploop with a single statement i */
    Ploop *loop = pluto_loop_alloc();
    loop->nstmts = 1;
    loop->depth = band->loop->depth;
    loop->stmts = (Stmt **)malloc(sizeof(Stmt *));
    memcpy(loop->stmts, &band->loop->stmts[i], sizeof(Stmt *));

    per_stmt_bands[i] = pluto_band_alloc(loop, band->width);
    per_stmt_bands[i]->post_tile_dist_hyp_in_band =
        band->post_tile_dist_hyp_in_band;
    per_stmt_bands[i]->post_tile_dist_hyp_out_band =
        band->post_tile_dist_hyp_out_band;
    nbands++;
    pluto_loop_free(loop);
  }
  *nstmt_bands = nbands;
  return per_stmt_bands;
}

/// Assumes that all bands in per_stmt_bands have same widths and begin at same
/// depth. Returns true if the statements in band b1 and b2 are fused in the
/// tile space.
bool are_stmts_fused_in_band(Band **per_stmt_bands, int b1, int b2,
                             int num_tiled_levels) {

  /* Assumptions of this routine. */
  assert(per_stmt_bands[b1]->loop->depth == per_stmt_bands[b2]->loop->depth);
  assert(per_stmt_bands[b1]->width == per_stmt_bands[b2]->width);

  unsigned band_begin = per_stmt_bands[b1]->loop->depth;
  unsigned band_end = per_stmt_bands[b1]->loop->depth +
                      num_tiled_levels * per_stmt_bands[b1]->width +
                      per_stmt_bands[b1]->post_tile_dist_hyp_in_band;

  Stmt *stmt1 = per_stmt_bands[b1]->loop->stmts[0];
  Stmt *stmt2 = per_stmt_bands[b2]->loop->stmts[0];

  for (unsigned i = band_begin; i < band_end; i++) {
    if (!pluto_is_depth_scalar(per_stmt_bands[b1]->loop, i)) {
      continue;
    }
    int col1 = stmt1->trans->ncols - 1;
    int col2 = stmt2->trans->ncols - 1;
    if (stmt1->trans->val[i][col1] != stmt2->trans->val[i][col2]) {
      return false;
    }
  }
  return true;
}

/// The routine fuses per statement bands of statements that are not distributed
/// in the tile space. The intra tile loop iterators of all these statements
/// which are not distributed have to be permuted to the inner levels together.
Band **fuse_per_stmt_bands(Band **per_stmt_bands, unsigned nbands,
                           int num_tiled_levels, unsigned *num_fused_bands) {

  /* Band map is map that maps each band to a band identifier. All bands that
   * are completely fused in the inter tile space will have the same band map.
   * Each band initialized to a unique identifier to start with and is updated
   * when fusion decisions are made.*/
  unsigned *band_map = (unsigned *)malloc(nbands * sizeof(unsigned));
  for (unsigned i = 0; i < nbands; i++) {
    band_map[i] = i;
  }

  /* Compute bands that need to be fused and update band_map. */
  unsigned total_bands = 0;
  for (int i = 0; i < nbands; i++) {
    /* If the current band is not fused with the band that was seen earlier,
     * then it has to be put in a new band. If it is fused with a band k such
     * that k<i, then all bands j>i would also be fused with k. Hence it is not
     * necessary to check for any band that is greater than i that can be fused
     * with i. */
    if (band_map[i] < i)
      continue;

    /* This is a new band that has not been fused with anything that was seen
     * before. Update band_map so that a new band is created for the band i. */
    band_map[i] = total_bands++;

    for (int j = i + 1; j < nbands; j++) {
      /* If the statements in which were originally in bands i and j have been
       * fused, then skip. */
      if (band_map[i] == band_map[j]) {
        continue;
      }
      /* If the statements in the bands i and j are not fused, then do not merge
       * the bands. */
      if (!are_stmts_fused_in_band(per_stmt_bands, i, j, num_tiled_levels)) {
        IF_DEBUG(
            printf("Stmts in bands %d and %d are distributed in tile space\n",
                   i, j););
        continue;
      }
      /* Statements are fused in the band*/
      band_map[j] = band_map[i];
    }
  }

  if (total_bands == nbands) {
    *num_fused_bands = total_bands;
    free(band_map);
    return per_stmt_bands;
  }

  /* Fuse bands based on band map. Note that Band map will be sorted, and
   * band_map[i] will be the id of the smallest band with which band[i] has to
   * be fused. */
  IF_DEBUG(printf("Total number of fused bands %d\n", total_bands););
  Band **fused_bands = (Band **)malloc(total_bands * sizeof(Band *));
  unsigned nfbands = 0;

  for (unsigned i = 0; i < nbands; i++) {
    IF_DEBUG(printf("Band %d to be fused with band %d\n", i, band_map[i]););
    if (band_map[i] >= nfbands) {
      fused_bands[nfbands++] = pluto_band_dup(per_stmt_bands[i]);
      continue;
    }
    unsigned band_id = band_map[i];
    Ploop *loop = fused_bands[band_id]->loop;
    int new_num_stmts = loop->nstmts + per_stmt_bands[i]->loop->nstmts;
    loop->stmts = (Stmt **)realloc(loop->stmts, new_num_stmts * sizeof(Stmt *));
    memcpy(loop->stmts + loop->nstmts, per_stmt_bands[i]->loop->stmts,
           per_stmt_bands[i]->loop->nstmts * sizeof(Stmt *));
    fused_bands[band_id]->loop->nstmts = new_num_stmts;
  }
  *num_fused_bands = nfbands;
  free(band_map);
  return fused_bands;
}

/// Optimize the intra-tile loop order for locality and vectorization.
bool pluto_intra_tile_optimize_band(Band *band, int num_tiled_levels,
                                    PlutoProg *prog) {
  unsigned num_new_levels =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;
  /* Band has to be the innermost band as well */
  if (!pluto_is_band_innermost(band, num_tiled_levels, num_new_levels)) {
    return false;
  }

  /* Band may have been distributed in the inter tile space.  If the band is
   * distributed then there can be multiple innner permutable bands. We need to
   * find optimize for each of them separately. First find bands with a single
   * statement in it and then fuse the bands with statements that are fused in
   * the tile space. */
  /* TODO: We need an early bailout condition here. If the number of statements
   * in the band is 1 or the loop nest is not tiled, then we can skip the
   * procedure to find inner most bands corresponding to the input band. */
  unsigned nstmt_bands = 0;
  Band **per_stmt_bands = get_per_stmt_band(band, &nstmt_bands);

  unsigned num_fused_bands = 0;
  Band **ibands = fuse_per_stmt_bands(per_stmt_bands, nstmt_bands,
                                      num_tiled_levels, &num_fused_bands);

  if (num_fused_bands != nstmt_bands) {
    pluto_bands_free(per_stmt_bands, nstmt_bands);
  }
  if (options->debug) {
    printf("Bands for intra tile optimiztion \n");
    pluto_bands_print(ibands, num_fused_bands);
  }

  bool retval = false;
  for (unsigned i = 0; i < num_fused_bands; i++) {
    Band *band = ibands[i];
    int depth =
        band->loop->depth + num_tiled_levels * band->width + num_new_levels;
    IF_DEBUG(printf("Getting loop at depth  %d\n", depth););

    unsigned nloops;
    Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                          depth, prog, &nloops);
    Ploop *best_loop = get_best_vectorizable_loop(loops, nloops, prog);

    if (best_loop && !pluto_loop_is_innermost(best_loop, prog)) {
      IF_DEBUG(printf("[pluto] intra_tile_opt: loop to be made innermost: "););
      IF_DEBUG(pluto_loop_print(best_loop););

      /* The last level in the innermost permutable band. This is true only if
       * the outermost permutable band and innermost permutable band are the
       * same. */
      unsigned last_level = band->loop->depth +
                            (num_tiled_levels + 1) * band->width +
                            num_new_levels;
      bool move_across_scalar_hyperplanes = false;

      /* Move loop across scalar hyperplanes only if the loop nest is tiled.*/
      if (num_tiled_levels > 0) {
        move_across_scalar_hyperplanes = true;
      }
      pluto_make_innermost_loop(best_loop, last_level,
                                move_across_scalar_hyperplanes, prog);
      retval = true;
    }
    pluto_loops_free(loops, nloops);
  }
  pluto_bands_free(ibands, num_fused_bands);
  return retval;
}

/// This routine is called only when the band is not tiled ?
// is_tiled: is the band tiled
bool pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled) {
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

/// Do stamtements s1 and s2 have any reuse at depth 'depth'. Returns true, if
/// there is reuse else returns false.
static bool has_reuse(Stmt *s1, Stmt *s2, int depth, PlutoProg *prog) {
  for (int i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (((dep->src == s1->id && dep->dest == s2->id) ||
         (dep->src == s2->id && dep->dest == s1->id)) &&
        dep->satvec[depth]) {
      return true;
    }
  }

  return false;
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
/// different dimensionalities.
unsigned get_depth_with_different_scc_dims(Graph *new_ddg, Band *band,
                                           PlutoProg *prog) {

  /* Assumes that there is only 1 level of tiling. */
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
   * this level. Correspondingly move the intra tile dimension. */
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
      free(scc_list);
      return inner_band_level;
    }

    /* Move the intra tile hyperplane to the next hyperplane of type loop. */
    for (int j = inner_band_level + 1; j < last_loop_depth; j++) {
      if (!pluto_is_depth_scalar(band->loop, j)) {
        break;
      }
      inner_band_level++;
    }
    dep_satisfaction_update(prog, i);
    free(scc_list);
  }
  prog->ddg = ddg;
  return last_loop_depth + 1;
}

/// Returns true if Scc given by scc_id is present in the list of SCCs.
bool is_scc_in_list(int *scc_list, int nsccs, int scc_id) {
  for (int i = 0; i < nsccs; i++) {
    if (scc_id == scc_list[i])
      return true;
  }
  return false;
}

/// Distribute SCCs in the band based on dimensionalities at level 'level'.
void distribute_sccs_in_band(Graph *new_ddg, Band *band, int level,
                             PlutoProg *prog) {
  int nstmts = prog->nstmts;

  IF_DEBUG(printf("Distributing SCCs at level %d\n", level););
  /* Add scalar hyperplane for all statements in the program. */
  for (int i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
  }
  int nsccs;
  int *scc_list = get_sccs_in_band(new_ddg, band, prog, &nsccs);
  unsigned offset = 0;
  Scc *sccs = new_ddg->sccs;
  int curr_dim = sccs[scc_list[0]].max_dim;

  /* Iterate over SCCs in the ddg.*/
  for (int i = 0; i < new_ddg->num_sccs; i++) {
    if (!is_scc_in_list(scc_list, nsccs, i))
      continue;
    if (curr_dim != sccs[i].max_dim) {
      offset++;
    }
    /* Iterate over each statement in the program. */
    Stmt **stmts = prog->stmts;
    for (int j = 0; j < prog->nstmts; j++) {
      if (stmts[j]->scc_id == i) {
        int col_num = stmts[j]->trans->ncols - 1;
        stmts[j]->trans->val[level][col_num] = offset;
      }
    }
  }
  pluto_prog_add_hyperplane(prog, level, H_SCALAR);
  prog->hProps[level].dep_prop = SEQ;
  return;
}

/// The routine updates feilds of all bands that count the number of post tile
/// distribution hyperplanes inside and outside the band.
void update_bands(Band **bands, int nbands, int depth) {
  for (int i = 0; i < nbands; i++) {
    int band_end_level = bands[i]->loop->depth + bands[i]->width +
                         bands[i]->post_tile_dist_hyp_in_band;
    if (depth >= band_end_level) {
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
}

/// Distribute the band based on dimensionalities. Returns true if the band was
/// distributed; else returns false. The band is distributed at the outermost
/// level(in the inter tile space) at which there are sccs of different
/// dimensionalities. The cut is introduced in the intra tile space.
bool distribute_band_dim_based(Band *band, PlutoProg *prog,
                               int num_tiled_levels, Band **bands, int nbands) {
  unsigned new_levels_introduced =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;

  int last_loop_depth = band->loop->depth +
                        (num_tiled_levels + 1) * band->width - 1 +
                        new_levels_introduced;

  Graph *new_ddg = get_ddg_for_outermost_non_scalar_level(prog, band);
  IF_DEBUG(pluto_matrix_print(stdout, new_ddg->adj););
  unsigned depth = get_depth_with_different_scc_dims(new_ddg, band, prog);

  if (depth == last_loop_depth + 1) {
    return false;
  }
  IF_DEBUG(printf("Cutting band at depth %d\n", depth););
  distribute_sccs_in_band(new_ddg, band, depth, prog);
  update_bands(bands, nbands, depth);
  graph_free(new_ddg);
  return true;
}

/// See comments for pluto_post_tile_distribute. num_levels_introduced indicates
/// the number of scalar dimensions introduced by the distributing the intratile
/// iterators of the previous bands. Returns true if the band is distributed.
bool pluto_post_tile_distribute_band(Band *band, PlutoProg *prog,
                                     int num_tiled_levels, Band **bands,
                                     int nbands) {
  if (band->loop->nstmts == 1)
    return false;

  /* Distribute based on SCC dimensionalities. */
  int retval =
      distribute_band_dim_based(band, prog, num_tiled_levels, bands, nbands);
  if (retval) {
    return true;
  }

  /* Distribute if there is no reuse. */
  unsigned new_levels_introduced =
      band->post_tile_dist_hyp_in_band + band->post_tile_dist_hyp_out_band;

  int last_loop_depth = band->loop->depth +
                        (num_tiled_levels + 1) * band->width - 1 +
                        new_levels_introduced;

  /* Find depth to distribute statements. */
  int depth = last_loop_depth;

  /* This doesn't strictly check for validity of distribution, but only finds a
   * set of loops from innermost which do not satisfy any inter-statement
   * dependences -- this is sufficient to distribute on the outermost among
   * those. Strictly speaking, should have checked whether there is a cycle of
   * dependences between statements at a given depth and all deeper depths (both
   * loops and scalar dims) */

  /* TODO:Assumption that the innermost and outerost bands are the same ?  */
  for (depth = band->loop->depth + band->width + new_levels_introduced;
       depth <= last_loop_depth; depth++) {
    if (!pluto_satisfies_inter_stmt_dep(prog, band->loop, depth))
      break;
  }

  IF_DEBUG(pluto_band_print(band););
  IF_DEBUG(printf("band_loop_depth %d, depth %d, last_loop_depth %d\n",
                  band->loop->depth, depth, last_loop_depth););

  /* There are no dimensions such that there inter statment deps are not
   * satisfied. */
  if (depth == last_loop_depth + 1) {
    return false;
  }

  /* Look for the first loop with reuse score = 0. */
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
    return false;

  IF_DEBUG(printf("[pluto] post_tile_distribute on band\n\t"););
  IF_DEBUG(pluto_band_print(band););
  IF_DEBUG(printf("distributing at depth %d\n", depth););

  /* Distribute statements. */
  pluto_separate_stmts(prog, band->loop->stmts, band->loop->nstmts, depth, 0);

  /* Update bands */
  update_bands(bands, nbands, depth);

  return true;
}

/// Distribute statements that are fused in the innermost level of a tile if
/// there is no reuse between them (when valid); this is mainly to take care of
/// cache capacity misses / pollution after index set splitting has been
/// performed using mid-point cutting. Returns true if post tile distribute
/// introduces a cut after intra tile iterators for any band. If not, then it
/// returns false.
bool pluto_post_tile_distribute(PlutoProg *prog, Band **bands, int nbands,
                                int num_tiled_levels) {
  if (num_tiled_levels == 0) {
    return false;
  }

  bool retval = false;
  for (int i = 0; i < nbands; i++) {
    IF_DEBUG(printf("Distributing band %d of width %d\n", i, bands[i]->width););
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
