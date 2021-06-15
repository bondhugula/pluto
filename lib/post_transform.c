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
#include "post_transform.h"

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "ddg.h"
#include "math_support.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "program.h"
#include "transforms.h"

/// Used to determine if an access has a short stride.  This is used to favour
/// loops with short strides at the innermost level. Loops with short stride
/// have exhibit spatial locality when they are present at the innermost level.
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

int getDeepestNonScalarLoop(PlutoProg *prog) {
  int loop;

  for (loop = prog->num_hyperplanes - 1; loop >= 0; loop--) {
    if (prog->hProps[loop].type != H_SCALAR) {
      break;
    }
  }

  return loop;
}

/// Dimensional reuse of a loop is the total number of accesses that have either
/// spatial or temporal or group temporal reuse. Choosing a large tile size for
/// the vector dimension removed the need for including accesses that have
/// spatial reuse. Group temporal reuse was not considered in the experiments
/// that were performed. Hence, we return the number of accesses that have
/// temporal reuse.
unsigned get_dimensional_reuse(Ploop *loop) {
  return get_num_invariant_accesses(loop);
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
/// of loops marked for register tiling. The support for unroll jam using this
/// routine is deprectated as hprops has become obsolete and orio does not
/// support unroll jamming of loop nests with vector pragmas .
int pluto_detect_mark_register_tile_loops(PlutoProg *prog) {
  int bandStart, bandEnd;
  int numRegTileLoops;
  int loop, i;
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

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
  PlutoOptions *options = prog->context->options;

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

/// Returns the best loop that can be moved to the innermost
/// level. The loop is treated as vectorizable if 1) it is parallel 2) it has
/// all stride 0 or stride 1 accesses (exhibits spatial locality). A weighted
/// sum of these factors is taken into account.
Ploop *get_best_vectorizable_loop(Ploop **loops, int nloops, PlutoProg *prog) {

  PlutoContext *context = prog->context;

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
/// is used to find statements in the band that are share the same intra-tile
/// iterators.
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
    nbands++;
    pluto_loop_free(loop);
  }
  *nstmt_bands = nbands;
  return per_stmt_bands;
}

/// Assumes that all bands in per_stmt_bands have the same width and begin at
/// the same depth. Returns true if the statements in band b1 and b2 are fused
/// in the tile space.
bool are_stmts_fused_in_band(Band **per_stmt_bands, int b1, int b2,
                             int num_tiled_levels) {

  /* Assumptions of this routine. */
  assert(per_stmt_bands[b1]->loop->depth == per_stmt_bands[b2]->loop->depth &&
         per_stmt_bands[b1]->width == per_stmt_bands[b2]->width);

  unsigned band_begin = per_stmt_bands[b1]->loop->depth;
  unsigned band_end = per_stmt_bands[b1]->loop->depth +
                      num_tiled_levels * per_stmt_bands[b1]->width;

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
                           int num_tiled_levels, unsigned *num_fused_bands,
                           PlutoContext *context) {

  if (nbands == 0) {
    *num_fused_bands = 0;
    return NULL;
  }
  /* Band map is a map that maps each band to a band identifier. All bands that
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
      if (band_map[i] == band_map[j])
        continue;
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

/// The routine looks for a parallel loop with the best cost and makes it the
/// innermost loop. This is called when the loop nest is tiled and outer and
/// inner permutable bands are not the same.
static bool make_parallel_loop_innermost(Band *band, unsigned num_tiled_levels,
                                         PlutoProg *prog) {
  PlutoContext *context = prog->context;

  int depth = band->loop->depth + num_tiled_levels * band->width;
  unsigned nloops;
  Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                        depth, prog, &nloops);
  pluto_loops_print(loops, nloops);
  unsigned nploops = 0;
  Ploop **par_loops = (Ploop **)malloc(nloops * sizeof(Ploop *));
  unsigned innermost_loop_depth = band->loop->depth;
  for (unsigned i = 0; i < nloops; i++) {
    if (pluto_loop_is_parallel(prog, loops[i])) {
      par_loops[nploops++] = loops[i];
    }
    if (loops[i]->depth > innermost_loop_depth) {
      innermost_loop_depth = loops[i]->depth;
    }
  }
  if (nploops == 0) {
    return false;
  }

  Ploop *best_loop = get_best_vectorizable_loop(par_loops, nploops, prog);
  if (best_loop && !pluto_loop_is_innermost(best_loop, prog)) {
    IF_DEBUG(printf("[pluto] intra_tile_opt: loop to be made innermost: "););
    pluto_loop_print(best_loop);

    bool move_across_scalar_hyperplanes = true;
    pluto_make_innermost_loop(best_loop, innermost_loop_depth + 1,
                              move_across_scalar_hyperplanes, prog);
    return true;
  }
  return false;
}

/// Optimize the intra-tile loop order for locality and vectorization.
bool pluto_intra_tile_optimize_band(Band *band, int num_tiled_levels,
                                    PlutoProg *prog) {
  /* Band has to be the innermost band as well */
  if (!pluto_is_band_innermost(band, num_tiled_levels, 0)) {
    /* If the loop nest is tiled (but not completely), then move an intratile
     * parallel loop in this band to the innermost level. */
    if (num_tiled_levels > 0) {
      return make_parallel_loop_innermost(band, num_tiled_levels, prog);
    }
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
  unsigned nstmt_bands;
  Band **per_stmt_bands = get_per_stmt_band(band, &nstmt_bands);

  unsigned num_fused_bands;
  Band **ibands =
      fuse_per_stmt_bands(per_stmt_bands, nstmt_bands, num_tiled_levels,
                          &num_fused_bands, prog->context);

  if (num_fused_bands != nstmt_bands) {
    pluto_bands_free(per_stmt_bands, nstmt_bands);
  }
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;
  if (options->debug) {
    printf("Bands for intra tile optimization \n");
    pluto_bands_print(ibands, num_fused_bands);
  }

  bool retval = false;
  for (unsigned i = 0; i < num_fused_bands; i++) {
    Band *band = ibands[i];
    int depth = band->loop->depth + num_tiled_levels * band->width;
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
      unsigned last_level =
          band->loop->depth + (num_tiled_levels + 1) * band->width;
      /* Move loop across scalar hyperplanes only if the loop nest is tiled. If
       * the loop nest is not tiled, then moving across scalar hyperplanes might
       * not be valid in all cases. */
      /* FIXME: Requires further exploration to find the exact reason for doing
       * this. */
      bool move_across_scalar_hyperplanes = false;
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

/// FIXME:This routine is called only when the band is not tiled ?
/// is_tiled: is the band tiled.
bool pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled) {
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
    if (!prog->context->options->silent) {
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
