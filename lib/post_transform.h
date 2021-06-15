/*
 * Pluto: An automatic parallelier and locality optimizer
 *
 * Copyright (C) 2007 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE
 * in the top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#ifndef _POST_TRANSFORM_H
#define _POST_TRANSFORM_H

typedef struct plutoProg PlutoProg;
typedef struct band Band;
typedef struct pLoop Ploop;
typedef struct plutoContext PlutoContext;

#if defined(__cplusplus)
extern "C" {
#endif

int getDeepestNonScalarLoop(PlutoProg *prog);
int pluto_pre_vectorize_band(Band *band, int num_tiling_levels,
                             PlutoProg *prog);
int pluto_detect_mark_register_tile_loops(PlutoProg *prog);

int gen_reg_tile_file(PlutoProg *prog);

int pluto_loop_is_vectorizable(Ploop *loop, PlutoProg *prog);
Band **get_per_stmt_band(Band *band, unsigned *nstmt_bands);

Band **fuse_per_stmt_bands(Band **per_stmt_bands, unsigned nbands,
                           int num_tiled_levels, unsigned *num_fused_bands,
                           PlutoContext *context);
Ploop *get_best_vectorizable_loop(Ploop **loops, int nloops, PlutoProg *prog);
unsigned get_dimensional_reuse(Ploop *loop);
#if defined(__cplusplus)
}
#endif

#endif // _POST_TRANSFORM_H
