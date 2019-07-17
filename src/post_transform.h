/*
 * PLuTo: An automatic parallelier and locality optimizer
 *
 * Copyright (C) 2007 Uday Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the
 * top-level directory of this program (`COPYING')
 *
 */

#ifndef _POST_TRANSFORM_H
#define _POST_TRANSFORM_H

#include "pluto.h"

int getDeepestNonScalarLoop(PlutoProg *prog);
int pluto_pre_vectorize_band(Band *band, int num_tiling_levels,
                             PlutoProg *prog);
int pluto_detect_mark_register_tile_loops(PlutoProg *prog);

int gen_reg_tile_file(PlutoProg *prog);
int pluto_post_tile_distribute(PlutoProg *prog, Band **bands, int nbands,
                               int n_tiled_levels);
#endif
