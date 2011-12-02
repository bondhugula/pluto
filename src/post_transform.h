/*
 * PLuTo: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
int detect_mark_unrollable_loops(PlutoProg *prog);
int pre_vectorize (PlutoProg *prog);
int gen_unroll_file(PlutoProg *prog);

#endif
