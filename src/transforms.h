/*
 * PLUTO: An automatic parallelizer and locality optimizer
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
#ifndef _TRANSFORMS_H_
#define _TRANSFORMS_H_

#include "pluto.h"

void pluto_sink_statement(Stmt *stmt, int depth, int val, PlutoProg *prog);
void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode, PlutoProg *prog);
void pluto_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2);
void pluto_reschedule_tile(PlutoProg *prog);
void pluto_interchange(PlutoProg *prog, int level1, int level2);
void pluto_sink_transformation(Stmt *stmt, int pos, PlutoProg *prog);
void pluto_make_innermost(Ploop *loop, PlutoProg *prog);

#endif
