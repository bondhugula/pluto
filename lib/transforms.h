/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE.MIT
 * in the top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#ifndef _TRANSFORMS_H_
#define _TRANSFORMS_H_

#include <stdbool.h>

typedef struct band Band;
typedef struct plutoProg PlutoProg;
typedef struct statement Stmt;
typedef struct pLoop Ploop;

void pluto_sink_statement(Stmt *stmt, int depth, int val, PlutoProg *prog);
void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode,
                     PlutoProg *prog);
void pluto_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands,
                                bool l2);
void pluto_reschedule_tile(PlutoProg *prog);
void pluto_interchange(PlutoProg *prog, int level1, int level2);
void pluto_sink_transformation(Stmt *stmt, unsigned pos);
void pluto_make_innermost_loop(Ploop *loop, unsigned last_level,
                               bool move_across_scalar_hyperplanes,
                               PlutoProg *prog);
void pluto_stmt_loop_interchange(Stmt *stmt, int level1, int level2);

#endif
