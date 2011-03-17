/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
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
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
 *
 */
#include "pluto.h"
#include "constraints.h"
#include "transforms.h"

#include "assert.h"

/* Sink statement; depth: 0-indexed */
void pluto_sink_statement(Stmt *stmt, int npar, int depth, int val)
{
    int d;

    assert(stmt->dim == stmt->domain->ncols-npar-1);

    pluto_constraints_add_dim(stmt->domain, depth);
    stmt->is_orig_loop = realloc(stmt->is_orig_loop, (stmt->dim+1)*sizeof(bool));
    for (d=depth; d<stmt->dim; d++) {
        stmt->is_orig_loop[d+1] = stmt->is_orig_loop[d];
    }
    stmt->is_orig_loop[depth] = false;
    pluto_constraints_set_var(stmt->domain, depth, val);

    stmt->dim++;
}
