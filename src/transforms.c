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

/* Sink statement; depth: 0-indexed */
void pluto_sink_statement(Stmt *stmt, int depth, int val)
{
    pluto_constraints_add_dim(stmt->domain, depth);
    stmt->is_orig_loop[depth] = 0;
    pluto_constraints_set_var(stmt->domain, depth, val);
}
