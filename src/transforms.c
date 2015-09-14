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
#include "program.h"
#include "pluto.h"
#include "constraints.h"
#include "transforms.h"

#include "assert.h"

/* Sink statement; depth: 0-indexed */
void pluto_sink_statement(Stmt *stmt, int depth, int val, PlutoProg *prog)
{
    assert(stmt->dim == stmt->domain->ncols-prog->npar-1);

    char iter[3];
    sprintf(iter, "d%d", stmt->dim);

    pluto_stmt_add_dim(stmt, depth, -1, iter, prog);

    pluto_constraints_set_var(stmt->domain, depth, val);
    stmt->is_orig_loop[depth] = false;
}


void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode, PlutoProg *prog)
{
    pluto_stmt_add_dim(stmt, 0, dim, supernode, prog);

    PlutoConstraints *domain = stmt->domain;

    pluto_constraints_add_inequality(domain, domain->nrows);
    domain->val[domain->nrows-1][0] = -factor;
    assert(stmt->trans->ncols == domain->ncols);
    int i;
    for (i=1; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i] = stmt->trans->val[dim+1][i];
    }

    pluto_constraints_add_inequality(domain, domain->nrows);
    domain->val[domain->nrows-1][0] = factor;
    assert(stmt->trans->ncols == domain->ncols);
    for (i=1; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i] = -stmt->trans->val[dim+1][i];
    }
    domain->val[domain->nrows-1][stmt->trans->ncols] += factor;

}
