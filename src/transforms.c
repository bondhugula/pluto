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
#include "program.h"
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


void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode)
{
    pluto_stmt_add_dim(stmt, 0, dim, supernode);

    PlutoConstraints *domain = stmt->domain;

    int pos;

    pos = domain->nrows;
    pluto_constraints_add_inequality(domain, domain->nrows);
    domain->val[domain->nrows-1][0] = -factor;
    assert(stmt->trans->ncols == domain->ncols);
    int i;
    for (i=0; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i+1] = stmt->trans->val[dim][i];
    }

    pluto_constraints_add_inequality(domain, domain->nrows);
    domain->val[domain->nrows-1][0] = factor;
    assert(stmt->trans->ncols == domain->ncols);
    for (i=0; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i+1] = -stmt->trans->val[dim][i];
    }
    domain->val[domain->nrows-1][stmt->trans->ncols] += factor;

}
