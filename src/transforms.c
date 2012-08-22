/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2012 Uday Bondhugula
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

/* Sink statement (domain); depth: 0-indexed */
void pluto_sink_statement(Stmt *stmt, int depth, int val, PlutoProg *prog)
{
    assert(stmt->dim == stmt->domain->ncols-prog->npar-1);

    char iter[3];
    sprintf(iter, "d%d", stmt->dim);

    pluto_stmt_add_dim(stmt, depth, -1, iter, H_SCALAR, prog);

    pluto_constraints_set_var(stmt->domain, depth, val);
    stmt->is_orig_loop[depth] = false;
}

/* Stripmine 'dim'th time dimension of stmt by stripmine factor; use
 * 'supernode' as the name of the supernode in the domain */
void pluto_stripmine(Stmt *stmt, int dim, int factor, char *supernode, PlutoProg *prog)
{
    pluto_stmt_add_dim(stmt, 0, dim, supernode, H_TILE_SPACE_LOOP, prog);

    PlutoConstraints *domain = stmt->domain;

    pluto_constraints_add_inequality(domain);
    domain->val[domain->nrows-1][0] = -factor;
    assert(stmt->trans->ncols == domain->ncols);
    int i;
    for (i=1; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i] = stmt->trans->val[dim+1][i];
    }

    pluto_constraints_add_inequality(domain);
    domain->val[domain->nrows-1][0] = factor;
    assert(stmt->trans->ncols == domain->ncols);
    for (i=1; i<stmt->trans->ncols-1; i++)   {
        domain->val[domain->nrows-1][i] = -stmt->trans->val[dim+1][i];
    }
    domain->val[domain->nrows-1][stmt->trans->ncols] += factor;

}

void pluto_interchange_stmt(Stmt *stmt, int level1, int level2,
        PlutoProg *prog)
{
    int j, tmp;

    for (j=0; j<stmt->trans->ncols; j++)   {
        tmp = stmt->trans->val[level1][j];
        stmt->trans->val[level1][j] = stmt->trans->val[level2][j];
        stmt->trans->val[level2][j] = tmp;
    }
}


void pluto_interchange(PlutoProg *prog, int level1, int level2)
{
    int k;
    HyperplaneProperties hTmp;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    for (k=0; k<nstmts; k++)    {
        pluto_interchange_stmt(stmts[k], level1, level2, prog);
    }

    hTmp = prog->hProps[level1]; 
    prog->hProps[level1] = prog->hProps[level2];
    prog->hProps[level2] = hTmp;
}


void pluto_sink_transformation(Stmt *stmt, int pos, PlutoProg *prog)
{
    int i, npar;

    npar = stmt->domain->ncols - stmt->dim - 1;

    assert(pos <= stmt->trans->nrows);
    assert(stmt->dim + npar + 1 == stmt->domain->ncols);

    /* Stmt should always have a transformation */
    assert(stmt->trans != NULL);

    pluto_matrix_add_row(stmt->trans, pos);

    stmt->hyp_types = realloc(stmt->hyp_types, 
            sizeof(int)*stmt->trans->nrows);
    for (i=stmt->trans->nrows-2; i>=pos; i--) {
        stmt->hyp_types[i+1] = stmt->hyp_types[i];
    }
    stmt->hyp_types[pos] = H_SCALAR;
}


void pluto_make_innermost(Ploop *loop, PlutoProg *prog)
{
    int i, d; 

    for (i=0; i<loop->nstmts; i++) {
        Stmt *stmt = loop->stmts[i];
        for (d=loop->depth; d<stmt->trans->nrows-1; d++) {
            pluto_interchange_stmt(stmt, d, d+1, prog);
        }
    }
}
