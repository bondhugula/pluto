/******************************************************************************
 *               libpluto -  A library version of Pluto                       *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2012 Uday Bondhugula                                         *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License as published by the Free Software Foundation; either               *
 * version 2 of the License, or (at your option) any later version.         *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * A copy of the GNU Lesser General Public Licence can be found in the file
 * `LICENSE.LGPL2' in the top-level directory of this distribution. 
 *
 */
#include "pluto.h"
#include "constraints.h"
#include "candl/candl.h"
#include "program.h"
#include "pluto/libpluto.h"
#include "isl/map.h"

PlutoOptions *options;

void normalize_schedule(PlutoConstraints *sched, Stmt *stmt)
{
    int del, r, c, j, snodes, nrows;
    PlutoConstraints *domain;

    domain = stmt->domain;
    snodes = stmt->dim - stmt->dim_orig;

    while (domain != NULL) {
        del = 0;
        nrows = domain->nrows;
        for (r=0; r<nrows; r++) {
            for (c=0; c<snodes; c++) {
                if (domain->val[r-del][c] != 0) break;
            }
            if (c < snodes) {
                PlutoConstraints *cut = pluto_constraints_select_row(domain, r-del);
                for (j=0; j<stmt->trans->nrows; j++) {
                    pluto_constraints_add_dim(cut, 0);
                }
                pluto_constraints_add(sched, cut);
                pluto_constraints_remove_row(domain, r-del);
                del++;
            }
        }
        domain = domain->next;
    }
}

/*
 * Output schedules are isl relations that have dims in the order
 * isl_dim_out, isl_dim_in, div, param, const
 */
__isl_give isl_union_map *pluto_schedule(isl_union_set *domains, 
        isl_union_map *dependences, 
        PlutoOptions *options_l)
{
    int i;
    isl_ctx *ctx;

    ctx = isl_ctx_alloc();

    options = options_l;

    PlutoProg *prog = pluto_prog_alloc();

    prog->nvar = -1;

    prog->nstmts = isl_union_set_n_set(domains);

    if (prog->nstmts >= 1) {
        prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
    }else{
        prog->stmts = NULL;
    }

    for (i=0; i<prog->nstmts; i++) {
        prog->stmts[i] = NULL;
    }

    extract_stmts(domains, prog->stmts);

    for (i=0; i<prog->nstmts; i++) {
        prog->nvar = PLMAX(prog->nvar, prog->stmts[i]->dim);
    }

    if (prog->nstmts >= 1) {
        Stmt *stmt = prog->stmts[0];
        prog->npar = stmt->domain->ncols - stmt->dim - 1;
    }else prog->npar = 0;

    prog->ndeps = isl_union_map_n_map(dependences);

    prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i] = pluto_dep_alloc();
    }
    extract_deps(prog->deps, prog->ndeps, prog->stmts,
            dependences, CANDL_RAW);

    pluto_auto_transform(prog);

    pluto_transformations_print(prog);
    pluto_detect_transformation_properties(prog);

    pluto_print_hyperplane_properties(prog);

    if (options->tile) {
        pluto_tile(prog);
    }

    /* Extract schedules */
    isl_space *space = isl_space_alloc(ctx, prog->npar, 0, 0);

    isl_union_map *schedules = isl_union_map_empty(space);

    // pluto_stmts_print(stdout, prog->stmts, prog->nstmts);

    for (i=0; i<prog->nstmts; i++) {
        isl_basic_map *bmap;
        isl_map *map;
        Stmt *stmt = prog->stmts[i];
        PlutoConstraints *sched = pluto_stmt_get_schedule(stmt);
        normalize_schedule(sched, stmt);

        bmap = isl_basic_map_from_pluto_constraints(ctx, sched, 
                stmt->domain->ncols-1, stmt->trans->nrows, prog->npar);
        map = isl_map_from_basic_map(bmap);
        schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));

        pluto_constraints_free(sched);
    }

    // pluto_stmts_print(stdout, prog->stmts, prog->nstmts);

    pluto_prog_free(prog);
    isl_ctx_free(ctx);

    return schedules;
}
