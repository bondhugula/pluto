/******************************************************************************
 *               libpluto -  A library version of Pluto                       *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2012 Uday Bondhugula                                         *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License as published by the Free Software Foundation; either               *
 * version 2 of the License, or (at your option) any later version.           *
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

/*
 * Pluto tiling modifies the domain; move the information from the
 * domain to the schedules (inequalities will be added to the schedules)
 * Polly for LLVM expects domains to remain unchanged
 * (space/dimensionality-wise)
 */
PlutoConstraints *normalize_domain_schedule(Stmt *stmt, PlutoProg *prog)
{
    int del, r, c, j, snodes, nrows;
    PlutoConstraints *domain;

    PlutoConstraints *sched = pluto_stmt_get_schedule(stmt);

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

    for (c=0; c<snodes; c++) {
        pluto_stmt_remove_dim(stmt, 0, prog);
    }
    return sched;
}

/*
 * Output schedules are isl relations that have dims in the order
 * isl_dim_out, isl_dim_in, div, param, const
 */
__isl_give isl_union_map *pluto_schedule(isl_union_set *domains, 
        isl_union_map *dependences, 
        PlutoOptions *options_l)
{
    int i, j, nbands, n_ibands, retval;
    isl_ctx *ctx;
    isl_space *space;

    ctx = isl_union_set_get_ctx(domains);
    space = isl_union_set_get_space(domains);

    // isl_union_set_dump(domains);
    // isl_union_map_dump(dependences);

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
        prog->params = (char **) malloc(sizeof(char *)*prog->npar);
    }else prog->npar = 0;

    for (i=0; i<prog->npar; i++) {
        prog->params[i] = NULL;
    }

    prog->ndeps = 0;
    isl_union_map_foreach_map(dependences, &isl_map_count, &prog->ndeps);

    prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i] = pluto_dep_alloc();
    }
    extract_deps(prog->deps, 0, prog->stmts,
            dependences, CANDL_RAW);

    IF_DEBUG(pluto_prog_print(prog););

    retval = pluto_auto_transform(prog);

    if (retval) {
        /* Failure */
        pluto_prog_free(prog);
        isl_space_free(space);

        if (!options->silent) {
            printf("[libpluto] failure, returning NULL schedules\n");
        }

        return NULL;
    }

    pluto_detect_transformation_properties(prog);

    if (!options->silent) {
        pluto_transformations_print(prog);
        pluto_print_hyperplane_properties(prog);
    }

    Band **bands, **ibands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
    printf("Outermost tilable bands: %d bands\n", nbands);
    pluto_bands_print(bands, nbands);
    printf("Innermost tilable bands: %d bands\n", n_ibands);
    pluto_bands_print(ibands, n_ibands);

    if (options->tile) {
        pluto_tile(prog);
    }else{
        int retval = pluto_intra_tile_optimize(prog, 0); 
        if (retval) {
            /* Detect properties again */
            pluto_detect_transformation_properties(prog);
            if (!options->silent) {
                printf("[Pluto] after intra tile opt\n");
                pluto_transformations_pretty_print(prog);
            }
        }
    }

    if ((options->parallel) && !options->tile && !options->identity)   {
        /* Obtain wavefront/pipelined parallelization by skewing if
         * necessary */
        int nbands;
        Band **bands;
        bands = pluto_get_outermost_permutable_bands(prog, &nbands);
        bool retval = create_tile_schedule(prog, bands, nbands);
        pluto_bands_free(bands, nbands);

        /* If the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user */
        if (retval && !options->silent)   {
            printf("[Pluto] WARNING: pipelined parallelism exists and --tile is not used.\n");
            printf("use --tile for better parallelization \n");
            fprintf(stdout, "[Pluto] After skewing:\n");
            pluto_transformations_pretty_print(prog);
            pluto_print_hyperplane_properties(prog);
        }

    }

    if (options->tile && !options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);
    }


    if (options->parallel && !options->silent) {
        int nploops;
        Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);
        printf("[pluto_mark_parallel] %d parallel loops\n", nploops);
        pluto_loops_print(ploops, nploops);
        printf("\n");
        pluto_loops_free(ploops, nploops);
    }

    // pluto_stmts_print(stdout, prog->stmts, prog->nstmts);

    /* Construct isl_union_map for pluto schedules */
    isl_union_map *schedules = isl_union_map_empty(isl_space_copy(space));
    for (i=0; i<prog->nstmts; i++) {
        isl_basic_map *bmap;
        isl_map *map;
        Stmt *stmt = prog->stmts[i];
        PlutoConstraints *sched = normalize_domain_schedule(stmt, prog);
        // pluto_constraints_print(stdout, sched);

        bmap = isl_basic_map_from_pluto_constraints(ctx, sched,
                sched->ncols - stmt->trans->nrows - prog->npar - 1, 
                stmt->trans->nrows, prog->npar);
        bmap = isl_basic_map_project_out(bmap, isl_dim_in, 0, sched->ncols - stmt->trans->nrows - stmt->domain->ncols);
        char name[20];
        snprintf(name, sizeof(name), "S_%d", i);
        map = isl_map_from_basic_map(bmap);

        /* Copy ids of the original parameter dimensions  */
        for (j=0; j<isl_space_dim(space, isl_dim_param); j++) {
            isl_id *id = isl_space_get_dim_id(space, isl_dim_param, j);
            map = isl_map_set_dim_id(map, isl_dim_param, j, id);
        }

        map = isl_map_set_tuple_name(map, isl_dim_in, name);
        schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));

        pluto_constraints_free(sched);
    }

    pluto_prog_free(prog);
    isl_space_free(space);

    return schedules;
}
