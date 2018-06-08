/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"


/* Read tile sizes from file tile.sizes */
static int read_tile_sizes(int *tile_sizes, int *l2_tile_size_ratios,
                           int num_tile_dims, Stmt **stmts, int nstmts, int firstLoop) {
    int i, j;

    FILE *tsfile = fopen("tile.sizes", "r");

    if (!tsfile)    return 0;

    IF_DEBUG(printf("[pluto] Reading %d tile sizes\n", num_tile_dims););

    if (options->ft >= 0 && options->lt >= 0)   {
        num_tile_dims = options->lt - options->ft + 1;
    }

    for (i=0; i < num_tile_dims && !feof(tsfile); i++)   {
        for (j=0; j<nstmts; j++) {
            if (pluto_is_hyperplane_loop(stmts[j], firstLoop+i)) break;
        }
        int loop = (j<nstmts);
        if (loop) {
            fscanf(tsfile, "%d", &tile_sizes[i]);
        } else {
            /* Size set for scalar dimension doesn't matter */
            tile_sizes[i] = 42;
        }
    }

    if (i < num_tile_dims)  {
        printf("WARNING: not enough tile sizes provided\n");
        fclose(tsfile);
        return 0;
    }

    i=0;
    while (i < num_tile_dims && !feof(tsfile))   {
        fscanf(tsfile, "%d", &l2_tile_size_ratios[i++]);
    }

    if (i < num_tile_dims)  {
        if (options->l2tile) printf("WARNING: not enough L2 tile sizes provided; using default\n");
        for (i=0; i<num_tile_dims; i++) {
            l2_tile_size_ratios[i] = 8;
        }
    }

    fclose(tsfile);
    return 1;
}

/*
 * Reschedule a diamond tile
 */
int pluto_diamond_tile_reschedule(PlutoProg *prog) {
    int i, j, tmp, retval;

    retval = 0;

    for (i=0; i<prog->nstmts; i++) {
        if (prog->stmts[i]->evicted_hyp) {
            int evicted_hyp_pos = prog->stmts[i]->evicted_hyp_pos;
            int fl = prog->stmts[i]->dim - prog->stmts[i]->dim_orig;
            PlutoMatrix *evicted_hyp = prog->stmts[i]->evicted_hyp;
            assert(fl + evicted_hyp->ncols == prog->stmts[i]->trans->ncols);
            for (j=0; j<evicted_hyp->ncols; j++) {
                tmp = evicted_hyp->val[0][j];
                evicted_hyp->val[0][j] = prog->stmts[i]->trans->val[fl+evicted_hyp_pos][fl+j];
                prog->stmts[i]->trans->val[fl+evicted_hyp_pos][fl+j] = tmp;
            }
            retval = 1;
        }
    }

    return retval;

}


/* Manipulates statement domain and transformation to tile scattering
 * dimensions from firstD to lastD */
void pluto_tile_band(PlutoProg *prog, Band *band, int *tile_sizes) {
    int j, s;
    int depth, npar;

    npar = prog->npar;

    int firstD = band->loop->depth;
    int lastD = band->loop->depth+band->width-1;

    int num_domain_supernodes[band->loop->nstmts];

    for (s=0; s<band->loop->nstmts; s++) {
        num_domain_supernodes[s] = 0;
    }

    for (depth=firstD; depth<=lastD; depth++) {
        for (s=0; s<band->loop->nstmts; s++) {
            Stmt *stmt = band->loop->stmts[s];
            /* 1. Specify tiles in the original domain.
             * NOTE: tile shape info comes in here */

            /* 1.1 Add additional dimensions */
            char iter[6];
            sprintf(iter, "zT%d", stmt->dim);

            int hyp_type = (stmt->hyp_types[depth + depth - firstD] == H_SCALAR)? H_SCALAR:
                           H_TILE_SPACE_LOOP;

            /* 1.2 Specify tile shapes in the original domain */
            // pluto_constraints_print(stdout, stmt->domain);
            if (hyp_type != H_SCALAR) {
                assert(tile_sizes[depth-firstD] >= 1);
                /* Domain supernodes aren't added for scalar dimensions */
                // printf("S%d dim: %d %d\n", stmt->id+1, stmt->dim, depth-firstD);
                pluto_stmt_add_dim(stmt, num_domain_supernodes[s], depth, iter, hyp_type, prog);
                /* Add relation b/w tile space variable and intra-tile variables like
                 * 32*xt <= 2t+i <= 32xt + 31 */
                /* Lower bound */
                pluto_constraints_add_inequality(stmt->domain);

                for (j=num_domain_supernodes[s]+1; j<stmt->dim+npar; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] =
                        stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }

                stmt->domain->val[stmt->domain->nrows-1][num_domain_supernodes[s]] =
                    -tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] =
                    stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar];

                PlutoConstraints *lb = pluto_constraints_select_row(stmt->domain,
                                       stmt->domain->nrows-1);
                pluto_update_deps(stmt, lb, prog);
                pluto_constraints_free(lb);

                /* Upper bound */
                pluto_constraints_add_inequality(stmt->domain);
                for (j=num_domain_supernodes[s]+1; j<stmt->dim+npar; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] =
                        -stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }

                stmt->domain->val[stmt->domain->nrows-1][num_domain_supernodes[s]]
                    = tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] =
                    -stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar]
                    +tile_sizes[depth-firstD]-1;

                PlutoConstraints *ub = pluto_constraints_select_row(stmt->domain,
                                       stmt->domain->nrows-1);
                pluto_update_deps(stmt, ub, prog);
                pluto_constraints_free(ub);

                num_domain_supernodes[s]++;

                // printf("after adding tile constraints\n");
                // pluto_constraints_print(stdout, stmt->domain);

                // printf("Stmt %d: depth: %d\n", stmt->id+1,depth);
                // pluto_matrix_print(stdout, stmt->trans);

            } else {
                /* Scattering function for tile space iterator is set the
                 * same as its associated domain iterator
                 * Dimension is not a loop; tile it trivially
                 */
                pluto_stmt_add_hyperplane(stmt, H_SCALAR, depth);
                for (j=0; j<stmt->dim+npar+1; j++) {
                    stmt->trans->val[depth][j] =
                        stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }
            }
            stmt->num_tiled_loops++;
            stmt->first_tile_dim = firstD ;
            stmt->last_tile_dim = lastD;
        } /* all statements */
    } /* all scats to be tiled */
}



/*
 * Updates statement domains and transformations to represent the new
 * tiled code. A schedule of tiles is created for parallel execution if
 * --parallel is on. Intra-tile optimization is done as part of this as well.
 */
void pluto_tile(PlutoProg *prog) {
    int nbands, i, j, n_ibands, num_tiled_levels, nloops;
    Band **bands, **ibands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
    IF_DEBUG(printf("[pluto_tile] Outermost tilable bands\n"););
    IF_DEBUG(pluto_bands_print(bands, nbands););
    IF_DEBUG(printf("[pluto_tile] Innermost tilable bands\n"););
    IF_DEBUG(pluto_bands_print(ibands, n_ibands););

    num_tiled_levels = 0;

    /*
     * Create bands for innermost parallel loops to be 1-d tiled
     * if they are not dominated by any band
     */
    Ploop **loops = pluto_get_parallel_loops(prog, &nloops);
    for (i=0; i<nloops; i++) {
        if (pluto_loop_is_innermost(loops[i], prog)) {
            for (j=0; j<nbands; j++) {
                if (is_loop_dominated(loops[i], bands[j]->loop, prog)) break;
            }
            if (j==nbands) {
                bands = realloc(bands, (nbands+1)*sizeof(Band *));
                bands[nbands++] = pluto_band_alloc(loops[i], 1);
            }
        }
    }
    pluto_loops_free(loops, nloops);

    /* Now, we are ready to tile */
    if (options->lt >= 0 && options->ft >= 0)   {
        /* User option specified tiling */

        assert(options->ft <= prog->num_hyperplanes-1);
        assert(options->lt <= prog->num_hyperplanes-1);
        assert(options->ft <= options->lt);

        /* L1 tiling */
        pluto_tile_scattering_dims(prog, bands, nbands, 0);
        num_tiled_levels++;

        if (options->l2tile)    {
            pluto_tile_scattering_dims(prog, bands, nbands, 1);
            num_tiled_levels++;
        }
    } else {
        /* L1 tiling */
        pluto_tile_scattering_dims(prog, bands, nbands, 0);
        num_tiled_levels++;
        if (options->l2tile)    {
            /* L2 tiling */
            pluto_tile_scattering_dims(prog, bands, nbands, 1);
            num_tiled_levels++;
        }
    }


    /* Detect properties after tiling */
    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);

    if (!options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        pluto_transformations_pretty_print(prog);
        /* pluto_print_hyperplane_properties(prog); */
    }

    if (options->lbtile) {
        int retval;
        retval = pluto_diamond_tile_reschedule(prog);

        if (retval) {
            pluto_compute_dep_directions(prog);
            pluto_compute_dep_satisfaction(prog);
            if (!options->silent) {
                printf("[Pluto] After intra_tile reschedule\n");
                pluto_transformations_pretty_print(prog);
            }
        }
    }

    if (options->intratileopt) {
        int retval = 0;
        for (i=0; i<nbands; i++) {
            retval |= pluto_intra_tile_optimize_band(bands[i], num_tiled_levels, prog);
        }
        if (retval) {
            pluto_compute_dep_directions(prog);
            pluto_compute_dep_satisfaction(prog);
            if (!options->silent) {
                printf("[Pluto] After intra-tile optimize\n");
                pluto_transformations_pretty_print(prog);
            }
        }
    }


    /* DEPRECATED: now taken care of in intra_tile_optimize */
#if 0
    if (options->prevector) {
        int retval = 0;
        for (i=0; i<nbands; i++) {
            int num_tiling_levels = options->tile + options->l2tile;
            retval |= pluto_pre_vectorize_band(bands[i], num_tiling_levels, prog);
        }
        if (retval) pluto_detect_transformation_properties(prog);
        if (retval && !options->silent) {
            printf("[Pluto] After pre-vectorize:\n");
            pluto_transformations_pretty_print(prog);
            pluto_print_hyperplane_properties(prog);
        }
    }
#endif

    if (options->parallel) {
        int retval = pluto_create_tile_schedule(prog, bands, nbands);
        if (retval && !options->silent) {
            printf("[Pluto] After tile scheduling:\n");
            pluto_transformations_pretty_print(prog);
        }
    }
    pluto_bands_free(bands, nbands);
    pluto_bands_free(ibands, n_ibands);
}


/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void pluto_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2) {
    int i, j, b;
    int depth;
    int tile_sizes[prog->num_hyperplanes];
    int l2_tile_size_ratios[prog->num_hyperplanes];

    Stmt **stmts = prog->stmts;

    for (j=0; j<prog->num_hyperplanes; j++)   {
        tile_sizes[j] = DEFAULT_L1_TILE_SIZE;
        /* L2 cache is around 64 times L1 cache */
        /* assuming 2-d - this tile size has to be eight
         * times the L1 tile size; NOTE: 8 and NOT
         * 8*default_tile_size -- there is a cumulative multiply
         * involved */
        l2_tile_size_ratios[j] = 8;
    }

    for (b=0; b<nbands; b++) {
        read_tile_sizes(tile_sizes, l2_tile_size_ratios, bands[b]->width,
                        bands[b]->loop->stmts, bands[b]->loop->nstmts, bands[b]->loop->depth);

        if (l2) {
            pluto_tile_band(prog, bands[b], l2_tile_size_ratios);
        } else {
            pluto_tile_band(prog, bands[b], tile_sizes);
        }
    } /* all bands */

    /* Sink everything to the same depth */
    int max = 0, curr;
    for (i=0; i<prog->nstmts; i++) {
        max = PLMAX(stmts[i]->trans->nrows, max);
    }
    for (i=0; i<prog->nstmts; i++) {
        curr = stmts[i]->trans->nrows;
        for (j=curr; j < max; j++) {
            pluto_sink_transformation(stmts[i], stmts[i]->trans->nrows, prog);
        }
    }

    curr = prog->num_hyperplanes;
    for (depth=curr; depth<max; depth++)    {
        pluto_prog_add_hyperplane(prog, depth, H_UNKNOWN);
    }
    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    pluto_detect_hyperplane_types_stmtwise(prog);
}


/* Transform a band of dimensions to get a wavefront
 * (a wavefront of tiles typically)
 *
 * Return: true if something was done, false otherwise
 */
bool pluto_create_tile_schedule_band(PlutoProg *prog, Band *band) {
    int i, j, k, depth;

    /* No need to create tile schedule */
    if (pluto_loop_is_parallel(prog, band->loop))  return false;

    /* A band can have scalar dimensions; it starts from a loop */
    int loop_depths[band->width];
    /* Number of dimensions which are loops for all statements in this
     * band */
    int nloops;

    loop_depths[0] = band->loop->depth;
    nloops = 1;
    for (depth=band->loop->depth+1;
            depth < band->loop->depth + band->width; depth++) {
        for (j=0; j<band->loop->nstmts; j++) {
            if (pluto_is_hyperplane_scalar(band->loop->stmts[j], depth))
                break;
        }
        if (j==band->loop->nstmts) {
            /* All of them are loops */
            loop_depths[nloops++] = depth;
        }
    }

    if (nloops <= 1) {
        /* Band doesn't have at least two dimensions for which all
         * statements have loops at those dimensions */
        return false;
    }

    /* Number of inner parallel dims the wavefront will yield */
    int nip_dims;

    /* loop_depths[0...nloops-1] are the depths for which a tile schedule
     * can be created */
    if (prog->options->multipar) {
        /* Full multi-dimensional wavefront */
        nip_dims = nloops-1;
    } else {
        /* just use the first two to create a tile schedule */
        nip_dims = 1;
    }

    /* Wavefront satisfies all deps, all inner loops in the band will
     * become parallel */
    int first = band->loop->depth;

    if (!options->innerpar) {
        /* Create the wavefront */
        for (i=0; i<band->loop->nstmts; i++)    {
            Stmt *stmt = band->loop->stmts[i];
            for (k=1; k<=nip_dims; k++) {
                for (j=0; j<stmt->trans->ncols; j++)    {
                    stmt->trans->val[first][j] +=
                        stmt->trans->val[loop_depths[k]][j];
                }
            }
        }
    }

    IF_DEBUG(printf("[pluto_create_tile_schedule] Created tile schedule for "););
    IF_DEBUG(printf("t%d to t%d\n", first+1, loop_depths[nip_dims]+1));

    /* Update dependence satisfaction levels (better to do this instead of
     * a complete complex dep satisfaction check since we know that the tile
     * schedule will satisfy the dependence satisfied by all the dimensions
     * that is a sum of) */
    for (i=0; i<prog->ndeps; i++) {
        Dep *dep = prog->deps[i];
        /* satvec s should have been computed */
        if (IS_RAR(dep->type)) continue;
        if (pluto_stmt_is_member_of(prog->stmts[dep->src]->id,
                                    band->loop->stmts, band->loop->nstmts)
                && pluto_stmt_is_member_of(prog->stmts[dep->dest]->id,
                                           band->loop->stmts, band->loop->nstmts)) {
            for (k=1; k<=nip_dims; k++) {
                dep->satvec[first] |= dep->satvec[loop_depths[k]];
                dep->satvec[loop_depths[k]] = 0;
            }
        }
    }
    /* Recompute dep directions ? not needed */

    return true;
}


bool pluto_create_tile_schedule(PlutoProg *prog, Band **bands, int nbands) {
    int i;
    bool retval = 0;

    IF_DEBUG(printf("creating tile schedule for bands: \n"););
    IF_DEBUG(pluto_bands_print(bands, nbands););

    for (i=0; i<nbands; i++) {
        retval |= pluto_create_tile_schedule_band(prog, bands[i]);
    }

    return retval;
}


/* Find the innermost permutable nest (at least two tilable hyperplanes) */
void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd) {
    int loop, j, lastloop;

    HyperplaneProperties *hProps = prog->hProps;

    lastloop = getDeepestNonScalarLoop(prog);

    if (hProps[lastloop].dep_prop == SEQ)   {
        *bandStart = *bandEnd = lastloop;
        return;
    }

    for (loop=prog->num_hyperplanes-1; loop>=0; loop--) {
        if (hProps[loop].type == H_SCALAR)   continue;
        if (hProps[loop].dep_prop == PIPE_PARALLEL
                || hProps[loop].dep_prop == PARALLEL)    {
            j=loop-1;
            while (j >= 0
                    && (hProps[j].dep_prop == PIPE_PARALLEL
                        || hProps[j].dep_prop == PARALLEL)
                    && hProps[j].band_num == hProps[loop].band_num
                    && hProps[j].type == H_LOOP) {
                j--;
            }

            if (j<=loop-2)  {
                *bandEnd = loop;
                *bandStart = j+1;
                return;
            }
        }
    }
    *bandStart = *bandEnd = lastloop;
}



