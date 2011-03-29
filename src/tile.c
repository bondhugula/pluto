#include <stdio.h>
#include <assert.h>

#include "pluto.h"
#include "post_transform.h"

/* Read tile sizes from file tile.sizes */
static int read_tile_sizes(int *tile_sizes, int *l2_tile_size_ratios,
        int num_tile_dims, HyperplaneProperties *hProps, int firstLoop)
{
    FILE *tsfile = fopen("tile.sizes", "r");

    if (!tsfile)    {
        return 0;
    }
    IF_DEBUG(printf("Reading %d tile sizes\n", num_tile_dims););

    if (options->ft >= 0 && options->lt >= 0)   {
        num_tile_dims = options->lt - options->ft + 1;
    }

    int i=0;
    while (i < num_tile_dims && !feof(tsfile))   {
        if (hProps[firstLoop+i].type != H_SCALAR) {
            fscanf(tsfile, "%d", &tile_sizes[i++]);
        }else{
            /* Size set for scalar dimension doesn't matter */
            tile_sizes[i++] = 42;
        }
    }

    if (i < num_tile_dims)  {
        printf("WARNING: not enough tile sizes provided\n");
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

    return 1;
}


/* Updates the statement domains and transformations to represent the new
 * tiled code. A schedule of tiles is created for parallel execution if
 * --parallel is on */
void pluto_tile(PlutoProg *prog)
{
    int tile_sizes[prog->num_hyperplanes];
    int l2_tile_size_ratios[prog->num_hyperplanes];
    int j;

    /* Tiling */

    int outermostBandStart, outermostBandEnd;
    int innermostBandStart, innermostBandEnd;

    getOutermostTilableBand(prog, &outermostBandStart, &outermostBandEnd);
    getInnermostTilableBand(prog, &innermostBandStart, &innermostBandEnd);

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Outermost tilable band: t%d--t%d\n", 
                outermostBandStart, outermostBandEnd);
    }

    int count = outermostBandEnd - outermostBandStart + 1;

    if (!read_tile_sizes(tile_sizes, l2_tile_size_ratios, count, prog->hProps, outermostBandStart)){
        for (j=0; j<prog->num_hyperplanes; j++)   {
            tile_sizes[j] = DEFAULT_L1_TILE_SIZE;
        }
        for (j=0; j<prog->nvar; j++)   {
            /* L2 cache is around 64 times L1 cache */
            /* assuming 2-d - this tile size has to be eight
             * times the L1 tile size; NOTE: 8 and NOT
             * 8*default_tile_size -- there is a cumulative multiply
             * involved here */
            l2_tile_size_ratios[j] = 8;
        }
    }

    /* Now, we are ready to tile */
    if (options->lt >= 0 && options->ft >= 0)   {
        /* User option specified tiling */

        assert(options->ft <= prog->num_hyperplanes-1);
        assert(options->lt <= prog->num_hyperplanes-1);
        assert(options->ft <= options->lt);

        /* L1 tiling */
        tile_scattering_dims(prog, options->ft, 
                options->lt, tile_sizes);

        if (options->l2tile)    {
            tile_scattering_dims(prog, options->ft, 
                    options->lt, tile_sizes);
        }
    }else{
        /* L1 tiling */
        tile_scattering_dims(prog, outermostBandStart, outermostBandEnd, tile_sizes);
        if (options->l2tile)    {
            /* L2 tiling */
            tile_scattering_dims(prog, 
                    outermostBandStart, 
                    outermostBandEnd, 
                    l2_tile_size_ratios);
        }
    }
}


/* Manipulates statement domain and transformation to tile scattering 
 * dimensions from firstD to lastD */
void tile_scattering_dims(PlutoProg *prog, int firstD, int lastD, int *tile_sizes)
{
    int i, j, k, s;
    int depth;

    HyperplaneProperties *hProps = prog->hProps;

    assert(lastD-firstD+1 <= prog->num_hyperplanes);

    int num_tiled_scat_dims = lastD - firstD + 1;

    for (s=0; s<prog->nstmts; s++)    {

        Stmt *stmt = &prog->stmts[s];

        if (stmt->tile) {

            int num_supernodes = num_tiled_scat_dims;

            /* 1. Specify tiles in the original domain. 
             * NOTE: tile shape info comes in here */

            /* 1.2 make space for the supernode in the domain */
            for (k=0; k<num_tiled_scat_dims; k++)    {
                pluto_constraints_add_dim(stmt->domain, 0);
            }

            /* 1.3 specify tile shapes in the original domain */
            for (depth=firstD; depth<=lastD; depth++)    {

                assert(tile_sizes[depth-firstD] != 0);

                // printf("before adding\n");
                // pluto_constraints_print(stdout, stmt->domain);

                /* Add relation b/w tile space variable and intra-tile variables like
                 * 32*xt <= 2t+i <= 32xt + 31 */
                pluto_constraints_add_inequality(stmt->domain, stmt->domain->nrows);

                /* Lower bound */
                for (j=num_supernodes; j<num_supernodes+stmt->dim; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] = 
                        stmt->trans->val[depth][j-num_supernodes];
                }

                stmt->domain->val[stmt->domain->nrows-1][depth-firstD] = 
                    -tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] = 
                    stmt->trans->val[depth][stmt->dim];

                /* Upper bound */
                pluto_constraints_add_inequality(stmt->domain, stmt->domain->nrows);
                for (j=num_supernodes; j<num_supernodes+stmt->dim; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] = 
                        -stmt->trans->val[depth][j-num_supernodes];
                }

                stmt->domain->val[stmt->domain->nrows-1][depth-firstD] 
                    = tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] = 
                    -stmt->trans->val[depth][stmt->dim]+tile_sizes[depth-firstD]-1;

                // printf("after adding\n");
                // pluto_constraints_print(stdout, stmt->domain);
            }


            /* 2. Update the transformation / scattering functions */

            /* 2.1 make space for the tile space scatterings */
            /* 2.1.1 new columns */
            for (k=0; k<num_supernodes; k++)  {
                pluto_matrix_add_col(&stmt->trans, 0);
            }
            /* 2.1.2 new rows - as many scattering dims we are tiling */
            stmt->is_supernode = (bool *) realloc(stmt->is_supernode,
                    (stmt->trans->nrows + num_supernodes)*sizeof(bool));
            for (k=0; k<num_tiled_scat_dims; k++)  {
                for (i=stmt->trans->nrows-1; i>=firstD; i--) {
                    for (j=0; j<stmt->trans->ncols; j++) {
                        stmt->trans->val[i+1][j] = stmt->trans->val[i][j];
                    }
                    // stmt->trans_loop_type[i+1] = stmt->trans_loop_type[i];
                    stmt->is_supernode[i+1] = stmt->is_supernode[i];
                    // stmt->trans_band_level[i+1] = stmt->trans_band_level[i];
                }
                for (j=0; j<stmt->trans->ncols; j++) {
                    stmt->trans->val[firstD][j] = 0;
                }
                /* This tile space loop has the same property (parallel, fwd dep, or
                 * seq as the original one */
                // stmt->trans_loop_type[firstD] = stmt->trans_loop_type[firstD+num_tiled_scat_dims];
                // stmt->trans_band_level[firstD] = stmt->trans_band_level[firstD+num_tiled_scat_dims];
                stmt->is_supernode[firstD] = true;
                stmt->trans->nrows++;
            }
            /* 2.2 Duplicate the scatterings for the tile space */
            for (depth=firstD; depth < firstD+num_tiled_scat_dims; depth++) {
                for (j=0; j<stmt->trans->ncols; j++)    {
                    if (depth - firstD == j)   {
                        stmt->trans->val[depth][j] = 1;
                    }else{
                        stmt->trans->val[depth][j] = 0;
                    }
                }
            }

            /* 3. Update is_outer_loop */
            stmt->is_orig_loop = realloc(stmt->is_orig_loop, 
                    sizeof(bool)*(stmt->dim+num_supernodes));
            for (k=0; k<num_supernodes; k++)  {
                for (i=0; i<stmt->dim; i++) {
                    stmt->is_orig_loop[i+1] = stmt->is_orig_loop[i];
                }
                stmt->is_orig_loop[0] = true;
                stmt->dim++;
            }

            stmt->num_tiled_loops += lastD-firstD+1;
        }else{
            /* stmt->tile is zero */

            /* Make space for the tile space scatterings */
            /* no need of new columns */

            /* 2.1.2 new rows - as many scattering dims we are tiling */
            for (k=0; k<num_tiled_scat_dims; k++)  {
                pluto_matrix_add_row(&stmt->trans, stmt->trans->nrows);

                stmt->is_supernode[stmt->trans->nrows+k] = false;
            }
        }
    }

    for (k=0; k<num_tiled_scat_dims; k++)  {
        for (i=prog->num_hyperplanes-1; i>=firstD; i--) {
            hProps[i+1] = hProps[i];
        }
        /* This tile space loop has the same property (parallel, fwd dep, or
         * seq as the original one */
        hProps[firstD] = hProps[firstD+num_tiled_scat_dims];
        hProps[firstD].type = (hProps[firstD+num_tiled_scat_dims].type == H_SCALAR)?
            H_SCALAR:H_TILE_SPACE_LOOP;
        prog->num_hyperplanes++;
    }

    for (i=firstD; i<firstD+num_tiled_scat_dims; i++)   {
        if (hProps[i].type == H_SCALAR)    {
            /* Fix it */
            for (k=0; k<prog->nstmts; k++)    {
                for (j=0; j<prog->stmts[k].trans->ncols; j++)    {
                    prog->stmts[k].trans->val[i][j] = prog->stmts[k].trans->val[i+num_tiled_scat_dims][j];
                }
            }

        }

    }
    // dump_hyperplanes(prog->stmts);

    //print_hyperplane_properties(prog->stmts);
}


/* Transform a band of dimensions to get a wavefront
 * (a wavefront of tiles typically)
 *
 * Return: true if something was done, false otherwise
 */
bool create_tile_schedule(PlutoProg *prog, int firstD, int lastD)
{
    int i, j, depth;
    Stmt *stmt;

    if (firstD == lastD)    return false;

    HyperplaneProperties *hProps = prog->hProps;

    if (hProps[firstD].dep_prop == PIPE_PARALLEL) {
        /* If the first one is PIPE_PARALLEL, we are guaranteed to
         * have at least one more pipe_parallel, otherwise the first
         * one would have been SEQ */
        assert(hProps[firstD].dep_prop == PIPE_PARALLEL || 
                hProps[firstD].dep_prop == PARALLEL);
        for (i=0; i<prog->nstmts; i++)    {
            stmt = &prog->stmts[i];
            /* Parallel loop - do a final tile space transformation for a tile schedule 
             * if necessary */

            /* Multiple degrees of pipelined parallelism? */
            /* take sum of all hyperplanes in the bind */
            if (options->multipipe && lastD-firstD >= 2)    {
                for (j=0; j<stmt->trans->ncols; j++)    {
                    for (depth=firstD+1; depth < firstD+3; depth++) {
                        stmt->trans->val[firstD][j] += stmt->trans->val[depth][j];
                    }
                }

            }else{
                /* Just one degree of pipelined parallelism */
                /* -- take sum of first two hyperplanes */
                for (j=0; j<stmt->trans->ncols; j++)    {
                    stmt->trans->val[firstD][j] += stmt->trans->val[firstD+1][j];
                }
            }
        }

        IF_DEBUG(printf("Created tile schedule "););
        if (options->multipipe && lastD-firstD >= 2)    {
            hProps[firstD].dep_prop = SEQ;
            hProps[firstD+1].dep_prop = PARALLEL;
            hProps[firstD+2].dep_prop = PARALLEL;
            IF_DEBUG(printf("for t%d, t%d, t%d\n", firstD, firstD+1, firstD+2););
        }else{
            hProps[firstD].dep_prop = SEQ;
            hProps[firstD+1].dep_prop = PARALLEL;
            IF_DEBUG(printf("for t%d, t%d\n", firstD, firstD+1););
        }
        return true;
    }else{
        return false;
    }
}



/* Find a *non-trivial* outermost band (at least two tilable hyperplanes) */
/* if no non-trivial band is found - it returns the first loop as the trivial
 * band */
void getOutermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd)
{
    int loop, j;

    HyperplaneProperties *hProps = prog->hProps;

    for (loop=0; loop<prog->num_hyperplanes; loop++) {
        /* Skip scalar dimensions at the top */
        if (hProps[loop].type == H_SCALAR) continue;
        if (hProps[loop].dep_prop == PIPE_PARALLEL 
                || hProps[loop].dep_prop == PARALLEL)    {
            j=loop+1;
            while (j<prog->num_hyperplanes
                    && (hProps[j].dep_prop == PIPE_PARALLEL
                        || hProps[j].dep_prop == PARALLEL)
                    && hProps[j].band_num == hProps[loop].band_num)
                //&& hProps[j].type == H_LOOP)
                j++;

            *bandStart = loop;
            *bandEnd = j-1;

            /* Peel off the scalar dimensions from the end */
            while (*bandEnd >= *bandStart && hProps[*bandEnd].type == H_SCALAR)  
                *bandEnd = *bandEnd - 1;
            /* At least two loops */
            if (*bandEnd - *bandStart >= 1)   return;
        }
    }

    /* If there is no tilable band (at least two consecutive loops with fwd or
     * no deps; using one parallel loop as the tilable band is useful so that
     * its intra-tile loop can be vectorized and the tile space can be OpenMP
     * parallel */
    for (loop=0; loop<prog->num_hyperplanes; loop++) {
        if (hProps[loop].dep_prop == PARALLEL)   { 
            *bandStart = loop;
            *bandEnd = loop;
            return;
        }
    }

    /* No non-trivial band; just return the first loop */
    *bandStart = *bandEnd = 0;
    for (loop=0; loop<prog->num_hyperplanes; loop++)    {
        if (hProps[loop].type != H_SCALAR)  {
            *bandStart = *bandEnd = loop;
            break;
        }
    }
}


/* Find the innermost permutable nest (at least two tilable hyperplanes) */
void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd)
{
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
                    && hProps[j].type == H_LOOP)
                j--;

            if (j<=loop-2)  {
                *bandEnd = loop;
                *bandStart = j+1;
                return;
            }
        }
    }
    *bandStart = *bandEnd = lastloop;
}
