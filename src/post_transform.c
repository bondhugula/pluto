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

#include "pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"


int is_invariant(Stmt *stmt, PlutoAccess *acc, int depth)
{
    int i, *divs;
    PlutoMatrix *newacc = pluto_get_new_access_func(stmt, acc->mat, &divs);
    assert(depth <= newacc->ncols-1);
    for (i=0; i<newacc->nrows; i++) {
        if (newacc->val[i][depth] != 0) break;
    }
    int is_invariant = (i==newacc->nrows);
    pluto_matrix_free(newacc);
    return is_invariant;
}

#define SHORT_STRIDE 4
int has_spatial_reuse(Stmt *stmt, PlutoAccess *acc, int depth)
{
    int i, *divs;
    PlutoMatrix *newacc = pluto_get_new_access_func(stmt, acc->mat, &divs);
    assert(depth <= newacc->ncols-1);

    /* Scalars */
    if (newacc->nrows == 0) return 0;

    for (i=0; i<newacc->nrows-1; i++) {
        /* No spatial reuse when the func is varying at a non-innermost dim */
        if (newacc->val[i][depth] != 0) {
            pluto_matrix_free(newacc);
            return 0;
        }
    }

    if (newacc->val[newacc->nrows-1][depth] >= 1 &&
            newacc->val[newacc->nrows-1][depth] <= SHORT_STRIDE) {
        pluto_matrix_free(newacc);
        return 1;
    }

    pluto_matrix_free(newacc);

    return 0;
}



int get_num_invariant_accesses(Ploop *loop, PlutoProg *prog)
{
    int i, j, ni;

    /* All statements under the loop, all accesses for the statement */
    ni = 0;
    for (i=0; i<loop->nstmts; i++) {
        Stmt *stmt = loop->stmts[i];
        for (j=0; j<stmt->nreads; j++) {
            ni += is_invariant(stmt, stmt->reads[j], loop->depth);
        }
        for (j=0; j<stmt->nwrites; j++) {
            ni += is_invariant(stmt, stmt->writes[j], loop->depth);
        }
    }
    return ni;
}

int get_num_spatial_accesses(Ploop *loop, PlutoProg *prog)
{
    int i, j, ns;

    /* All statements under the loop, all accesses for the statement */
    ns = 0;
    for (i=0; i<loop->nstmts; i++) {
        Stmt *stmt = loop->stmts[i];
        for (j=0; j<stmt->nreads; j++) {
            ns += has_spatial_reuse(stmt, stmt->reads[j], loop->depth);
        }
        for (j=0; j<stmt->nwrites; j++) {
            ns += has_spatial_reuse(stmt, stmt->writes[j], loop->depth);
        }
    }
    return ns;
}


int get_num_accesses(Ploop *loop, PlutoProg *prog)
{
    int i, ns;

    /* All statements under the loop, all accesses for the statement */
    ns = 0;
    for (i=0; i<loop->nstmts; i++) {
        ns += loop->stmts[i]->nreads + loop->stmts[i]->nwrites;
    }

    return ns;
}


int getDeepestNonScalarLoop(PlutoProg *prog)    
{
    int loop;

    for (loop=prog->num_hyperplanes-1; loop>=0; loop--) {
        if (prog->hProps[loop].type != H_SCALAR)    {
            break;
        } 
    }

    return loop;
}


/* Vectorize first loop in band that meets criteria */
int pluto_pre_vectorize_band(Band *band, int num_tiling_levels, PlutoProg *prog)
{
    int num, l;

    /* Band has to be the innermost band as well */
    if (!pluto_is_band_innermost(band, num_tiling_levels)) return 0;

    Ploop **loops;

    loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts, 
            band->loop->depth + num_tiling_levels*band->width, prog, &num);

    for (l=0; l<num; l++) {
        if (!pluto_loop_is_parallel(prog, loops[l])) continue;
        int s, t, a;
        a = get_num_accesses(loops[l], prog);
        s = get_num_spatial_accesses(loops[l], prog);
        t = get_num_invariant_accesses(loops[l], prog);
        /* Vectorize only if each access has either spatial or temporal
         * reuse */
        /* if accesses haven't been provided, a would be 0 */
        if (a >= 1 && a == s + t) break;
    }

    if (l < num) {
        pluto_make_innermost(loops[l], prog);
        IF_DEBUG(printf("[Pluto] Loop to be vectorized: "););
        IF_DEBUG(pluto_loop_print(loops[l]););
        return 1;
    }

    return 0;
}


int pluto_pre_vectorize(PlutoProg *prog)
{
    int nbands, i;
    Band **bands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    int retval = 0;
    for (i=0; i<nbands; i++) {
        retval |= pluto_pre_vectorize_band(bands[i], 0, prog); 
    }
    if (retval) pluto_transformations_pretty_print(prog);
    pluto_bands_free(bands, nbands);
    return 0;
}


/* Detect upto two loops to register tile (unroll-jam) */
int pluto_detect_mark_unrollable_loops(PlutoProg *prog)   
{
    int bandStart, bandEnd;
    int numUnrollableLoops;
    int loop, i;

    HyperplaneProperties *hProps = prog->hProps;

    /* Loops to be unroll-jammed come from the innermost tilable band; there
     * is trivially always one hyperplane in this band; discount the last one
     * in this band if it's vectorizable. If the innermost tilable doesn't
     * give two loops to unroll-jam, look for parallel loops from inner to
     * outer to fill up the quota of two */

    getInnermostTilableBand(prog, &bandStart, &bandEnd);

    numUnrollableLoops=0;

    int lastloop = getDeepestNonScalarLoop(prog);

    IF_DEBUG(fprintf(stdout, "[Pluto post transform] Innermost tilable band: t%d--t%d\n", 
                bandStart+1, bandEnd+1));

    for (i=0; i<prog->num_hyperplanes; i++) {
        prog->hProps[i].unroll = NO_UNROLL;
    }

    /* NOTE: CLooG iterators are t0 to t<num>-1 */

    if (bandEnd == lastloop && bandStart < bandEnd)   {
        /* Leave alone the vectorizable loop */
        if (hProps[bandEnd].dep_prop == PARALLEL && options->prevector == 1)  {
            for (i=PLMAX(bandEnd-2, bandStart); i<=bandEnd-1; i++)    {
                if (hProps[i].type == H_TILE_SPACE_LOOP)   continue;
                prog->hProps[i].unroll = UNROLLJAM;
                numUnrollableLoops++;
            }
        }else{
            if (hProps[bandEnd-1].type != H_TILE_SPACE_LOOP) {
                hProps[bandEnd-1].unroll = UNROLLJAM;
                numUnrollableLoops++;
            }
            if (hProps[bandEnd].type != H_TILE_SPACE_LOOP) {
                hProps[bandEnd].unroll = UNROLL;
                numUnrollableLoops++;
            }
        }
    }else{
        /* Can unroll only the last loop of course - leave alone if it's
         * vectorizable  */
        if (hProps[lastloop].dep_prop != PARALLEL || options->prevector == 0) {
            hProps[lastloop].unroll = UNROLL;
            numUnrollableLoops=1;
        }
    }

    if (numUnrollableLoops < 2) {
        /* Any parallel loop at any level can be unrolled */
        for (loop=bandStart-1; loop>=0; loop--)    {
            if (hProps[loop].dep_prop == PARALLEL && hProps[loop].type != H_TILE_SPACE_LOOP) {
                hProps[loop].unroll = UNROLLJAM;
                numUnrollableLoops++;
                if (numUnrollableLoops == UNROLLJAM) break;
            }
        }
    }

    IF_DEBUG(fprintf(stdout, 
                "[Pluto post transform] Detected %d unroll/jammable loops\n\n", 
                numUnrollableLoops));

    return numUnrollableLoops;
}


/* Create a .unroll - empty .unroll if no unroll-jammable loops */
int gen_unroll_file(PlutoProg *prog)
{
    int i;

    HyperplaneProperties *hProps = prog->hProps;
    FILE *unrollfp = fopen(".unroll", "w");

    if (!unrollfp)  {
        printf("Error opening .unroll file for writing\n");
        return -1;
    }

    for (i=0; i<prog->num_hyperplanes; i++) {
        if (hProps[i].unroll == UNROLL)  {
            fprintf(unrollfp, "t%d Unroll %d\n", i+1, options->ufactor);
        }else if (hProps[i].unroll == UNROLLJAM)    {
            fprintf(unrollfp, "t%d UnrollJam %d\n", i+1, options->ufactor);
        }
    }

    fclose(unrollfp);
    return 0;
}


/*
 * is_tiled: is band tiled?
 */
int pluto_intra_tile_optimize_band(Band *band, int is_tiled, PlutoProg *prog)
{
    int num, l, max_score;
    Ploop *maxloc;

    /* Band has to be the innermost band as well */
    if (is_tiled) {
        if (prog->num_hyperplanes != band->loop->depth + 2*band->width)
            return 0;
    }else{
        if (prog->num_hyperplanes != band->loop->depth + band->width)
            return 0;
    }

    Ploop **loops;

    if (is_tiled) { 
        loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts, 
                band->loop->depth + band->width, prog, &num);
    }else{
        loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts, 
                band->loop->depth, prog, &num);
    }

    max_score = 0;
    maxloc = NULL;
    for (l=0; l<num; l++) {
        int s, t, score;
        s = get_num_spatial_accesses(loops[l], prog);
        t = get_num_invariant_accesses(loops[l], prog);
        score = (2*s + 4*t)*loops[l]->nstmts;
        /* Using >= since we'll take the last one if all else is the same */
        if (score >= max_score) {
            max_score = score;
            maxloc = loops[l];
        }
    }

    if (max_score >= 1) {
        pluto_make_innermost(maxloc, prog);
        IF_DEBUG(printf("[Pluto] intratile_loc: loop to be made innermost:"););
        IF_DEBUG(pluto_loop_print(maxloc););
        pluto_loops_free(loops, num);
        return 1;
    }

    pluto_loops_free(loops, num);
    return 0;
}


/* is_tiled: is the band tiled */
int pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled)
{
    int i, nbands, retval;
    Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

    retval = 0;
    for (i=0; i<nbands; i++) {
        retval |= pluto_intra_tile_optimize_band(bands[i], is_tiled, prog); 
    }
    pluto_bands_free(bands, nbands);
    return retval;
}




int get_outermost_parallel_loop(const PlutoProg *prog)
{
    int parallel_loop, loop;
    HyperplaneProperties *hProps = prog->hProps;

    parallel_loop = -1;
    for (loop=0; loop<prog->num_hyperplanes; loop++) {
        if (hProps[loop].dep_prop == PARALLEL && hProps[loop].type != H_SCALAR)   {
            parallel_loop = loop;

            // Just the outermost parallel one
            break;
        }
    }

    return parallel_loop;
}


/* Unroll scattering functions - incomplete / not used */
void unroll_phis(PlutoProg *prog, int unroll_dim, int ufactor)
{
    int i, j, k;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int npar = prog->npar;

    /*
     * Change the 'stmts' array from this
     *
     * | stmt 0 | stmt 1 | ... | stmt k |
     *
     * to
     *
     * | stmt 0 | stmt 0 | stmt 0 | stmt 0 | stmt 1 | stmt 1 ... |
     *
     * if you are unrolling by four
     */
    for (i=nstmts-1; i>=0; i--)   {

        Stmt *stmt = stmts[i];

        int unroll[prog->nvar];
        int num_unroll = 0;

        for (j=0; j<prog->nvar; j++)    {
            unroll[j] = 0;
        }

        /* 1.1 which original dimensions to unroll */
        for (j=0; j<stmt->dim; j++)    {
            if (stmt->trans->val[unroll_dim][j] != 0) {
                if (unroll[j] == 0)   {
                    num_unroll++;
                }
                unroll[j] = 1;
            }
        }

        for (k=ufactor-1; k>=0; k--)   {
            // Stmt *zstmt = pluto_stmt_dup(stmts[i]);
            Stmt *zstmt = pluto_stmt_alloc(stmts[i]->dim, stmts[i]->domain, stmts[i]->trans);

            for (j=0; j<zstmt->dim; j++) {
                if (unroll[j])  {

                    pluto_constraints_add_dim(zstmt->domain, zstmt->dim, NULL);
                    /* Just put a dummy iterator name since Cloog will
                     * generate a remapping for this too */
                    sprintf(zstmt->iterators[zstmt->dim], "zU%d", j); // where is memory for iterators allocated?

                    /* Now add rows */

                    /* i = ufactor*x + k */
                    pluto_constraints_zero_row(zstmt->domain, zstmt->domain->nrows);
                    pluto_constraints_zero_row(zstmt->domain, zstmt->domain->nrows+1);

                    zstmt->domain->val[zstmt->domain->nrows][zstmt->dim] = -ufactor;
                    zstmt->domain->val[zstmt->domain->nrows][j] = 1;
                    zstmt->domain->val[zstmt->domain->nrows][zstmt->dim+num_unroll+npar] = -k;

                    zstmt->domain->val[zstmt->domain->nrows+1][zstmt->dim] = ufactor;
                    zstmt->domain->val[zstmt->domain->nrows+1][j] = -1;
                    zstmt->domain->val[zstmt->domain->nrows+1][zstmt->dim+num_unroll+npar] = k;
                    zstmt->domain->nrows += 2;

                    /* 0 <= i - ufactor*x <= ufactor - 1 */
                    pluto_constraints_zero_row(zstmt->domain, zstmt->domain->nrows);
                    pluto_constraints_zero_row(zstmt->domain, zstmt->domain->nrows+1);

                    zstmt->domain->val[zstmt->domain->nrows][zstmt->dim] = -ufactor;
                    zstmt->domain->val[zstmt->domain->nrows][j] = 1;
                    zstmt->domain->val[zstmt->domain->nrows][zstmt->dim+num_unroll+npar] = 0;

                    zstmt->domain->val[zstmt->domain->nrows+1][zstmt->dim] = ufactor;
                    zstmt->domain->val[zstmt->domain->nrows+1][j] = -1;
                    zstmt->domain->val[zstmt->domain->nrows+1][zstmt->dim+num_unroll+npar] = ufactor-1;

                    zstmt->domain->nrows += 2;
                }
            }
            zstmt->dim += num_unroll;

            /* Now update the scatterings */
            /* for the new variable added to the domain */
            pluto_matrix_add_col(zstmt->trans, zstmt->trans->ncols-1);

            /* 1. Align the sub-domains (z-polyhedra) */
            zstmt->trans->val[unroll_dim][zstmt->trans->ncols-1] -= k;

            /* Add a new dimension */
            pluto_matrix_add_row(zstmt->trans, zstmt->trans->nrows);

            zstmt->trans->val[zstmt->trans->nrows-1][zstmt->trans->ncols-1] = k;

            // printf("%d %d \n", i*ufactor+k, zstmt->trans->ncols);


            /* Add the statement to the list of statements */

            stmts[i*ufactor + k] =  zstmt;
        }
    }
    nstmts = nstmts*ufactor;
}
