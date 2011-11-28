/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2008 Uday Kumar Bondhugula
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
#include <stdio.h>
#include <assert.h>

#include "pluto.h"
#include "post_transform.h"
#include "program.h"


void interchange_scattering_dims(PlutoProg *prog, int level1, int level2)
{
    int k, j, tmp;
    HyperplaneProperties hTmp;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    for (k=0; k<nstmts; k++)    {
        for (j=0; j<stmts[k]->trans->ncols; j++)   {
            tmp = stmts[k]->trans->val[level1][j];
            stmts[k]->trans->val[level1][j] = stmts[k]->trans->val[level2][j];
            stmts[k]->trans->val[level2][j] = tmp;
        }

        // tmp = trans_loop_type[level1];
        // stmts[k].trans_loop_type[level1] = stmts[k].trans_loop_type[level2];
        // stmts[k].trans_loop_type[level2] = tmp;
    }

    hTmp = prog->hProps[level1]; 
    prog->hProps[level1] = prog->hProps[level2];
    prog->hProps[level2] = hTmp;
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


int pre_vectorize(PlutoProg *prog)
{
    int lastloop, loop;

    HyperplaneProperties *hProps = prog->hProps;

    /* find the deepest parallel loop NOT belonging to the outermost band that
     * has been identified for parallelization */
    lastloop = getDeepestNonScalarLoop(prog);

    FILE *vfp = fopen(".vectorize", "w");
    for (loop=lastloop; loop>=0; loop--)    {
        /* This loop will not be a tile space loop */
        // printf("%d\n",  loop);
        if ((hProps[loop].dep_prop == PARALLEL 
                    || hProps[loop].dep_prop == PIPE_PARALLEL_INNER_PARALLEL) 
                && (hProps[loop].type != H_TILE_SPACE_LOOP)) {
            if (!options->silent)   {
                fprintf(stdout, "[Pluto post transform] pre-vectorize: moving dimension t%d in\n", 
                        loop+1);
            }

            /* Move this loop inside to possibly enable compiler's (ICC) 
             * auto-vectorization */
            interchange_scattering_dims(prog,loop,lastloop);

            fprintf(vfp, "t%d\n", lastloop+1);

            break;
        }
    }
    fclose(vfp);

    return 1;
}


/* Detect upto two loops to register tile (unroll-jam) */
int detect_unrollable_loops(PlutoProg *prog)   
{
    int bandStart, bandEnd;
    int numUnrollableLoops;
    int loop, i;

    FILE *unrollfp = fopen(".unroll", "w");

    Stmt **stmts = prog->stmts;
    HyperplaneProperties *hProps = prog->hProps;

    if (!unrollfp)  {
        return -1;
    }

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

    for (i=0; i<stmts[0]->trans->nrows; i++) {
        if (hProps[i].unroll == UNROLL)  {
            fprintf(unrollfp, "t%d Unroll %d\n", i+1, options->ufactor);
        }else if (hProps[i].unroll == UNROLLJAM)    {
            fprintf(unrollfp, "t%d UnrollJam %d\n", i+1, options->ufactor);
        }
    }

    IF_DEBUG(fprintf(stdout, 
                "[Pluto post transform] Detected %d unroll/jammable loops\n\n", 
                numUnrollableLoops));

    fclose(unrollfp);

    return numUnrollableLoops;
}


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

                    pluto_constraints_add_dim(zstmt->domain, zstmt->dim);
                    /* Just put a dummy iterator name since Cloog will
                     * generate a remapping for this too */
                    sprintf(zstmt->iterators[zstmt->dim], "zU%d", j);

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
