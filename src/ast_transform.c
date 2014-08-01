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
#include <stdlib.h>
#include <stdio.h>

#include "pluto.h"
#include "program.h"
#include "ast_transform.h"

#include "cloog/cloog.h"

/*
 * Clast-based parallel/distributed loop marking */
void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, nloops, nstmts, nploops;
    int *stmts;
    assert(root != NULL);

    // int filter[1] = {1};

    Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);
    /* Loops should be in the increasing order of their depths */
    qsort(ploops, nploops, sizeof(Ploop *), pluto_loop_compar);

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);

    IF_DEBUG(printf("[pluto_mark_parallel] parallel loops\n"););
    IF_DEBUG(pluto_loops_print(ploops, nploops););

    FILE *pidefs = NULL;
    if (options->distmem) {
        pidefs = fopen("pi_defs.h", "w");
    }

    // clast_pprint(stdout, root, 0, cloogOptions);

    for (i=0; i<nploops; i++) {
        char iter[5];
        struct clast_for **loops;
        sprintf(iter, "t%d", ploops[i]->depth+1);
        int *stmtids = malloc(ploops[i]->nstmts*sizeof(int));
        int max_depth = 0;
        for (j=0; j<ploops[i]->nstmts; j++) {
            Stmt *stmt = ploops[i]->stmts[j];
            if (stmt->trans->nrows > max_depth) max_depth = stmt->trans->nrows;
            stmtids[j] = stmt->id+1;
        }

        IF_DEBUG(printf("Looking for loop\n"););
        IF_DEBUG(pluto_loop_print(ploops[i]););
        // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

        ClastFilter filter = {iter, stmtids, ploops[i]->nstmts, subset};
        clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

        if (options->distmem) {
            /* The parallel loop shouldn't get separated with distmem */
            if (nloops >= 2)  clast_pprint(stdout, root, 0, cloogOptions);
            assert(nloops<=1);
        }

        /* There should be at least one */
        if (nloops==0) {
            /* Sometimes loops may disappear (1) tile size larger than trip count
             * 2) it's a scalar dimension but can't be determined from the
             * trans matrix */
            printf("Warning: parallel poly loop not found in AST\n");
            continue;
        }else{
            for (j=0; j<nloops; j++) {
                loops[j]->parallel = CLAST_PARALLEL_NOT;
                char *private_vars = malloc(128);
                strcpy(private_vars, "lbv,ubv");

                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_comp");
                }

                if (options->distmem) {
                    IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
                    loops[j]->parallel = CLAST_PARALLEL_MPI;
                    if (options->mpiomp) loops[j]->parallel += CLAST_PARALLEL_OMP;
                }

                if (options->parallel && !options->distmem) {
                    IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
                    loops[j]->parallel = CLAST_PARALLEL_OMP;
                    int depth = ploops[i]->depth+1;
                    for (depth++;depth<=max_depth;depth++) {
                        sprintf(private_vars+strlen(private_vars), ",t%d", depth);
                    }
                    loops[j]->private_vars = strdup(private_vars);
                    free(private_vars);
                }

                if (options->distmem) {
                    fprintf(pidefs, "#define _LB_REPLACE_ME_DISTLOOG%d ", i);
                    pprint_expr(cloogOptions, pidefs, loops[j]->LB);
                    fprintf(pidefs, "\n");
                    fprintf(pidefs, "#define _UB_REPLACE_ME_DISTLOOG%d ", i);
                    pprint_expr(cloogOptions, pidefs, loops[j]->UB);
                    fprintf(pidefs, "\n");
                }
                loops[j]->private_vars = strdup("lbv,ubv");
            }
        }
        free(stmtids);
        free(loops);
        free(stmts);

        if (options->distmem) {
            /* MPI parallelize the copy_out statements */
            /* Time the copy_out statements */
            int count = 0;
            stmtids = malloc(prog->nstmts*sizeof(int));
            for (j=0; j<prog->nstmts; j++) {
                if ((prog->stmts[j]->type == COPY_OUT
                            || prog->stmts[j]->type == FOIFI_COPY_OUT)
                        && pluto_stmt_is_member_of(prog->stmts[j]->parent_compute_stmt->id,
                            ploops[i]->stmts, ploops[i]->nstmts)
                   ) {
                    stmtids[count++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[count-1]);
                }
            }
            /* Get all loops at that depth with copy_out statements */
            ClastFilter filter1 = {iter, stmtids, count, subset};
            clast_filter(root, filter1, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_pack");
                }
                loops[j]->parallel = CLAST_PARALLEL_MPI;
            }
            free(stmtids);
            free(loops);
            free(stmts);

            /* MPI parallelize the sigma statements */
            /* Time the sigma statements */
            count = 0;
            stmtids = malloc(prog->nstmts*sizeof(int));
            for (j=0; j<prog->nstmts; j++) {
                if ((prog->stmts[j]->type == SIGMA)
                        && pluto_stmt_is_member_of(prog->stmts[j]->parent_compute_stmt->id,
                            ploops[i]->stmts, ploops[i]->nstmts)
                   ) {
                    stmtids[count++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[count-1]);
                }
            }
            /* Get all loops at that depth with sigma statements */
            ClastFilter filter2 = {iter, stmtids, count, subset};
            clast_filter(root, filter2, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing and marking sigma %s parallel\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_comm");
                }
                loops[j]->parallel = CLAST_PARALLEL_MPI;
            }
            free(stmtids);
            free(loops);
            free(stmts);

            /* Time the copy_in statements */
            count = 0;
            stmtids = malloc(prog->nstmts*sizeof(int));
            for (j=0; j<prog->nstmts; j++) {
                if ((prog->stmts[j]->type == COPY_IN
                            || prog->stmts[j]->type == FOIFI_COPY_IN)
                        && pluto_stmt_is_member_of(prog->stmts[j]->parent_compute_stmt, 
                            ploops[i]->stmts, ploops[i]->nstmts)
                   ) {
                    stmtids[count++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[count-1]);
                }
            }
            /* Get all loops at that depth with copy_in statements */
            ClastFilter filter3 = {iter, stmtids, count, subset};
            clast_filter(root, filter3, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing copy_in %s\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_unpack");
                }
            }
            free(stmtids);
            free(loops);
            free(stmts);

            /* MPI parallelize the write copy_out statements */
            /* Time the write copy_out statements */
            count = 0;
            stmtids = malloc(prog->nstmts*sizeof(int));
            for (j=0; j<prog->nstmts; j++) {
                if (prog->stmts[j]->type == LW_COPY_OUT
                        && pluto_stmt_is_member_of(prog->stmts[j]->parent_compute_stmt->id,
                            ploops[i]->stmts, ploops[i]->nstmts)
                   ) {
                    stmtids[count++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[count-1]);
                }
            }
            /* Get all loops at that depth with write copy_out statements */
            ClastFilter filter4 = {iter, stmtids, count, subset};
            clast_filter(root, filter4, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing write copy_out %s\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_writeout");
                }
                loops[j]->parallel = CLAST_PARALLEL_MPI;
            }
            free(stmtids);
            free(loops);
            free(stmts);

            /* Time the write copy_in statements */
            count = 0;
            stmtids = malloc(prog->nstmts*sizeof(int));
            for (j=0; j<prog->nstmts; j++) {
                if (prog->stmts[j]->type == LW_COPY_IN
                        && pluto_stmt_is_member_of(prog->stmts[j]->parent_compute_stmt->id,
                            ploops[i]->stmts, ploops[i]->nstmts)
                   ) {
                    stmtids[count++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[count-1]);
                }
            }
            /* Get all loops at that depth with write copy_in statements */
            ClastFilter filter5 = {iter, stmtids, count, subset};
            clast_filter(root, filter5, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing write copy_in %s\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->commreport) {
                    loops[j]->time_var_name = strdup("t_writeout");
                }
            }
            free(stmtids);
            free(loops);
            free(stmts);

        }
    }

    if (options->distmem) {
        fclose(pidefs);
    }

    pluto_loops_free(ploops, nploops);
}


/*
 * Clast-based vector loop marking */
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, nloops, nstmts, nploops;
    struct clast_for **loops;
    int *stmts;
    assert(root != NULL);

    Ploop **ploops = pluto_get_parallel_loops(prog, &nploops);

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);
    // clast_pprint(stdout, root, 0, cloogOptions);

    for (i=0; i<nploops; i++) {
        /* Only the innermost ones */
        if (!pluto_is_loop_innermost(ploops[i], prog)) continue;

        IF_DEBUG(printf("[pluto_mark_vector] marking loop\n"););
        IF_DEBUG(pluto_loop_print(ploops[i]););
        char iter[5];
        sprintf(iter, "t%d", ploops[i]->depth+1);
        int *stmtids = malloc(ploops[i]->nstmts*sizeof(int));
        for (j=0; j<ploops[i]->nstmts; j++) {
            stmtids[j] = ploops[i]->stmts[j]->id+1;
        }

        // IF_DEBUG(printf("Looking for loop\n"););
        // IF_DEBUG(pluto_loop_print(ploops[i]););
        // IF_DEBUG(printf("%s in \n", iter););
        // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

        ClastFilter filter = {iter, stmtids, ploops[i]->nstmts, subset};
        clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

        /* There should be at least one */
        if (nloops==0) {
            /* Sometimes loops may disappear (1) tile size larger than trip count
             * 2) it's a scalar dimension but can't be determined from the
             * trans matrix */
            printf("Warning: parallel poly loop not found in AST\n");
            continue;
        }
        for (j=0; j<nloops; j++) {
            // printf("Marking %s ivdep\n", loops[j]->iterator);
            loops[j]->parallel += CLAST_PARALLEL_VEC;
        }
        free(stmtids);
        free(loops);
        free(stmts);
    }

    pluto_loops_free(ploops, nploops);
}
