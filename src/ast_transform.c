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

void get_dist_dims(const PlutoProg *prog, Ploop *loop, int *dist_dims, int *num_dist_dims)
{
    assert(!options->dynschedule);
    *num_dist_dims = 1;
    dist_dims[0] = loop->depth;
    if (options->multi_level_distribution) {
        int i, num_inner_loops = 0;
        Ploop **inner_loops = pluto_get_loops_under(loop->stmts, loop->nstmts, loop->depth, 
                prog, &num_inner_loops);
        Stmt *stmt = loop->stmts[0];
        int dim;
        assert(stmt->last_tile_dim >= loop->depth);
        for (dim=loop->depth+1; dim<=stmt->last_tile_dim; dim++) {
            if (pluto_is_hyperplane_loop(stmt, dim)) {
                for (i=0; i<num_inner_loops; i++) {
                    if (inner_loops[i]->depth == dim) {
                        if (!pluto_loop_is_parallel(prog, inner_loops[i])) {
                            break;
                        }
                    }
                }
                if (i==num_inner_loops) {
                    dist_dims[(*num_dist_dims)++] = dim;
                    break; // for now, distributing only 2 loops at the most
                }
                else {
                    break; // should not distribute this loop, and loops inner to it
                }
            }
        }
        pluto_loops_free(inner_loops,num_inner_loops);
    }
}

void pluto_mark_parallel_writeout(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions )
{
    int i, j, nloops, nstmts;
    int *stmtids = malloc(prog->nstmts * sizeof(int));
    int *stmts;
    struct clast_for **loops;

    FILE *pidefs = NULL;
    pidefs = fopen("pi_defs.h", "a");

    int max_ploop_id = 0;
    for (i=0; i<prog->nstmts; i++) {
        if (prog->stmts[i]->ploop_id > max_ploop_id) {
            max_ploop_id = prog->stmts[i]->ploop_id;
        }
    }

    int pidefs_added[max_ploop_id+1];
    for (i=0; i<=max_ploop_id; i++) pidefs_added[i] = 0;

    /* MPI parallelize the write copy_out statements */
    j = 0;
    Stmt *stmt = NULL;
    while (j<prog->nstmts) {
        for (; j<prog->nstmts; j++) {
            if (prog->stmts[j]->type == LW_COPY_OUT) {
                stmt = prog->stmts[j];
                stmtids[0] = stmt->id+1;
                assert(stmt->id != -1);
                // printf("S%d\n", stmtids[0]);
                break;
            }
        }
        if (j<prog->nstmts) { // lw_copy_out statement found
            j++;
            /* Get all loops containing exactly the lw_copy_out statement */
            ClastFilter filter = {NULL, stmtids, 1, exact};
            clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);
            if (pidefs_added[stmt->ploop_id] == 0) {
                for (i=0; i<nloops; i++) {
                    fprintf(pidefs, "#define _LB_REPLACE_ME_DISTLOOG%d%s ", stmt->ploop_id, loops[i]->iterator);
                    clast_pprint_expr(cloogOptions, pidefs, loops[i]->LB);
                    fprintf(pidefs, "\n");
                    fprintf(pidefs, "#define _UB_REPLACE_ME_DISTLOOG%d%s ", stmt->ploop_id, loops[i]->iterator);
                    clast_pprint_expr(cloogOptions, pidefs, loops[i]->UB);
                    fprintf(pidefs, "\n");

                }
                pidefs_added[stmt->ploop_id] = 1;
            }
            free(loops);
            free(stmts);
        }
    }

    fclose(pidefs);
//        char iter[5];
//        sprintf(iter, "t%d", ploops[i]->depth+1);
//    /* MPI parallelize the DATA INIT statements */
//      /* Time the data_inti statements */
//      int count = 0;
//      for (j=0; j<prog->nstmts; j++) {
//          if ((prog->stmts[j]->type == DATA_DIST_INIT || prog->stmts[j]->type == DATA_DIST_COPY)
//             ) {
//              stmtids[count++] = prog->stmts[j]->id+1;
//              assert(prog->stmts[j]->id != -1);
//              // printf("S%d\n", stmtids[count-1]);
//          }
//      }
//      /* Get all loops at that depth with copy_out statements */
//      ClastFilter filter7 = {iter, stmtids, count, subset};
//      clast_filter(root, filter7, &loops, &nloops, &stmts, &nstmts);
//      for (j=0; j<nloops; j++) {
//          // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
//          // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
//          if (options->commreport) {
//              loops[j]->time = 1;
//              loops[j]->time_name = strdup("t_data_init");
//          }
//          loops[j]->parallel = CLAST_PARALLEL_MPI;
//      }
//      free(loops);
//      free(stmts);
//    }

    free(stmtids);
}

/*
 * Clast-based parallel/distributed loop marking */
void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, d, nloops, nstmts, nploops;
    int *stmts;
    assert(root != NULL);

    // int filter[1] = {1};
    FILE *pidefs = NULL;
	if(options->distmem && options->data_dist && !options->verify_output)
		pidefs = fopen("pi_defs.h", "a");

    Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);
    /* Loops should be in the increasing order of their depths */
    qsort(ploops, nploops, sizeof(Ploop *), pluto_loop_compar);

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);

    IF_DEBUG(printf("[pluto_mark_parallel] parallel loops\n"););
    IF_DEBUG(pluto_loops_print(ploops, nploops););

    // clast_pprint(stdout, root, 0, cloogOptions);

    int *stmtids = malloc(prog->nstmts*sizeof(int));
    for (i=0; i<nploops; i++) {
        char iter[5];
        struct clast_for **loops;
        sprintf(iter, "t%d", ploops[i]->depth+1);
        int max_depth = 0;
        for (j=0; j<ploops[i]->nstmts; j++) {
            Stmt *stmt = ploops[i]->stmts[j];
            if (options->distmem) assert(stmt->ploop_id == i);
            if (stmt->trans->nrows > max_depth) max_depth = stmt->trans->nrows;
            stmtids[j] = stmt->id+1;
        }
        int nstmtids = ploops[i]->nstmts;

        IF_DEBUG(printf("Looking for loop\n"););
        IF_DEBUG(pluto_loop_print(ploops[i]););
        // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

        if (options->distmem) {
            int num_dist_dims = 0;
            int dist_dims[max_depth];
            get_dist_dims(prog, ploops[i], dist_dims, &num_dist_dims);

            for (d = 0; d<num_dist_dims; d++) {
                int dim = dist_dims[d]+1;
                sprintf(iter, "t%d", dim);
                /* Get all loops at that depth with compute statements */
                ClastFilter filter = {iter, stmtids, nstmtids, subset};
                clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

                /* The parallel loop shouldn't get separated with distmem */
                if (nloops >= 2)  clast_pprint(stdout, root, 0, cloogOptions);
                assert(nloops<=1);

                /* There should be at least one */
                if (nloops==0) {
                    /* Sometimes loops may disappear (1) tile size larger than trip count
                     * 2) it's a scalar dimension but can't be determined from the
                     * trans matrix */
                    PLUTO_MESSAGE(printf("[Pluto] Warning: parallel poly loop not found in AST\n"););
                }else{
                    for (j=0; j<nloops; j++) {
                        IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
                        if (num_dist_dims > 1) loops[j]->loop_id = i;
                        loops[j]->parallel = CLAST_PARALLEL_MPI;
                        if (d == 0) { // only for the first distributed dimension
                            if (options->mpiomp) {
                                loops[j]->parallel += CLAST_PARALLEL_OMP;
                                char *private_vars = malloc(512);
                                strcpy(private_vars, "lbv,ubv,_lb_dist,_ub_dist");
                                int depth = ploops[i]->depth+1;
                                for (depth++;depth<=max_depth;depth++) {
                                    sprintf(private_vars+strlen(private_vars), 
                                            ",lbd_t%d,ubd_t%d,t%d", depth, depth, depth);
                                }
                                loops[j]->private_vars = strdup(private_vars);
                                free(private_vars);
                            }
                            if (options->timereport) {
                                loops[j]->time_var_name = strdup("t_comp");
                            }
                        }

                        if(options->distmem && options->data_dist && !options->verify_output){
							fprintf(pidefs, "#define _LB_REPLACE_ME_DISTLOOG%dt%d ", i, dim);
							clast_pprint_expr(cloogOptions, pidefs, loops[j]->LB);
							fprintf(pidefs, "\n");
							fprintf(pidefs, "#define _UB_REPLACE_ME_DISTLOOG%dt%d ", i, dim);
							clast_pprint_expr(cloogOptions, pidefs, loops[j]->UB);
							fprintf(pidefs, "\n");
                        }

                    }
                }
                free(loops);
                free(stmts);
            }


            for (d = 0; d<num_dist_dims; d++) {
                sprintf(iter, "t%d", dist_dims[d]+1);
                /* MPI parallelize the DATA SETUP statements */
                /* Time the data_inti statements */
                int count = 0;
                for (j=0; j<prog->nstmts; j++) {
                    if ((prog->stmts[j]->type == DATA_SETUP)
                            && prog->stmts[j]->ploop_id == i
                       ) {
                        stmtids[count++] = prog->stmts[j]->id+1;
                        assert(prog->stmts[j]->id != -1);
                        // printf("S%d\n", stmtids[count-1]);
                    }
                }
                /* Get all loops at that depth with copy_out statements */
                ClastFilter filter6 = {iter, stmtids, count, subset};
                clast_filter(root, filter6, &loops, &nloops, &stmts, &nstmts);
                for (j=0; j<nloops; j++) {
                    // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
                    // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                    if (options->timereport) {
                        loops[j]->time_var_name = strdup("t_data_init");
                    }
                    // loops[j]->parallel = CLAST_PARALLEL_MPI;
                }
                free(loops);
                free(stmts);

                /* MPI parallelize the DATA INIT statements */
                /* Time the data_inti statements */
                count = 0;
                for (j=0; j<prog->nstmts; j++) {
                    if ((prog->stmts[j]->type == DATA_DIST_INIT || prog->stmts[j]->type == DATA_DIST_COPY)
                            && prog->stmts[j]->ploop_id == i
                       ) {
                        stmtids[count++] = prog->stmts[j]->id+1;
                        assert(prog->stmts[j]->id != -1);
                        // printf("S%d\n", stmtids[count-1]);
                    }
                }
                /* Get all loops at that depth with copy_out statements */
                ClastFilter filter7 = {iter, stmtids, count, subset};
                clast_filter(root, filter7, &loops, &nloops, &stmts, &nstmts);
                for (j=0; j<nloops; j++) {
                    // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
                    // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                    if (options->timereport && d == 0) {
                        loops[j]->time_var_name = strdup("t_data_init");
                    }
                    loops[j]->loop_id = i;
                    loops[j]->parallel = CLAST_PARALLEL_MPI;
                }
                free(loops);
                free(stmts);

                /* MPI parallelize the DATA MANG statements */
                /* Time the data_inti statements */
                count = 0;
                for (j=0; j<prog->nstmts; j++) {
                    if ((prog->stmts[j]->type == DATA_DIST_MANG)
                            && prog->stmts[j]->ploop_id == i
                       ) {
                        stmtids[count++] = prog->stmts[j]->id+1;
                        assert(prog->stmts[j]->id != -1);
                        // printf("S%d\n", stmtids[count-1]);
                    }
                }
                /* Get all loops at that depth with copy_out statements */
                ClastFilter filter8 = {iter, stmtids, count, subset};
                clast_filter(root, filter8, &loops, &nloops, &stmts, &nstmts);
                for (j=0; j<nloops; j++) {
                    // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
                    // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                    if (options->timereport && d == 0) {
                        loops[j]->time_var_name = strdup("t_data_mang");
                    }
                    loops[j]->loop_id = i;
                    loops[j]->parallel = CLAST_PARALLEL_MPI;
                }
                free(loops);
                free(stmts);
            }

            /* MPI parallelize the copy_out statements */
            /* Time the copy_out statements */
            nstmtids = 0;
            for (j=0; j<prog->nstmts; j++) {
                if (prog->stmts[j]->type == COPY_OUT
                        && prog->stmts[j]->ploop_id == i
                   ) {
                    stmtids[nstmtids++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[nstmtids-1]);
                }
            }
            for (d = 0; d<num_dist_dims; d++) {
                sprintf(iter, "t%d", dist_dims[d]+1);
                /* Get all loops at that depth with copy_out statements */
                ClastFilter filter1 = {iter, stmtids, nstmtids, subset};
                clast_filter(root, filter1, &loops, &nloops, &stmts, &nstmts);
                for (j=0; j<nloops; j++) {
                    // printf("Timing and marking copy_out %s parallel\n", loops[j]->iterator);
                    // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                    if (num_dist_dims > 1) loops[j]->loop_id = i;
                    loops[j]->parallel = CLAST_PARALLEL_MPI;
                    if (d == 0) { // only for the first distributed dimension
                        if (options->timereport) {
                            loops[j]->time_var_name = strdup("t_pack");
                        }
                    }
                }
                free(loops);
                free(stmts);
            }

            /* MPI parallelize the sigma statements */
            /* Time the sigma statements */
            nstmtids = 0;
            for (j=0; j<prog->nstmts; j++) {
                if ((prog->stmts[j]->type == SIGMA)
                        && prog->stmts[j]->ploop_id == i
                   ) {
                    stmtids[nstmtids++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[nstmtids-1]);
                }
            }
            for (d = 0; d<num_dist_dims; d++) {
                sprintf(iter, "t%d", dist_dims[d]+1);
                /* Get all loops at that depth with sigma statements */
                ClastFilter filter2 = {iter, stmtids, nstmtids, subset};
                clast_filter(root, filter2, &loops, &nloops, &stmts, &nstmts);
                for (j=0; j<nloops; j++) {
                    // printf("Timing and marking sigma %s parallel\n", loops[j]->iterator);
                    // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                    if (num_dist_dims > 1) loops[j]->loop_id = i;
                    loops[j]->parallel = CLAST_PARALLEL_MPI;
                    if (d == 0) { // only for the first distributed dimension
                        if (options->timereport) {
                            loops[j]->time_var_name = strdup("t_comm");
                        }
                    }
                }
                free(loops);
                free(stmts);
            }

            /* Time the copy_in statements */
            nstmtids = 0;
            for (j=0; j<prog->nstmts; j++) {
                if (prog->stmts[j]->type == COPY_IN
                        && prog->stmts[j]->ploop_id == i
                   ) {
                    stmtids[nstmtids++] = prog->stmts[j]->id+1;
                    assert(prog->stmts[j]->id != -1);
                    // printf("S%d\n", stmtids[nstmtids-1]);
                }
            }
            sprintf(iter, "t%d", ploops[i]->depth+1);
            /* Get all loops at that depth with copy_in statements */
            ClastFilter filter3 = {iter, stmtids, nstmtids, subset};
            clast_filter(root, filter3, &loops, &nloops, &stmts, &nstmts);
            for (j=0; j<nloops; j++) {
                // printf("Timing copy_in %s\n", loops[j]->iterator);
                // clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                if (options->timereport) {
                    loops[j]->time_var_name = strdup("t_unpack");
                }
            }
            free(loops);
            free(stmts);
        }else {
            ClastFilter filter = {iter, stmtids, nstmtids, subset};
            clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

            /* There should be at least one */
            if (nloops==0) {
                /* Sometimes loops may disappear (1) tile size larger than trip count
                 * 2) it's a scalar dimension but can't be determined from the
                 * trans matrix */
                PLUTO_MESSAGE(printf("Warning: parallel poly loop not found in AST\n"););
            }else{
                for (j=0; j<nloops; j++) {
                    loops[j]->parallel = CLAST_PARALLEL_NOT;
                    char *private_vars = malloc(128);
                    strcpy(private_vars, "lbv,ubv");
                    if (options->parallel) {
                        IF_DEBUG(printf("Marking %s parallel\n", loops[j]->iterator););
                        loops[j]->parallel = CLAST_PARALLEL_OMP;
                        int depth = ploops[i]->depth+1;
                        for (depth++;depth<=max_depth;depth++) {
                            sprintf(private_vars+strlen(private_vars), ",t%d", depth);
                        }
                    }
                    loops[j]->private_vars = strdup(private_vars);
                    free(private_vars);
                }
            }
            free(loops);
            free(stmts);
        }
    }
    free(stmtids);

	if(options->distmem && options->data_dist && !options->verify_output){
		fclose(pidefs);
	}
    if (options->distmem) {
        pluto_mark_parallel_writeout(root, prog, cloogOptions);
    }

    pluto_loops_free(ploops, nploops);
}

/*
 * Clast-based parallel/distributed loop marking for distmem-dynschedule code*/
void pluto_mark_parallel_dynschedule(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, l, nloops, nstmts;
    int stmtids[1];
    int *stmts;
    struct clast_for **loops;

    FILE *pidefs = NULL;
    if (!options->distmem) {
        pidefs = fopen("pi_defs.h", "a");
    }

    /* OMP parallelize the all_tasks statements */
    j = 0;
    while (j<prog->nstmts) {
        char iter[5];
        int depth = 0;
        Stmt *stmt = NULL;
        for (; j<prog->nstmts; j++) {
            if (prog->stmts[j]->type == ALL_TASKS) {
                stmt = prog->stmts[j];
                stmtids[0] = stmt->id+1;
                assert(stmt->id != -1);
                // printf("S%d\n", stmtids[0]);
                for (i=0; i<stmt->trans->nrows; i++) {
                    if (stmt->hyp_types[i] != H_SCALAR) break;
                }
                depth = i+1;
                sprintf(iter, "t%d", i+1);
                break;
            }
        }
        if (j<prog->nstmts) { // all_tasks statement found
            j++;
            /* Get all loops at that depth with exactly the all_tasks statement */
            ClastFilter filter = {iter, stmtids, 1, exact};
            clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);
            for (l=0; l<nloops; l++) {
//                 printf("Marking all_tasks %s parallel\n", loops[j]->iterator);
//                 clast_pprint(stdout, loops[j]->body, 0, cloogOptions);
                loops[l]->parallel = CLAST_PARALLEL_OMP;
                char *private_vars = malloc(512);
                strcpy(private_vars, "");
                if (options->distmem) {
                    sprintf(private_vars+strlen(private_vars), "remote_dep_tasks, ");
                }
                sprintf(private_vars+strlen(private_vars), "local_dep_tasks, affinity");
                for (depth++;depth<=stmt->trans->nrows;depth++) {
                    sprintf(private_vars+strlen(private_vars), ", t%d", depth);
                }
                loops[l]->private_vars = strdup(private_vars);
                free(private_vars);
                char *reduction_vars = malloc(128);
                strcpy(reduction_vars, "+:_num_tasks_to_execute");
                if (options->distmem) {
                    sprintf(reduction_vars+strlen(reduction_vars), ",_num_tasks_to_unpack");
                }
                loops[l]->reduction_vars = strdup(reduction_vars);
                free(reduction_vars);
            }
//            if ((nloops>0)) {
//                assert(nloops == 1);
//            }
            free(loops);
            free(stmts);

            if (!options->distmem) {
                /* Get all loops with exactly the all_tasks statement */
                ClastFilter filter1 = {NULL, stmtids, 1, exact};
                clast_filter(root, filter1, &loops, &nloops, &stmts, &nstmts);
                for (l=0; l<nloops; l++) {
                    fprintf(pidefs, "#define _LB_REPLACE_ME_DISTLOOG%d%s ", stmt->ploop_id, loops[l]->iterator);
                    clast_pprint_expr(cloogOptions, pidefs, loops[l]->LB);
                    fprintf(pidefs, "\n");
                    fprintf(pidefs, "#define _UB_REPLACE_ME_DISTLOOG%d%s ", stmt->ploop_id, loops[l]->iterator);
                    clast_pprint_expr(cloogOptions, pidefs, loops[l]->UB);
                    fprintf(pidefs, "\n");
                }
                free(loops);
                free(stmts);
            }
        }
    }

    pluto_mark_parallel_writeout(root, prog, cloogOptions);

    if (!options->distmem) {
        fclose(pidefs);
    }
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

    Ploop **ploops;
//    if(options->data_dist && options->data_tile_opt )
//		ploops = pluto_get_all_loops(prog, &nploops);
//    else
    ploops = pluto_get_parallel_loops(prog, &nploops);

    IF_DEBUG(printf("[pluto_mark_vector] parallel loops\n"););
    IF_DEBUG(pluto_loops_print(ploops, nploops););

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);
    // clast_pprint(stdout, root, 0, cloogOptions);

    for (i=0; i<nploops; i++) {
        /* Only the innermost ones */
        if (!pluto_is_loop_innermost(ploops[i], prog)) continue;

        IF_DEBUG(printf("[pluto_mark_vector] marking loop vectorizable\n"););
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
            //loops[j]->pos_offset = ploops[i]->stmts[0]->pos_peel_offset;
            //loops[j]->neg_offset = ploops[i]->stmts[0]->neg_peel_offset;
        }

        for(j=0;j<ploops[i]->nstmts;j++)
        	ploops[i]->stmts[j]->inner_loop_vec = 1;

        free(stmtids);
        free(loops);
        free(stmts);
    }


    pluto_loops_free(ploops, nploops);
}
