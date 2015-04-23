#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"
#include "ddg.h"
#include "version.h"


#include <isl/map.h>
#include <isl/mat.h>
#include <isl/set.h>
#include <isl/flow.h>
#include <isl/union_map.h>


/*
 * 'inner_dist_loop_level' outer schedule rows will be treated as
 * parameters along with global parameters; \pi(those inner_dist_loop_level rows,
 * global params, nprocs) is what is constructed
 *
 * loop_num: parallel loop number
 */
void generate_pi(FILE *outfp, FILE *headerfp, int outer_dist_loop_level, int inner_dist_loop_level, 
        const PlutoProg *prog, Ploop *loop, int loop_num, int *dimensions_to_skip, int num_dimensions_to_skip)
{
    int i;
    Stmt *anchor_stmt = loop->stmts[0];

    assert(inner_dist_loop_level>=0);

    int j, num_inner_loops = 0;
    Ploop **inner_loops = NULL;
    if (!options->dynschedule) {
        inner_loops = pluto_get_loops_under(loop->stmts, loop->nstmts, loop->depth, 
                prog, &num_inner_loops);
    }

    int num_dims = 0;
    int dims[inner_dist_loop_level+1];
    int num_dist_dims = 0;
    int dist_dims[inner_dist_loop_level+1];
    int index_to_skip = 0;
    int is_higher_loop_dim_seq = 0;
    for (i=0; i<=inner_dist_loop_level; i++) {
        if ((index_to_skip >= num_dimensions_to_skip) || (i < dimensions_to_skip[index_to_skip])) {
            dims[num_dims++] = i;
            if (i>=outer_dist_loop_level) {
                if (options->dynschedule) {
                    // since the dimensions-to-skip have not been added to the anchor statement,
                    // is_loop should be checked on the original dimension.
                    // so get the offset to the original dimension
                    if (pluto_is_hyperplane_loop(anchor_stmt, i-index_to_skip)) {
                        // loop is distributed, even if it is not parallel
                        // since the dynamic scheduler will automatically handle
                        // the inter-task dependences
                        dist_dims[num_dist_dims++] = i;
                    }
                }
                else {
                    if (!is_higher_loop_dim_seq) {
                        if (pluto_is_hyperplane_loop(anchor_stmt, i)) {
                            // check if the loop is parallel
                            for (j=0; j<num_inner_loops; j++) {
                                if (inner_loops[j]->depth == i) {
                                    if (!pluto_loop_is_parallel(prog, inner_loops[j])) {
                                        break;
                                    }
                                }
                            }
                            if (j==num_inner_loops) {
                                dist_dims[num_dist_dims++] = i;
                            }
                            else {
                                is_higher_loop_dim_seq = 1; // should not distribute this loop, and loops inner to it
                            }
                        }
                    }
                }
            }
        } else {
            assert(i == dimensions_to_skip[index_to_skip]);
            index_to_skip++;
        }
    }
    assert(num_dims > 0);
    assert(num_dist_dims > 0);

    if (!options->dynschedule) {
        pluto_loops_free(inner_loops,num_inner_loops);
    }

    int index_offset = 1;
    if (options->dynschedule && !options->distmem) { 
        // for the scalar dimension added to separate out statements in init_tasks
        index_offset++;
    }


    fprintf(outfp, "int pi_%d(", loop_num);
    fprintf(headerfp, "int pi_%d(", loop_num);
    fprintf(outfp, " int t%d", dims[0]+index_offset);
    fprintf(headerfp, " int t%d", dims[0]+index_offset);
    for (i=1; i<num_dims; i++) {
        fprintf(outfp, ", int t%d", dims[i]+index_offset);
        fprintf(headerfp, ", int t%d", dims[i]+index_offset);
    }
    for (i=0; i<prog->npar; i++) {
        fprintf(outfp, ", int %s", prog->params[i]);
        fprintf(headerfp, ", int %s", prog->params[i]);
    }
    fprintf(outfp, ", int nprocs)\n{\n");
    fprintf(headerfp, ", int nprocs);\n");

    for (i=0; i<num_dist_dims; i++) {
        int dim = dist_dims[i] + index_offset;
        fprintf(outfp, "\
                int __lb%d, __ub%d, __p%d;\n\
                __lb%d = %s%dt%d;\n\
                __ub%d = %s%dt%d;\n", 
                i, i, i,
                i, "_LB_REPLACE_ME_DISTLOOG", loop_num, dim,
                i, "_UB_REPLACE_ME_DISTLOOG", loop_num, dim);
    }
    switch (num_dist_dims) {
        case 1:
            fprintf(outfp, "\
#ifdef __AUTO_COMPUTE_PI\n\
                    __p0 = polyrt_one_dim_pi(t%d, __lb0, __ub0, nprocs);\n\
                    return pi_mappings_%d[__p0];\n\
#endif\n", dist_dims[0]+index_offset, loop_num);
            fprintf(outfp, "\
#ifdef __USE_BLOCK_CYCLIC\n\
                    if (__is_block_cyclic[%d])\n \
                    return polyrt_one_dim_pi_block_cyclic(t%d, __lb0, __ub0, nprocs, __BLOCK_CYCLIC_BLOCK_SIZE);\n\
                    else\n \
#endif\n\
                    return polyrt_one_dim_pi(t%d, __lb0, __ub0, nprocs);\n",
                    loop_num, dist_dims[0]+index_offset,
                    dist_dims[0]+index_offset);
            break;
        case 2:
			fprintf(outfp, "\
#ifdef __AUTO_COMPUTE_PI\n\
                    __p0 = polyrt_one_dim_pi(t%d, __lb0, __ub0, nprocs);\n\
                    __p1 = polyrt_one_dim_pi(t%d, __lb1, __ub1, nprocs);\n\
                    return pi_mappings_%d[__p0][__p1];\n\
#endif\n", dist_dims[0]+index_offset, dist_dims[1]+index_offset, loop_num);

            fprintf(outfp, "\
#ifdef __USE_BLOCK_CYCLIC\n\
                    if (__is_block_cyclic[%d])\n \
                    return polyrt_two_dim_pi_block_cyclic(t%d, __lb0, __ub0, t%d, __lb1, __ub1, nprocs, __BLOCK_CYCLIC_BLOCK_SIZE);\n\
                    else\n \
#endif\n\
                    return polyrt_two_dim_pi(%d, t%d, __lb0, __ub0, t%d, __lb1, __ub1, nprocs);\n",
                    loop_num, dist_dims[0]+index_offset, dist_dims[1]+index_offset,
                    loop_num, dist_dims[0]+index_offset, dist_dims[1]+index_offset);
            break;
        case 3:
            fprintf(outfp, "\
                    return polyrt_three_dim_pi(t%d, __lb0, __ub0, t%d, __lb1, __ub1, t%d, __lb2, __ub2, nprocs);\n", 
                    dist_dims[0]+index_offset, dist_dims[1]+index_offset, dist_dims[2]+index_offset);
            break;
        default:
            assert("more than 3-dimensional distribution not supported");
    }
    fprintf(outfp, "}\n");

    if (options->distmem && options->dynschedule) {
        fprintf(outfp, "int pi_threads_%d(", loop_num);
        fprintf(headerfp, "int pi_threads_%d(", loop_num);
        fprintf(outfp, " int t%d", dims[0]+index_offset);
        fprintf(headerfp, " int t%d", dims[0]+index_offset);
        for (i=1; i<num_dims; i++) {
            fprintf(outfp, ", int t%d", dims[i]+index_offset);
            fprintf(headerfp, ", int t%d", dims[i]+index_offset);
        }
        for (i=0; i<prog->npar; i++) {
            fprintf(outfp, ", int %s", prog->params[i]);
            fprintf(headerfp, ", int %s", prog->params[i]);
        }
        fprintf(outfp, ", int my_rank, int nprocs, int num_threads)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs, int num_threads);\n");
        fprintf(outfp, "int lbp, ubp;\n");
        for (i=0; i<num_dist_dims; i++) {
            int dim = dist_dims[i] + index_offset;
            fprintf(outfp, "\
                    polyrt_multi_dim_loop_dist(%s%dt%d, %s%dt%d, nprocs, my_rank, %d, %d, &lbp, &ubp);\n",
                    "_LB_REPLACE_ME_DISTLOOG", loop_num, dim,
                    "_UB_REPLACE_ME_DISTLOOG", loop_num, dim,
                    num_dist_dims, i);
            fprintf(outfp, "\
                    int __lb%d, __ub%d;\n\
                    __lb%d = lbp;\n\
                    __ub%d = ubp;\n", 
                    i, i, i, i);
        }
        switch (num_dist_dims) {
            case 1:
                fprintf(outfp, "\
                        return polyrt_one_dim_pi(t%d, __lb0, __ub0, num_threads);\n", 
                        dist_dims[0]+index_offset);
                break;
            case 2:
                fprintf(outfp, "\
                        return polyrt_two_dim_pi(%d, t%d, __lb0, __ub0, t%d, __lb1, __ub1, num_threads);\n", 
                        loop_num, dist_dims[0]+index_offset, dist_dims[1]+index_offset);
                break;
            case 3:
                fprintf(outfp, "\
                        return polyrt_three_dim_pi(t%d, __lb0, __ub0, t%d, __lb1, __ub1, t%d, __lb2, __ub2, num_threads);\n", 
                        dist_dims[0]+index_offset, dist_dims[1]+index_offset, dist_dims[2]+index_offset);
                break;
            default:
                assert("more than 3-dimensional distribution not supported");
        }
        fprintf(outfp, "}\n");
    }

    FILE *pidefsfp = fopen("pi_defs.h", "a");
    for (i=0; i<num_dist_dims; i++) {
        fprintf(pidefsfp, "#define __DIM%dt%d %d\n", loop_num, dist_dims[i]+index_offset, i);
    }
    fprintf(pidefsfp, "#define __NUM_DIST_DIMS%d %d\n", loop_num, num_dist_dims);
    fclose(pidefsfp);
}

void pluto_gen_sigma_cloog_code(PlutoProg *sigma, FILE *cloogfp, FILE *outfp, int partition_id) {
    int i;
    for(i = 0; i < sigma->nstmts; i++){
        sigma->stmts[i]->id += partition_id+1;
    }

    pluto_gen_cloog_file(cloogfp, sigma);

    for(i = 0; i < sigma->nstmts; i++){
        sigma->stmts[i]->id -= partition_id+1;
    }
    rewind(cloogfp);

    generate_declarations(sigma, outfp);
    pluto_gen_cloog_code(sigma, -1, -1, cloogfp, outfp);
}

/* Create constraints equating first copy_level dimensions to copy_level
 * dimensions after the first src_nrows variables
 */
static PlutoConstraints *get_context_equality(int copy_level, int src_nrows, int ncols)
{
    int i;
    PlutoConstraints *cst = pluto_constraints_alloc(copy_level, ncols);
    for (i=0; i<copy_level; i++) {
        pluto_constraints_add_equality(cst);
        cst->val[cst->nrows-1][i] = 1;
        cst->val[cst->nrows-1][src_nrows+i] = -1;
    }
    return cst;
}

/*
 * generate_sigma_dep_split: Generates receivers for a given tile for a particular
 * stmt_access_pair
 * 'copy_level': outer schedule rows of all loops
 * loop_num: source loop
 * 'copy_level[loop_num]' outer schedule rows of source loop will be treated as
 * parameters along with global parameters
 */
void generate_sigma_dep_split(struct stmt_access_pair **wacc_stmts,
        int naccs, int *copy_level, PlutoProg *prog, PlutoConstraintsList *partitions,
        int loop_num, int *pi_mappings, FILE *outfp, FILE *headerfp)
{
    int i, j, partition_id;
    char *acc_name;
    int src_copy_level = copy_level[loop_num];
    PlutoConstraintsList *prev, *curr, *partitionj;

    assert(naccs != 0);

    acc_name = wacc_stmts[0]->acc->name;

    int num_partitions = 0;
    curr = partitions;
    while (curr != NULL) {
        num_partitions++;
        curr = curr->next;
    }
    PlutoConstraints *tdpoly_list_of_partition[num_partitions][prog->ndeps];
    int dep_loop_num_list_of_partition[num_partitions][prog->ndeps];
    int num_dest_loops_of_partition[num_partitions];
    int broadcast_of_partition[num_partitions];

    for (i=0; i<num_partitions; i++) broadcast_of_partition[i] = 0;

    num_partitions = 0;
    prev = NULL;
    curr = partitions;
    while (curr != NULL) {

        PlutoConstraints *tdpoly_list[prog->ndeps];
        int dep_loop_num_list[prog->ndeps];
        int num_dest_loops = 0;
        int broadcast = 0;

        // union of dependences
        PlutoDepList *curr_dep = curr->deps;
        while(curr_dep != NULL) {
            Dep *dep = curr_dep->dep;

            /* Only RAW deps matter */
            assert(dep->type == OSL_DEPENDENCE_RAW);

            Stmt *src = prog->stmts[dep->src];
            Stmt *dest = prog->stmts[dep->dest];

            int dep_loop_num = pi_mappings[dest->id];
            if (dep_loop_num == -1) { 
                // destination statement will be executed by all processors
                broadcast = 1;
            }

            int dest_copy_level = copy_level[dep_loop_num];
            int total_copy_level = src_copy_level + dest_copy_level;

            //tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
            PlutoConstraints *tdpoly = pluto_constraints_dup(dep->src_unique_dpolytope);
            //tdpoly = pluto_constraints_dup(dep->dpolytope);

            assert(!pluto_constraints_is_empty(tdpoly));

            assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

            //IF_DEBUG(printf("tile space dpoly before project_out\n"););
            IF_DEBUG(pluto_constraints_print(stdout, tdpoly););
            pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
            pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
            assert(tdpoly->ncols == total_copy_level + prog->npar+1);

            IF_DEBUG(print_polylib_visual_sets("Sigma", tdpoly));

            for (j=0; j<num_dest_loops; j++) {
                if (dep_loop_num_list[j] == dep_loop_num) {
                    pluto_constraints_unionize(tdpoly_list[j], tdpoly);
                    pluto_constraints_free(tdpoly);
                    break;
                }
            }
            if (j == num_dest_loops) {
                tdpoly_list[num_dest_loops] = tdpoly;
                dep_loop_num_list[num_dest_loops] = dep_loop_num;
                num_dest_loops++;
            }

            curr_dep = curr_dep->next;
        }
        
        for (i=num_dest_loops-1; i>=1; i--) {
            for (j=0; j<i; j++) {
                if (dep_loop_num_list[j] > dep_loop_num_list[j+1]) {
                    int temp = dep_loop_num_list[j];
                    dep_loop_num_list[j] = dep_loop_num_list[j+1];
                    dep_loop_num_list[j+1] = temp;
                    PlutoConstraints *temp_tdpoly = tdpoly_list[j];
                    tdpoly_list[j] = tdpoly_list[j+1];
                    tdpoly_list[j+1] = temp_tdpoly;
                }
            }
        }

        partitionj = partitions;
        for (j=0; j<num_partitions; j++) {
            if (num_dest_loops_of_partition[j] == num_dest_loops) {
                for (i=0; i<num_dest_loops; i++) {
                    if (dep_loop_num_list_of_partition[j][i] != dep_loop_num_list[i])
                        break;
                    int dep_loop_num = dep_loop_num_list[i];
                    int dest_copy_level = copy_level[dep_loop_num];
                    int total_copy_level = src_copy_level + dest_copy_level;
                    PlutoConstraints *tdpoly1, *tdpoly2;
                    tdpoly1 = pluto_constraints_dup(tdpoly_list_of_partition[j][i]);
                    tdpoly2 = pluto_constraints_dup(tdpoly_list[i]);
                    // do not consider constraints on parameters for equality
                    pluto_constraints_project_out(tdpoly1, total_copy_level, prog->npar);
                    pluto_constraints_project_out(tdpoly2, total_copy_level, prog->npar);
                    int constraints_are_equal = pluto_constraints_are_equal(tdpoly1, tdpoly2);
                    pluto_constraints_free(tdpoly1);
                    pluto_constraints_free(tdpoly2);
                    if (!constraints_are_equal)
                        break;
                }
                if (i==num_dest_loops) { // constraints of receiving tiles are equal for all dest loops
                    if (broadcast) broadcast_of_partition[j] = 1;
                    for (i=0; i<num_dest_loops; i++) { // because parameters are not considered for equality
                        pluto_constraints_unionize(tdpoly_list_of_partition[j][i], tdpoly_list[i]);
                    }
                    pluto_constraints_unionize(partitionj->constraints, curr->constraints);
                    PlutoDepList *dep_list = curr->deps; // !!!roshan may not need to copy/fuse deps
                    while (dep_list != NULL) {
                        pluto_deps_list_append(partitionj->deps, dep_list->dep);
                        dep_list = dep_list->next;
                    }

                    for (i=0; i<num_dest_loops; i++) pluto_constraints_free(tdpoly_list[i]);

                    prev->next = curr->next;
                    curr->next = NULL;
                    pluto_constraints_list_free(curr);
                    curr = prev->next;
                    break;
                }
            }
            partitionj = partitionj->next;
        }
        if (j == num_partitions) {
            assert(partitionj == curr);
            if (broadcast) broadcast_of_partition[num_partitions] = 1;
            num_dest_loops_of_partition[num_partitions] = num_dest_loops;
            for (i=0; i<num_dest_loops; i++) {
                dep_loop_num_list_of_partition[num_partitions][i] = dep_loop_num_list[i];
                tdpoly_list_of_partition[num_partitions][i] = tdpoly_list[i];
            }
            num_partitions++;
            prev = curr;
            curr = curr->next;
        }

    }

    curr = partitions;
    partition_id = 0;
    while (curr != NULL) {

        PlutoProg *sigma = pluto_prog_alloc();
        PlutoProg *is_receiver = pluto_prog_alloc();
        PlutoProg *sigma_check = NULL;
        if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) sigma_check = pluto_prog_alloc();

        // printf("Compute sigma\n");

        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(sigma, param, sigma->npar);
            pluto_prog_add_param(is_receiver, param, is_receiver->npar);
        }

        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(sigma, prog->params[i], sigma->npar);
            pluto_prog_add_param(is_receiver, prog->params[i], is_receiver->npar);
        }

#if 0
        for (i=0;i<src_copy_level; i++) {
            pluto_prog_add_hyperplane(sigma,0,H_LOOP);
            pluto_prog_add_hyperplane(is_receiver,0,H_LOOP);
        }
#endif

        if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) {
            for (i=0; i<src_copy_level; i++) {
                char param[6];
                sprintf(param, "ts%d",i+1);
                pluto_prog_add_param(sigma_check, param, sigma_check->npar);
            }

            for (i=0; i<prog->npar; i++) {
                pluto_prog_add_param(sigma_check, prog->params[i], sigma_check->npar);
            }

#if 0
            for (i=0;i<src_copy_level; i++) {
                pluto_prog_add_hyperplane(sigma_check,0,H_LOOP);
            }
#endif
        }

        if (!broadcast_of_partition[partition_id]) {
            for (i=0; i<num_dest_loops_of_partition[partition_id]; i++) {
                PlutoConstraints *tdpoly = tdpoly_list_of_partition[partition_id][i];
                int dep_loop_num = dep_loop_num_list_of_partition[partition_id][i];
                int dest_copy_level = copy_level[dep_loop_num];
                int total_copy_level = src_copy_level + dest_copy_level;

                if (src_copy_level == dest_copy_level) {
                    /* Interchange source and dest tile space iterators */
                    for (j=0; j<src_copy_level; j++)    {
                        pluto_constraints_interchange_cols(tdpoly, j, j+src_copy_level);
                    }
                }else{
                    // move dest_copy_level loop iterators to the beginning/top
                    for (j=0; j<dest_copy_level; j++) {
                        pluto_constraints_add_dim(tdpoly, 0, NULL);
                    }
                    for (j=0; j<dest_copy_level; j++) {
                        pluto_constraints_interchange_cols(tdpoly, j, j+total_copy_level);
                    }
                    for (j=0; j<dest_copy_level; j++) {
                        pluto_constraints_remove_dim(tdpoly, total_copy_level);
                    }
                }

                IF_MORE_DEBUG(printf("Tile space dep poly (dest, src, npar, 1)\n"););
                //IF_MORE_DEBUG(pluto_constraints_print(stdout, tdpoly););
                IF_MORE_DEBUG(print_polylib_visual_sets("tile_poly", tdpoly));

                char **iters = malloc(dest_copy_level*sizeof(char *));
                char *indices = malloc(512);
                strcpy(indices, "");

                for (j=0; j<dest_copy_level+prog->npar; j++) {
                    if (j>=1) sprintf(indices+strlen(indices), ", ");
                    if (j<=dest_copy_level-1) {
                        iters[j] = malloc(5);
                        sprintf(iters[j], "d%d", j+1);
                        sprintf(indices+strlen(indices), "%s", iters[j]);
                    }else sprintf(indices+strlen(indices), "%s", prog->params[j-dest_copy_level]);
                }

                PlutoMatrix *trans = pluto_matrix_identity(dest_copy_level);
                for (j=0; j<src_copy_level+prog->npar+1; j++) {
                    pluto_matrix_add_col(trans, trans->ncols);
                }

                char *sigma_text;
				sigma_text = malloc(strlen("__p = pi_(,nprocs); if (__p != my_rank) { receiver_list[__p] = 1; }")
						+ 5 + strlen(indices)+1);
				sprintf(sigma_text, "__p = pi_%d(%s,nprocs); if (__p != my_rank) { receiver_list[__p] = 1; }",
						dep_loop_num, indices);
                char *is_receiver_text;
                is_receiver_text = malloc(strlen("if (pi_(, nprocs) == my_rank) return 1;")
                        +5+strlen(indices)+1);
                sprintf(is_receiver_text, "if (pi_%d(%s, nprocs) == my_rank) return 1;",
                        dep_loop_num, indices);
                char *sigma_check_text = NULL;
                if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) {
                    sigma_check_text = malloc(strlen("recv_proc = pi_(, nprocs); if (recv_proc != my_rank) { if (receiver_list[recv_proc] > __FOP_UNICAST_RECV_LIMIT) return 0; else receiver_list[recv_proc]++; }")
                            +5+strlen(indices)+1);
                    sprintf(sigma_check_text,"recv_proc = pi_%d(%s, nprocs); if (recv_proc != my_rank) { if (receiver_list[recv_proc] > __FOP_UNICAST_RECV_LIMIT) return 0; else receiver_list[recv_proc]++; }",
                            dep_loop_num, indices);
                }

                pluto_add_stmt(sigma,tdpoly,trans,iters,sigma_text, IN_FUNCTION);
                pluto_add_stmt(is_receiver,tdpoly,trans,iters,is_receiver_text, IN_FUNCTION);
                if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) pluto_add_stmt(sigma_check,tdpoly,trans,iters,sigma_check_text, IN_FUNCTION);

                // pluto_stmt_print(stdout, sigma->stmts[0]);
                for (j=0; j<dest_copy_level; j++) {
                    free(iters[j]);
                }
                free(iters);

                free(sigma_text);
                free(is_receiver_text);
                if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) free(sigma_check_text);
                free(indices);
                pluto_matrix_free(trans);
                pluto_constraints_free(tdpoly);
            }

            //pluto_constraints_print(stdout, udcst);

            pluto_pad_stmt_transformations(sigma);
            pluto_pad_stmt_transformations(is_receiver);
            if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) pluto_pad_stmt_transformations(sigma_check);
            // pluto_prog_print(sigma);
            if (sigma->nstmts >= 1) {
                assert(sigma->stmts[0]->trans->nrows == sigma->num_hyperplanes);
            }

            pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts, 0, 0);
            pluto_separate_stmts(is_receiver, is_receiver->stmts, is_receiver->nstmts, 0, 0);
            if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) pluto_separate_stmts(sigma_check, sigma_check->stmts, sigma_check->nstmts, 0, 0);
        }

        FILE *cloogfp=NULL;

        fprintf(outfp, "void sigma_%s_%d_%d(", acc_name, partition_id+1, loop_num);
        fprintf(headerfp, "void sigma_%s_%d_%d(", acc_name, partition_id+1, loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs, int *receiver_list)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs, int *receiver_list);\n");

        fprintf(outfp, "int __p;\n");

        if (broadcast_of_partition[partition_id]) {
			fprintf(outfp, "\nfor (__p=0; __p<nprocs; __p++) if (__p != my_rank) receiver_list[__p] = 1;\n");
        }
        else {
            cloogfp = fopen("sigma_fop.cloog", "w+");
            assert(cloogfp != NULL);
            pluto_gen_sigma_cloog_code(sigma, cloogfp, outfp, partition_id);
            fclose(cloogfp);

            for (i=0; i<sigma->nstmts; i++) {
                fprintf(outfp, "#undef S%d\n", i+1);
            }
        }
        fprintf(outfp, "}\n\n");

        pluto_prog_free(sigma);

        fprintf(outfp, "int is_receiver_%s_%d_%d(", acc_name, partition_id+1, loop_num);
        fprintf(headerfp, "int is_receiver_%s_%d_%d(", acc_name, partition_id+1, loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs);\n");

        if (broadcast_of_partition[partition_id]) {
            fprintf(outfp, "\nreturn 1;\n");
        }
        else {
            cloogfp = fopen("is_receiver_fop.cloog", "w+");
            assert(cloogfp != NULL);
            pluto_gen_sigma_cloog_code(is_receiver, cloogfp, outfp, partition_id);
            fclose(cloogfp);

            for (i=0; i<is_receiver->nstmts; i++) {
                fprintf(outfp, "#undef S%d\n", i+1);
            }
            fprintf(outfp, "\nreturn 0;\n");
        }
        fprintf(outfp, "}\n\n");

        pluto_prog_free(is_receiver);

        if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) {
            fprintf(outfp, "int sigma_check_%s_%d_%d(", acc_name, partition_id+1, loop_num);
            fprintf(headerfp, "int sigma_check_%s_%d_%d(", acc_name, partition_id+1, loop_num);
            for (i=0; i<src_copy_level+prog->npar; i++)    {
                if (i!=0) {
                    fprintf(outfp, ", ");
                    fprintf(headerfp, ", ");
                }
                if (i<=src_copy_level-1) {
                    fprintf(outfp, "int ts%d", i+1);
                    fprintf(headerfp, "int ts%d", i+1);
                }
                else {
                    fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                    fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
                }
            }
            fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
            fprintf(headerfp, ", int my_rank, int nprocs);\n");

            fprintf(outfp, "int i, recv_proc;\n");
            fprintf(outfp, "\nint receiver_list[nprocs];\n for (i=0; i<nprocs; i++) receiver_list[i]=0;\n");

            cloogfp = fopen("sigma_check_fop.cloog", "w+");
            assert(cloogfp != NULL);
            pluto_gen_sigma_cloog_code(sigma_check, cloogfp, outfp, partition_id);
            fclose(cloogfp);

            for (i=0; i<sigma_check->nstmts; i++) {
                fprintf(outfp, "#undef S%d\n", i+1);
            }
            fprintf(outfp, "\nreturn 1;\n");
            fprintf(outfp, "}\n\n");

            pluto_prog_free(sigma_check);
        }

        curr = curr->next;
        partition_id++;
    }
    assert(partition_id == num_partitions);

}

void generate_sigma_common(struct stmt_access_pair **wacc_stmts,
        int naccs, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings, int isRemoteOnly,
		int *broadcast,
		int *num_dest_loops,
		PlutoConstraints *tdpoly_list[prog->ndeps],
		int dep_loop_num_list[prog->ndeps],
		PlutoMatrix *trans[prog->ndeps],
		char **iters[prog->ndeps],
		char *indices[prog->ndeps])
{
    int i, j;
    int src_copy_level = copy_level[loop_num];

    *num_dest_loops = 0;

    // nocommopt => no precise determination of receivers => broadcast
    *broadcast = !(options->commopt_fop || options->commopt_foifi || options->commopt);

    if (!*broadcast) {
        for (i=0; i<prog->ndeps; i++)   {
            Dep *dep = prog->deps[i];
            if (isRemoteOnly) {
                /* Only RAW deps matter */
                if (dep->type != OSL_DEPENDENCE_RAW) continue;
            }
            else {
                assert((dep->type == OSL_DEPENDENCE_RAW) || (dep->type == OSL_DEPENDENCE_WAW) || (dep->type == OSL_DEPENDENCE_WAR));
            }

            Stmt *src = prog->stmts[dep->src];
            Stmt *dest = prog->stmts[dep->dest];

            /* If the dependence doesn't originate from this access */
            for (j=0; j<naccs; j++) {
                if (dep->src_acc == wacc_stmts[j]->acc) break;
            }
            if (j==naccs)   continue;

            /* Add a statement for this dependence */
            int dep_loop_num = pi_mappings[dest->id];
            if (dep_loop_num == -1) {
                // destination statement will be executed by all processors
                *broadcast = 1;
                break;
            }

            int dest_copy_level = copy_level[dep_loop_num];
            int total_copy_level = src_copy_level + dest_copy_level;

            PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

            assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

            //IF_DEBUG(printf("tile space dpoly before project_out\n"););
            //IF_DEBUG(pluto_constraints_print(stdout, tdpoly););
            pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
            pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
            assert(tdpoly->ncols == total_copy_level + prog->npar+1);

            if (src_copy_level == dest_copy_level) {
                PlutoConstraints *param_eq = get_context_equality(src_copy_level, src_copy_level, tdpoly->ncols);
                pluto_constraints_subtract(tdpoly, param_eq);
            }

            for (j=0; j<*num_dest_loops; j++) {
                if (dep_loop_num_list[j] == dep_loop_num) {
                    pluto_constraints_unionize(tdpoly_list[j], tdpoly);
                    pluto_constraints_free(tdpoly);
                    break;
                }
            }
            if (j == *num_dest_loops) {
                tdpoly_list[*num_dest_loops] = tdpoly;
                dep_loop_num_list[*num_dest_loops] = dep_loop_num;
                (*num_dest_loops)++;
            }
        }
    }

    if (!*broadcast) {
        for (i=0; i<*num_dest_loops; i++) {
            PlutoConstraints *tdpoly = tdpoly_list[i];
            int dep_loop_num = dep_loop_num_list[i];
            int dest_copy_level = copy_level[dep_loop_num];
            int total_copy_level = src_copy_level + dest_copy_level;

            if (src_copy_level == dest_copy_level) {
                /* Interchange source and dest tile space iterators */
                for (j=0; j<src_copy_level; j++)    {
                    pluto_constraints_interchange_cols(tdpoly, j, j+src_copy_level);
                }
            }else{
                // move dest_copy_level loop iterators to the beginning/top
                for (j=0; j<dest_copy_level; j++) {
                    pluto_constraints_add_dim(tdpoly, 0, NULL);
                }
                for (j=0; j<dest_copy_level; j++) {
                    pluto_constraints_interchange_cols(tdpoly, j, j+total_copy_level);
                }
                for (j=0; j<dest_copy_level; j++) {
                    pluto_constraints_remove_dim(tdpoly, total_copy_level);
                }
            }

            IF_MORE_DEBUG(printf("Tile space dep poly (dest, src, npar, 1)\n"););
            IF_MORE_DEBUG(pluto_constraints_print(stdout, tdpoly););

            iters[i] = malloc(dest_copy_level*sizeof(char *));
            indices[i] = malloc(512);
            strcpy(indices[i], "");

            for (j=0; j<dest_copy_level+prog->npar; j++) {
                if (j>=1) sprintf(indices[i]+strlen(indices[i]), ", ");
                if (j<=dest_copy_level-1) { 
                    iters[i][j] = malloc(5);
                    sprintf(iters[i][j], "d%d", j+1);
                    sprintf(indices[i]+strlen(indices[i]), "%s", iters[i][j]);
                }else sprintf(indices[i]+strlen(indices[i]), "%s", prog->params[j-dest_copy_level]);
            }

            trans[i] = pluto_matrix_identity(dest_copy_level);
            for (j=0; j<src_copy_level+prog->npar+1; j++) {
                pluto_matrix_add_col(trans[i], trans[i]->ncols);
            }
        }
    }
}

void generate_tau_common(struct stmt_access_pair **racc_stmts,
        int naccs, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings,
		int *num_src_loops,
		PlutoConstraints *tdpoly_list[prog->ndeps],
		int src_loop_num_list[prog->ndeps],
		PlutoMatrix *trans[prog->ndeps],
		char **iters[prog->ndeps],
		char *indices[prog->ndeps])
{
    int i, j;
    int dest_copy_level = copy_level[loop_num];

    *num_src_loops = 0;

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* If the dependence isn't incident on these accesses */
        for (j=0; j<naccs; j++) {
            if (dep->dest_acc == racc_stmts[j]->acc) break;
        }
        if (j==naccs) continue;

        /* Add a statement for this dependence */
        int src_loop_num = pi_mappings[src->id];
        assert (src_loop_num != -1);

        int src_copy_level = copy_level[src_loop_num];
        int total_copy_level = src_copy_level + dest_copy_level;

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //IF_DEBUG(printf("tile space dpoly before project_out\n"););
        //IF_DEBUG(pluto_constraints_print(stdout, tdpoly););
        pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
        pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
        assert(tdpoly->ncols == total_copy_level + prog->npar+1);

        if (src_copy_level == dest_copy_level) {
            PlutoConstraints *param_eq = get_context_equality(src_copy_level, src_copy_level, tdpoly->ncols);
            pluto_constraints_subtract(tdpoly, param_eq);
        }

        for (j=0; j<*num_src_loops; j++) {
            if (src_loop_num_list[j] == src_loop_num) {
                pluto_constraints_unionize(tdpoly_list[j], tdpoly);
                pluto_constraints_free(tdpoly);
                break;
            }
        }
        if (j == *num_src_loops) {
            tdpoly_list[*num_src_loops] = tdpoly;
            src_loop_num_list[*num_src_loops] = src_loop_num;
            (*num_src_loops)++;
        }
    }

    for (i=0; i<*num_src_loops; i++) {
        PlutoConstraints *tdpoly = tdpoly_list[i];
        int src_loop_num = src_loop_num_list[i];
        int src_copy_level = copy_level[src_loop_num];

        IF_MORE_DEBUG(printf("Tile space dep poly (src, dest, npar, 1)\n"););
        IF_MORE_DEBUG(pluto_constraints_print(stdout, tdpoly););

        iters[i] = malloc(src_copy_level*sizeof(char *));
        indices[i] = malloc(512);
        strcpy(indices[i], "");

        for (j=0; j<src_copy_level+prog->npar; j++) {
            if (j>=1) sprintf(indices[i]+strlen(indices[i]), ", ");
            if (j<=src_copy_level-1) {
                iters[i][j] = malloc(5);
                sprintf(iters[i][j], "d%d", j+1);
                sprintf(indices[i]+strlen(indices[i]), "%s", iters[i][j]);
            }else sprintf(indices[i]+strlen(indices[i]), "%s", prog->params[j-src_copy_level]);
        }

        trans[i] = pluto_matrix_identity(src_copy_level);
        for (j=0; j<dest_copy_level+prog->npar+1; j++) {
            pluto_matrix_add_col(trans[i], trans[i]->ncols);
        }
    }
}

void generate_outgoing_common(Stmt **loop_stmts,
        int nstmts, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings, int isRemoteOnly,
		int *num_dest_loops,
		PlutoConstraints *tdpoly_list[prog->ndeps],
		int dep_loop_num_list[prog->ndeps],
		PlutoMatrix *trans[prog->ndeps],
		char **iters[prog->ndeps],
		char *indices[prog->ndeps],
		char *task_indices[prog->ndeps])
{
    int i, j, k;
    int src_copy_level = copy_level[loop_num];

    *num_dest_loops = 0;

    int ndeps;
    Dep **deps;
    if (options->dyn_trans_deps_tasks 
            && !isRemoteOnly && options->distmem && options->lastwriter) {
        ndeps = prog->ntransdeps;
        deps = prog->transdeps;
    }
    else {
        ndeps = prog->ndeps;
        deps = prog->deps;
    }

    Stmt *anchor_stmt_list[ndeps];

    for (i=0; i<ndeps; i++)   {
        Dep *dep = deps[i];
        if (isRemoteOnly) {
            /* Only RAW deps matter */
            if (dep->type != OSL_DEPENDENCE_RAW) continue;
        }
        else {
            assert((dep->type == OSL_DEPENDENCE_RAW) || (dep->type == OSL_DEPENDENCE_WAW) || (dep->type == OSL_DEPENDENCE_WAR));
        }

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* If the dependence doesn't originate from these statements */
        for (j=0; j<nstmts; j++) {
            if (src->id == loop_stmts[j]->id) break;
        }
        if (j==nstmts) continue;

        /* Add a statement for this dependence */
        int dep_loop_num = pi_mappings[dest->id];
        assert (dep_loop_num != -1);

        int dest_copy_level = copy_level[dep_loop_num];
        int total_copy_level = src_copy_level + dest_copy_level;

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //IF_DEBUG(printf("tile space dpoly before project_out\n"););
        //IF_DEBUG(pluto_constraints_print(stdout, tdpoly););
        pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
        pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
        assert(tdpoly->ncols == total_copy_level + prog->npar+1);

        if (src_copy_level == dest_copy_level) {
            PlutoConstraints *param_eq = get_context_equality(src_copy_level, src_copy_level, tdpoly->ncols);
            pluto_constraints_subtract(tdpoly, param_eq);
        }

        for (j=0; j<*num_dest_loops; j++) {
            if (dep_loop_num_list[j] == dep_loop_num) {
                pluto_constraints_unionize(tdpoly_list[j], tdpoly);
                pluto_constraints_free(tdpoly);
                break;
            }
        }
        if (j == *num_dest_loops) {
            tdpoly_list[*num_dest_loops] = tdpoly;
            anchor_stmt_list[*num_dest_loops] = dest;
            dep_loop_num_list[*num_dest_loops] = dep_loop_num;
            (*num_dest_loops)++;
        }
    }

    for (i=0; i<*num_dest_loops; i++) {
        PlutoConstraints *tdpoly = tdpoly_list[i];
        Stmt *anchor_stmt = anchor_stmt_list[i];
        int dep_loop_num = dep_loop_num_list[i];
        int dest_copy_level = copy_level[dep_loop_num];
        int total_copy_level = src_copy_level + dest_copy_level;

        if (src_copy_level == dest_copy_level) {
            /* Interchange source and dest tile space iterators */
            for (j=0; j<src_copy_level; j++)    {
                pluto_constraints_interchange_cols(tdpoly, j, j+src_copy_level);
            }
        }else{
            // move dest_copy_level loop iterators to the beginning/top
            for (j=0; j<dest_copy_level; j++) {
                pluto_constraints_add_dim(tdpoly, 0, NULL);
            }
            for (j=0; j<dest_copy_level; j++) {
                pluto_constraints_interchange_cols(tdpoly, j, j+total_copy_level);
            }
            for (j=0; j<dest_copy_level; j++) {
                pluto_constraints_remove_dim(tdpoly, total_copy_level);
            }
        }

        IF_MORE_DEBUG(printf("Tile space dep poly (dest, src, npar, 1)\n"););
        IF_MORE_DEBUG(pluto_constraints_print(stdout, tdpoly););

        iters[i] = malloc(dest_copy_level*sizeof(char *));
        indices[i] = malloc(512);
        strcpy(indices[i], "");

        for (j=0; j<dest_copy_level; j++) {
            if (j>=1) sprintf(indices[i]+strlen(indices[i]), ", ");
            iters[i][j] = malloc(5);
            sprintf(iters[i][j], "d%d", j+1);
            sprintf(indices[i]+strlen(indices[i]), "%s", iters[i][j]);
        }

        int dims[dest_copy_level];
        int num_dims = 0;
        for (j=0; j<dest_copy_level; j++) {
            if (anchor_stmt->hyp_types[j] != H_SCALAR) {
                dims[num_dims++] = j;
            }
        }

        task_indices[i] = malloc(1024);
        strcpy(task_indices[i],"[");
        for (j=0; j<num_dims; j++) {
            sprintf(task_indices[i]+strlen(task_indices[i]), "(%s-lb_tasks_loop%d_dim%d)",
                    iters[i][dims[j]], dep_loop_num, j);
            for (k=j+1; k<num_dims; k++) {
                sprintf(task_indices[i]+strlen(task_indices[i]), "*max_num_tasks_loop%d_dim%d",
                        dep_loop_num, k);
            }
            strcat(task_indices[i],"+");
        }
        strcat(task_indices[i],"0]");

        trans[i] = pluto_matrix_identity(dest_copy_level);
        for (j=0; j<src_copy_level+prog->npar+1; j++) {
            pluto_matrix_add_col(trans[i], trans[i]->ncols);
        }
    }
}

void generate_incoming_common(Stmt **loop_stmts,
        int nstmts, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings, int isRemoteOnly,
		int *num_src_loops,
		PlutoConstraints *tdpoly_list[prog->ndeps],
		int src_loop_num_list[prog->ndeps],
		PlutoMatrix *trans[prog->ndeps],
		char **iters[prog->ndeps],
		char *indices[prog->ndeps])
{
    int i, j;
    int dest_copy_level = copy_level[loop_num];

    *num_src_loops = 0;

    int ndeps;
    Dep **deps;
    if (options->dyn_trans_deps_tasks 
            && !isRemoteOnly && options->distmem && options->lastwriter) {
        ndeps = prog->ntransdeps;
        deps = prog->transdeps;
    }
    else {
        ndeps = prog->ndeps;
        deps = prog->deps;
    }

    for (i=0; i<ndeps; i++)   {
        Dep *dep = deps[i];
        if (isRemoteOnly) {
            /* Only RAW deps matter */
            if (dep->type != OSL_DEPENDENCE_RAW) continue;
        }
        else {
            assert((dep->type == OSL_DEPENDENCE_RAW) || (dep->type == OSL_DEPENDENCE_WAW) || (dep->type == OSL_DEPENDENCE_WAR));
        }

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* If the dependence isn't incident on these statements */
        for (j=0; j<nstmts; j++) {
            if (dest->id == loop_stmts[j]->id) break;
        }
        if (j==nstmts) continue;

        /* Add a statement for this dependence */
        int src_loop_num = pi_mappings[src->id];
        assert (src_loop_num != -1);

        int src_copy_level = copy_level[src_loop_num];
        int total_copy_level = src_copy_level + dest_copy_level;

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //IF_DEBUG(printf("tile space dpoly before project_out\n"););
        //IF_DEBUG(pluto_constraints_print(stdout, tdpoly););
        pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
        pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
        assert(tdpoly->ncols == total_copy_level + prog->npar+1);

        if (src_copy_level == dest_copy_level) {
            PlutoConstraints *param_eq = get_context_equality(src_copy_level, src_copy_level, tdpoly->ncols);
            pluto_constraints_subtract(tdpoly, param_eq);
        }

        for (j=0; j<*num_src_loops; j++) {
            if (src_loop_num_list[j] == src_loop_num) {
                pluto_constraints_unionize(tdpoly_list[j], tdpoly);
                pluto_constraints_free(tdpoly);
                break;
            }
        }
        if (j == *num_src_loops) {
            tdpoly_list[*num_src_loops] = tdpoly;
            src_loop_num_list[*num_src_loops] = src_loop_num;
            (*num_src_loops)++;
        }
    }

    for (i=0; i<*num_src_loops; i++) {
        PlutoConstraints *tdpoly = tdpoly_list[i];
        int src_loop_num = src_loop_num_list[i];
        int src_copy_level = copy_level[src_loop_num];

        IF_MORE_DEBUG(printf("Tile space dep poly (src, dest, npar, 1)\n"););
        IF_MORE_DEBUG(pluto_constraints_print(stdout, tdpoly););

        iters[i] = malloc(src_copy_level*sizeof(char *));
        indices[i] = malloc(512);
        strcpy(indices[i], "");

        for (j=0; j<src_copy_level+prog->npar; j++) {
            if (j>=1) sprintf(indices[i]+strlen(indices[i]), ", ");
            if (j<=src_copy_level-1) {
                iters[i][j] = malloc(5);
                sprintf(iters[i][j], "d%d", j+1);
                sprintf(indices[i]+strlen(indices[i]), "%s", iters[i][j]);
            }else sprintf(indices[i]+strlen(indices[i]), "%s", prog->params[j-src_copy_level]);
        }

        trans[i] = pluto_matrix_identity(src_copy_level);
        for (j=0; j<dest_copy_level+prog->npar+1; j++) {
            pluto_matrix_add_col(trans[i], trans[i]->ncols);
        }
    }
}

void generate_outgoing(Stmt **loop_stmts,
        int nstmts, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings,
        char *tasks_loops_decl, FILE *outfp, FILE *headerfp)
{
    int i, j;
    int src_copy_level = copy_level[loop_num];

    PlutoProg *is_receiver = NULL;
    PlutoProg *count_remote_dep_tasks = NULL;
    PlutoProg *remote_update_dep_tasks = NULL;
    if (options->distmem) {
        is_receiver = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(is_receiver, param, is_receiver->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(is_receiver, prog->params[i], is_receiver->npar);
        }

        count_remote_dep_tasks = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(count_remote_dep_tasks, param, count_remote_dep_tasks->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(count_remote_dep_tasks, prog->params[i], count_remote_dep_tasks->npar);
        }

        remote_update_dep_tasks = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(remote_update_dep_tasks, param, remote_update_dep_tasks->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(remote_update_dep_tasks, prog->params[i], remote_update_dep_tasks->npar);
        }
    }

    PlutoProg *count_local_dep_tasks = NULL;
    PlutoProg *local_update_dep_tasks = NULL;
    if (options->dynschedule) {
        count_local_dep_tasks = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(count_local_dep_tasks, param, count_local_dep_tasks->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(count_local_dep_tasks, prog->params[i], count_local_dep_tasks->npar);
        }

        local_update_dep_tasks = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(local_update_dep_tasks, param, local_update_dep_tasks->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(local_update_dep_tasks, prog->params[i], local_update_dep_tasks->npar);
        }
    }

    PlutoProg *add_outgoing_edges = NULL;
    if (options->dynschedule_graph) {
        add_outgoing_edges = pluto_prog_alloc();
        for (i=0; i<src_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(add_outgoing_edges, param, add_outgoing_edges->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(add_outgoing_edges, prog->params[i], add_outgoing_edges->npar);
        }
    }

    char *params = malloc(256);
    strcpy(params, "");
    if (prog->npar>=1) {
        sprintf(params+strlen(params), "%s", prog->params[0]);
        for (i=1; i<prog->npar; i++) {
            sprintf(params+strlen(params), ",%s", prog->params[i]);
        }
    }

    if (options->distmem) {
        int num_dest_loops;
        PlutoConstraints *tdpoly_list[prog->ndeps];
        int dep_loop_num_list[prog->ndeps];
        PlutoMatrix *trans[prog->ndeps];
        char **iters[prog->ndeps];
        char *indices[prog->ndeps];
        char *task_indices[prog->ndeps];

        // for remote dependences only i.e., RAW dependences only
        generate_outgoing_common(loop_stmts, nstmts, copy_level, prog, loop_num, pi_mappings, 1,
                &num_dest_loops,
                tdpoly_list,
                dep_loop_num_list,
                trans,
                iters,
                indices,
                task_indices);

        for (i=0; i<num_dest_loops; i++) {
            PlutoConstraints *tdpoly = tdpoly_list[i];
            int dep_loop_num = dep_loop_num_list[i];
            int dest_copy_level = copy_level[dep_loop_num];

            char *is_receiver_text = NULL;
            is_receiver_text = malloc(strlen("if (pi_(,, nprocs) == my_rank) return 1;")
                    +5+strlen(indices[i])+strlen(params)+1);
            sprintf(is_receiver_text, "if (pi_%d(%s,%s, nprocs) == my_rank) return 1;",
                    dep_loop_num, indices[i], params);
            pluto_add_stmt(is_receiver,tdpoly,trans[i],iters[i],is_receiver_text, IN_FUNCTION);

            char *count_remote_dep_tasks_text = NULL;
            count_remote_dep_tasks_text = malloc(strlen("if (pi_(,, nprocs) != my_rank) count++;")
                    +5+strlen(indices[i])+strlen(params)+1);
            sprintf(count_remote_dep_tasks_text, "if (pi_%d(%s,%s, nprocs) != my_rank) count++;",
                    dep_loop_num, indices[i], params);
            pluto_add_stmt(count_remote_dep_tasks,tdpoly,trans[i],iters[i],count_remote_dep_tasks_text, IN_FUNCTION);

            char *remote_update_dep_tasks_text = malloc(2048);
            sprintf(remote_update_dep_tasks_text, "if (pi_%d(%s,%s, nprocs) == my_rank) { \
                    __omp_atomic_capture captured_firing_count = --tasks_loop%d%s; \
                    if (captured_firing_count == 0) { \
                    remote_dep_tasks = get_num_remote_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                    local_dep_tasks = get_num_local_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                    if (__is_multi_partitioned[%d] || __is_block_cyclic[%d]) affinity = -1; \
                    else affinity = pi_threads_%d(%s,%s,my_rank,nprocs,num_threads); \
                    __Task task(%d,%s,remote_dep_tasks,local_dep_tasks,affinity); \
                    pqueue.push(task); ",
                    dep_loop_num, indices[i], params, dep_loop_num, task_indices[i], 
                    dep_loop_num, indices[i], params, dep_loop_num, indices[i], params,
                    dep_loop_num, dep_loop_num, 
                    dep_loop_num, indices[i], params, dep_loop_num, indices[i]);
            sprintf(remote_update_dep_tasks_text+strlen(remote_update_dep_tasks_text),
                    "IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                    fprintf(__debug_print_fp, \"added node %%d loop %d task ",
                    dep_loop_num);
            for (j=0; j<dest_copy_level; j++) {
                sprintf(remote_update_dep_tasks_text+strlen(remote_update_dep_tasks_text),  "%%d ");
            }
            sprintf(remote_update_dep_tasks_text+strlen(remote_update_dep_tasks_text),  "\
                    remote_dep_tasks %%d local_dep_tasks %%d ready tasks %%lu\\n\", \
                    my_rank, %s, remote_dep_tasks, local_dep_tasks, pqueue.size())); \
                    IF_DEBUG_FLUSH(fflush(__debug_print_fp)); } }",
                    indices[i]);
            pluto_add_stmt(remote_update_dep_tasks,tdpoly,trans[i],iters[i],remote_update_dep_tasks_text, IN_FUNCTION);

            for (j=0; j<dest_copy_level; j++) {
                free(iters[i][j]);
            }

            free(is_receiver_text);
            free(count_remote_dep_tasks_text);
            free(remote_update_dep_tasks_text);
            free(indices[i]);
            free(task_indices[i]);
            free(iters[i]);
            pluto_matrix_free(trans[i]);
            pluto_constraints_free(tdpoly);
        }

        pluto_pad_stmt_transformations(is_receiver);
        if (is_receiver->nstmts >= 1) {
            assert(is_receiver->stmts[0]->trans->nrows == is_receiver->num_hyperplanes);
        }
        pluto_separate_stmts(is_receiver, is_receiver->stmts, is_receiver->nstmts, 0, 0);

        pluto_pad_stmt_transformations(count_remote_dep_tasks);
        if (count_remote_dep_tasks->nstmts >= 1) {
            assert(count_remote_dep_tasks->stmts[0]->trans->nrows == count_remote_dep_tasks->num_hyperplanes);
        }
        pluto_separate_stmts(count_remote_dep_tasks, count_remote_dep_tasks->stmts, count_remote_dep_tasks->nstmts, 0, 0);

        pluto_pad_stmt_transformations(remote_update_dep_tasks);
        if (remote_update_dep_tasks->nstmts >= 1) {
            assert(remote_update_dep_tasks->stmts[0]->trans->nrows == remote_update_dep_tasks->num_hyperplanes);
        }
        pluto_separate_stmts(remote_update_dep_tasks, remote_update_dep_tasks->stmts, remote_update_dep_tasks->nstmts, 0, 0);
    }

    {
        int local_num_dest_loops;
        PlutoConstraints *local_tdpoly_list[prog->ndeps];
        int local_dep_loop_num_list[prog->ndeps];
        PlutoMatrix *local_trans[prog->ndeps];
        char **local_iters[prog->ndeps];
        char *local_indices[prog->ndeps];
        char *local_task_indices[prog->ndeps];

        // for all memory dependences only i.e., RAW, WAW, WAR dependences
        generate_outgoing_common(loop_stmts, nstmts, copy_level, prog, loop_num, pi_mappings, 0,
                &local_num_dest_loops,
                local_tdpoly_list,
                local_dep_loop_num_list,
                local_trans,
                local_iters,
                local_indices,
                local_task_indices);

        for (i=0; i<local_num_dest_loops; i++) {
            PlutoConstraints *tdpoly = local_tdpoly_list[i];
            int dep_loop_num = local_dep_loop_num_list[i];
            int dest_copy_level = copy_level[dep_loop_num];

            char *count_local_dep_tasks_text = NULL;
            if (options->dynschedule) {
                if (options->distmem) {
                    count_local_dep_tasks_text = malloc(strlen("if (pi_(,, nprocs) == my_rank) count++;")
                            +5+strlen(local_indices[i])+strlen(params)+1);
                    sprintf(count_local_dep_tasks_text, "if (pi_%d(%s,%s, nprocs) == my_rank) count++;",
                            dep_loop_num, local_indices[i], params);
                }
                else {
                    count_local_dep_tasks_text = malloc(16);
                    strcpy(count_local_dep_tasks_text, "count++;");
                }
                pluto_add_stmt(count_local_dep_tasks,tdpoly,local_trans[i],local_iters[i],count_local_dep_tasks_text, IN_FUNCTION);
                free(count_local_dep_tasks_text);
            }

            char *local_update_dep_tasks_text = NULL;
            if (options->dynschedule) {
                local_update_dep_tasks_text = malloc(2048);
                if (options->distmem) {
                    sprintf(local_update_dep_tasks_text,
                            "__omp_atomic_capture captured_firing_count = --tasks_loop%d%s; \
                            if (captured_firing_count == 0) { \
                            if (pi_%d(%s,%s, nprocs) == my_rank) { \
                            remote_dep_tasks = get_num_remote_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                            local_dep_tasks = get_num_local_dep_tasks_%d(%s,%s,my_rank,nprocs); \
                            if (__is_multi_partitioned[%d] || __is_block_cyclic[%d]) affinity = -1; \
                            else affinity = pi_threads_%d(%s,%s,my_rank,nprocs,num_threads); \
                            __Task task(%d,%s,remote_dep_tasks,local_dep_tasks,affinity); \
                            pqueue.push(task); ",
                            dep_loop_num, local_task_indices[i], dep_loop_num, local_indices[i], params,
                            dep_loop_num, local_indices[i], params, dep_loop_num, local_indices[i], params,
                            dep_loop_num, dep_loop_num, 
                            dep_loop_num, local_indices[i], params, dep_loop_num, local_indices[i]);
                    sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),
                            "IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                            fprintf(__debug_print_fp, \"added node %%d loop %d task ",
                            dep_loop_num);
                    for (j=0; j<dest_copy_level; j++) {
                        sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),  "%%d ");
                    }
                    sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),  "\
                            remote_dep_tasks %%d local_dep_tasks %%d ready tasks %%lu\\n\", \
                            my_rank, %s, remote_dep_tasks, local_dep_tasks, pqueue.size())); \
                            IF_DEBUG_FLUSH(fflush(__debug_print_fp)); } }",
                            local_indices[i]);
                }
                else {
                    sprintf(local_update_dep_tasks_text,
                            "__omp_atomic_capture captured_firing_count = --tasks_loop%d%s; \
                            if (captured_firing_count == 0) { \
                            local_dep_tasks = get_num_local_dep_tasks_%d(%s,%s); \
                            affinity = pi_%d(%s,%s,num_threads); \
                            __Task task(%d,%s,local_dep_tasks,affinity); \
                            pqueue.push(task); ",
                            dep_loop_num, local_task_indices[i], 
                            dep_loop_num, local_indices[i], params,
                            dep_loop_num, local_indices[i], params, 
                            dep_loop_num, local_indices[i]);
                    sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),
                            "IF_DYNSCHEDULER_MORE_DEBUG_PRINT(\
                            fprintf(__debug_print_fp, \"added node %%d loop %d task ",
                            dep_loop_num);
                    for (j=0; j<dest_copy_level; j++) {
                        sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),  "%%d ");
                    }
                    sprintf(local_update_dep_tasks_text+strlen(local_update_dep_tasks_text),  "\
                            local_dep_tasks %%d ready tasks %%lu\\n\", \
                            my_rank, %s, local_dep_tasks, pqueue.size())); \
                            IF_DEBUG_FLUSH(fflush(__debug_print_fp)); }",
                            local_indices[i]);
                }
                pluto_add_stmt(local_update_dep_tasks,tdpoly,local_trans[i],local_iters[i],local_update_dep_tasks_text, IN_FUNCTION);
                free(local_update_dep_tasks_text);
            }

            char *add_outgoing_edges_text = NULL;
            if (options->dynschedule_graph) {
                add_outgoing_edges_text = malloc(1024);
                sprintf(add_outgoing_edges_text, 
                        "tbb::flow::make_edge(*current_task, *tasks_loop%d%s);",
                        dep_loop_num, local_task_indices[i]);
                pluto_add_stmt(add_outgoing_edges,tdpoly,local_trans[i],local_iters[i],add_outgoing_edges_text, IN_FUNCTION);
                free(add_outgoing_edges_text);
            }

            free(local_indices[i]);
            free(local_task_indices[i]);
            free(local_iters[i]);
            pluto_matrix_free(local_trans[i]);
            pluto_constraints_free(tdpoly);
        }

        if (options->dynschedule) {
            pluto_pad_stmt_transformations(count_local_dep_tasks);
            if (count_local_dep_tasks->nstmts >= 1) {
                assert(count_local_dep_tasks->stmts[0]->trans->nrows == count_local_dep_tasks->num_hyperplanes);
            }
            pluto_separate_stmts(count_local_dep_tasks, count_local_dep_tasks->stmts, count_local_dep_tasks->nstmts, 0, 0);

            pluto_pad_stmt_transformations(local_update_dep_tasks);
            if (local_update_dep_tasks->nstmts >= 1) {
                assert(local_update_dep_tasks->stmts[0]->trans->nrows == local_update_dep_tasks->num_hyperplanes);
            }
            pluto_separate_stmts(local_update_dep_tasks, local_update_dep_tasks->stmts, local_update_dep_tasks->nstmts, 0, 0);
        }

        if (options->dynschedule_graph) {
            pluto_pad_stmt_transformations(add_outgoing_edges);
            if (add_outgoing_edges->nstmts >= 1) {
                assert(add_outgoing_edges->stmts[0]->trans->nrows == add_outgoing_edges->num_hyperplanes);
            }
            pluto_separate_stmts(add_outgoing_edges, add_outgoing_edges->stmts, add_outgoing_edges->nstmts, 0, 0);
        }
    }
    free(params);

    FILE *cloogfp = NULL;

    if (options->distmem) {
        fprintf(outfp, "int is_receiver_%d(", loop_num);
        fprintf(headerfp, "int is_receiver_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs);\n");
        cloogfp = fopen("is_receiver.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, is_receiver);
        rewind(cloogfp);
        generate_declarations(is_receiver, outfp);
        pluto_gen_cloog_code(is_receiver, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<is_receiver->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "\nreturn 0;\n");
        fprintf(outfp, "}\n\n");

        fprintf(outfp, "int get_num_remote_dep_tasks_%d(", loop_num);
        fprintf(headerfp, "int get_num_remote_dep_tasks_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs);\n");
        fprintf(outfp, "\nint count = 0;\n");
        cloogfp = fopen("count_remote_dep_tasks.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, count_remote_dep_tasks);
        rewind(cloogfp);
        generate_declarations(count_remote_dep_tasks, outfp);
        pluto_gen_cloog_code(count_remote_dep_tasks, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<count_remote_dep_tasks->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "\nreturn count;\n");
        fprintf(outfp, "}\n\n");

        fprintf(outfp, "void remote_update_dep_tasks_%d(", loop_num);
        fprintf(headerfp, "void remote_update_dep_tasks_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs, int num_threads%s)\n{\n", tasks_loops_decl);
        fprintf(headerfp, ", int my_rank, int nprocs, int num_threads%s);\n", tasks_loops_decl);
        fprintf(outfp, "int remote_dep_tasks, local_dep_tasks, affinity, captured_firing_count;\n");
        cloogfp = fopen("remote_update_dep_tasks.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, remote_update_dep_tasks);
        rewind(cloogfp);
        generate_declarations(remote_update_dep_tasks, outfp);
        pluto_gen_cloog_code(remote_update_dep_tasks, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<remote_update_dep_tasks->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "}\n\n");
    }

    if (options->dynschedule) {
        fprintf(outfp, "int get_num_local_dep_tasks_%d(", loop_num);
        fprintf(headerfp, "int get_num_local_dep_tasks_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        if (options->distmem) {
            fprintf(outfp, ", int my_rank, int nprocs");
            fprintf(headerfp, ", int my_rank, int nprocs");
        }
        fprintf(outfp, ")\n{\n");
        fprintf(headerfp, ");\n");
        fprintf(outfp, "\nint count = 0;\n");
        cloogfp = fopen("count_local_dep_tasks.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, count_local_dep_tasks);
        rewind(cloogfp);
        generate_declarations(count_local_dep_tasks, outfp);
        pluto_gen_cloog_code(count_local_dep_tasks, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<count_local_dep_tasks->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "\nreturn count;\n");
        fprintf(outfp, "}\n\n");

        fprintf(outfp, "void local_update_dep_tasks_%d(", loop_num);
        fprintf(headerfp, "void local_update_dep_tasks_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        if (options->distmem) {
            fprintf(outfp, ", int my_rank, int nprocs");
            fprintf(headerfp, ", int my_rank, int nprocs");
        }
        fprintf(outfp, ", int num_threads");
        fprintf(headerfp, ", int num_threads");
        fprintf(outfp, "%s)\n{\n", tasks_loops_decl);
        fprintf(headerfp, "%s);\n", tasks_loops_decl);
        fprintf(outfp, "int remote_dep_tasks, local_dep_tasks, affinity, captured_firing_count;\n");
        cloogfp = fopen("local_update_dep_tasks.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, local_update_dep_tasks);
        rewind(cloogfp);
        generate_declarations(local_update_dep_tasks, outfp);
        pluto_gen_cloog_code(local_update_dep_tasks, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<local_update_dep_tasks->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "}\n\n");
    }

    if (options->dynschedule_graph) {
        fprintf(outfp, "void add_outgoing_edges_%d(", loop_num);
        fprintf(headerfp, "void add_outgoing_edges_%d(", loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=src_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
            }
        }
        fprintf(outfp, "%s, tbb::flow::continue_node<tbb::flow::continue_msg> *current_task)\n{\n", tasks_loops_decl);
        fprintf(headerfp, "%s, tbb::flow::continue_node<tbb::flow::continue_msg> *current_task);\n", tasks_loops_decl);
        fprintf(outfp, "\nint count = 0;\n");
        cloogfp = fopen("add_outgoing_edges.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, add_outgoing_edges);
        rewind(cloogfp);
        generate_declarations(add_outgoing_edges, outfp);
        pluto_gen_cloog_code(add_outgoing_edges, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<add_outgoing_edges->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "}\n\n");
    }

    if (options->distmem) {
        pluto_prog_free(is_receiver);
        pluto_prog_free(count_remote_dep_tasks);
        pluto_prog_free(remote_update_dep_tasks);
    }
    if (options->dynschedule) {
        pluto_prog_free(count_local_dep_tasks);
        pluto_prog_free(local_update_dep_tasks);
    }
    if (options->dynschedule_graph) {
        pluto_prog_free(add_outgoing_edges);
    }
}

void generate_incoming(Stmt **loop_stmts,
        int nstmts, int *copy_level, PlutoProg *prog, int loop_num,
        int *pi_mappings, FILE *outfp, FILE *headerfp)
{
    int i, j;
    int dest_copy_level = copy_level[loop_num];

    PlutoProg *count_remote_src_tasks = NULL;
    if (options->distmem) {
        count_remote_src_tasks = pluto_prog_alloc();
        for (i=0; i<dest_copy_level; i++) {
            char param[6];
            sprintf(param, "ts%d",i+1);
            pluto_prog_add_param(count_remote_src_tasks, param, count_remote_src_tasks->npar);
        }
        for (i=0; i<prog->npar; i++) {
            pluto_prog_add_param(count_remote_src_tasks, prog->params[i], count_remote_src_tasks->npar);
        }
    }

    PlutoProg *count_local_src_tasks = NULL;
    count_local_src_tasks = pluto_prog_alloc();
    for (i=0; i<dest_copy_level; i++) {
        char param[6];
        sprintf(param, "ts%d",i+1);
        pluto_prog_add_param(count_local_src_tasks, param, count_local_src_tasks->npar);
    }
    for (i=0; i<prog->npar; i++) {
        pluto_prog_add_param(count_local_src_tasks, prog->params[i], count_local_src_tasks->npar);
    }

    if (options->distmem) {
        int num_src_loops;
        PlutoConstraints *tdpoly_list[prog->ndeps];
        int src_loop_num_list[prog->ndeps];
        PlutoMatrix *trans[prog->ndeps];
        char **iters[prog->ndeps];
        char *indices[prog->ndeps];

        // for remote dependences only i.e., RAW dependences only
        generate_incoming_common(loop_stmts, nstmts, copy_level, prog, loop_num, pi_mappings, 1,
                &num_src_loops,
                tdpoly_list,
                src_loop_num_list,
                trans,
                iters,
                indices);

        for (i=0; i<num_src_loops; i++) {
            PlutoConstraints *tdpoly = tdpoly_list[i];
            int src_loop_num = src_loop_num_list[i];
            int src_copy_level = copy_level[src_loop_num];

            char *count_remote_src_tasks_text = NULL;
            count_remote_src_tasks_text = malloc(strlen("if (pi_(, nprocs) != my_rank) count;")
                    +5+strlen(indices[i])+1);
            sprintf(count_remote_src_tasks_text, "if (pi_%d(%s, nprocs) != my_rank) count++;",
                    src_loop_num, indices[i]);
            pluto_add_stmt(count_remote_src_tasks,tdpoly,trans[i],iters[i],count_remote_src_tasks_text, IN_FUNCTION);

            for (j=0; j<src_copy_level; j++) {
                free(iters[i][j]);
            }

            free(count_remote_src_tasks_text);
            free(indices[i]);
            free(iters[i]);
            pluto_matrix_free(trans[i]);
            pluto_constraints_free(tdpoly);
        }

        pluto_pad_stmt_transformations(count_remote_src_tasks);
        if (count_remote_src_tasks->nstmts >= 1) {
            assert(count_remote_src_tasks->stmts[0]->trans->nrows == count_remote_src_tasks->num_hyperplanes);
        }
        pluto_separate_stmts(count_remote_src_tasks, count_remote_src_tasks->stmts, count_remote_src_tasks->nstmts, 0, 0);
    }

    {
        int local_num_src_loops;
        PlutoConstraints *local_tdpoly_list[prog->ndeps];
        int local_src_loop_num_list[prog->ndeps];
        PlutoMatrix *local_trans[prog->ndeps];
        char **local_iters[prog->ndeps];
        char *local_indices[prog->ndeps];

        // for all memory dependences only i.e., RAW, WAW, WAR dependences
        generate_incoming_common(loop_stmts, nstmts, copy_level, prog, loop_num, pi_mappings, 0,
                &local_num_src_loops,
                local_tdpoly_list,
                local_src_loop_num_list,
                local_trans,
                local_iters,
                local_indices);

        for (i=0; i<local_num_src_loops; i++) {
            PlutoConstraints *tdpoly = local_tdpoly_list[i];
            int src_loop_num = local_src_loop_num_list[i];
            int src_copy_level = copy_level[src_loop_num];

            char *count_local_src_tasks_text = NULL;
            if (options->distmem) {
                count_local_src_tasks_text = malloc(strlen("if (pi_(, nprocs) == my_rank) count;")
                        +5+strlen(local_indices[i])+1);
                sprintf(count_local_src_tasks_text, "if (pi_%d(%s, nprocs) == my_rank) count++;",
                        src_loop_num, local_indices[i]);
            }
            else {
                count_local_src_tasks_text = malloc(16);
                strcpy(count_local_src_tasks_text, "count++;");
            }
            pluto_add_stmt(count_local_src_tasks,tdpoly,local_trans[i],local_iters[i],count_local_src_tasks_text, IN_FUNCTION);

            for (j=0; j<src_copy_level; j++) {
                free(local_iters[i][j]);
            }

            free(count_local_src_tasks_text);
            free(local_indices[i]);
            free(local_iters[i]);
            pluto_matrix_free(local_trans[i]);
            pluto_constraints_free(tdpoly);
        }

        pluto_pad_stmt_transformations(count_local_src_tasks);
        if (count_local_src_tasks->nstmts >= 1) {
            assert(count_local_src_tasks->stmts[0]->trans->nrows == count_local_src_tasks->num_hyperplanes);
        }
        pluto_separate_stmts(count_local_src_tasks, count_local_src_tasks->stmts, count_local_src_tasks->nstmts, 0, 0);
    }


    FILE *cloogfp = NULL;

    if (options->distmem) {
        fprintf(outfp, "int get_num_remote_src_tasks_%d(", loop_num);
        fprintf(headerfp, "int get_num_remote_src_tasks_%d(", loop_num);
        for (i=0; i<dest_copy_level+prog->npar; i++)    {
            if (i!=0) {
                fprintf(outfp, ", ");
                fprintf(headerfp, ", ");
            }
            if (i<=dest_copy_level-1) {
                fprintf(outfp, "int ts%d", i+1);
                fprintf(headerfp, "int ts%d", i+1);
            }
            else {
                fprintf(outfp, "int %s", prog->params[i-dest_copy_level]);
                fprintf(headerfp, "int %s", prog->params[i-dest_copy_level]);
            }
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
        fprintf(headerfp, ", int my_rank, int nprocs);\n");
        fprintf(outfp, "\nint count = 0;\n");
        cloogfp = fopen("count_remote_src_tasks.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, count_remote_src_tasks);
        rewind(cloogfp);
        generate_declarations(count_remote_src_tasks, outfp);
        pluto_gen_cloog_code(count_remote_src_tasks, -1, -1, cloogfp, outfp);
        fclose(cloogfp);
        for (i=0; i<count_remote_src_tasks->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
        fprintf(outfp, "\nreturn count;\n");
        fprintf(outfp, "}\n\n");
    }

    fprintf(outfp, "int get_num_local_src_tasks_%d(", loop_num);
    fprintf(headerfp, "int get_num_local_src_tasks_%d(", loop_num);
    for (i=0; i<dest_copy_level+prog->npar; i++)    {
        if (i!=0) {
            fprintf(outfp, ", ");
            fprintf(headerfp, ", ");
        }
        if (i<=dest_copy_level-1) {
            fprintf(outfp, "int ts%d", i+1);
            fprintf(headerfp, "int ts%d", i+1);
        }
        else {
            fprintf(outfp, "int %s", prog->params[i-dest_copy_level]);
            fprintf(headerfp, "int %s", prog->params[i-dest_copy_level]);
        }
    }
    if (options->distmem) {
        fprintf(outfp, ", int my_rank, int nprocs");
        fprintf(headerfp, ", int my_rank, int nprocs");
    }
    fprintf(outfp, ")\n{\n");
    fprintf(headerfp, ");\n");
    fprintf(outfp, "\nint count = 0;\n");
    cloogfp = fopen("count_local_src_tasks.cloog", "w+");
    pluto_gen_cloog_file(cloogfp, count_local_src_tasks);
    rewind(cloogfp);
    generate_declarations(count_local_src_tasks, outfp);
    pluto_gen_cloog_code(count_local_src_tasks, -1, -1, cloogfp, outfp);
    fclose(cloogfp);
    for (i=0; i<count_local_src_tasks->nstmts; i++) {
        fprintf(outfp, "#undef S%d\n", i+1);
    }
    fprintf(outfp, "\nreturn count;\n");
    fprintf(outfp, "}\n\n");

    if (options->distmem) {
        pluto_prog_free(count_remote_src_tasks);
    }
    pluto_prog_free(count_local_src_tasks);
}


/*
 * generate_sigma: Generates receivers for a given tile for a particular
 * stmt_access_pair
 * 'copy_level': outer schedule rows of all loops
 * loop_num: source loop
 * 'copy_level[loop_num]' outer schedule rows of source loop will be treated as
 * parameters along with global parameters
 * pi_mappings: mapping from stmt id to parallel loop number or pi function #
 */
void generate_sigma(struct stmt_access_pair **wacc_stmts,
        int naccs, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings, FILE *outfp, FILE *headerfp)
{
    int i, j;
    char *acc_name;
    int src_copy_level = copy_level[loop_num];
    int is_sigma_required = !(options->commopt_fop || options->commopt_foifi);
    int is_receiver_required = options->dynschedule;

    assert(naccs != 0);

    acc_name = wacc_stmts[0]->acc->name;

    PlutoProg *sigma = NULL;
    if (is_sigma_required) {
		sigma = pluto_prog_alloc();

		// printf("Compute sigma\n");

		for (i=0; i<src_copy_level; i++) {
			char param[6];
			sprintf(param, "ts%d",i+1);
			pluto_prog_add_param(sigma, param, sigma->npar);
		}

		for (i=0; i<prog->npar; i++) {
			pluto_prog_add_param(sigma, prog->params[i], sigma->npar);
		}

#if 0
    for (i=0;i<src_copy_level; i++) {
        pluto_prog_add_hyperplane(sigma,0,H_LOOP);
    }
#endif
    }

    PlutoProg *is_receiver = NULL;
    if (is_receiver_required) {
		is_receiver = pluto_prog_alloc();

		for (i=0; i<src_copy_level; i++) {
			char param[6];
			sprintf(param, "ts%d",i+1);
			pluto_prog_add_param(is_receiver, param, is_receiver->npar);
		}

		for (i=0; i<prog->npar; i++) {
			pluto_prog_add_param(is_receiver, prog->params[i], is_receiver->npar);
		}
    }

    int broadcast;
    {
        int num_dest_loops;
        PlutoConstraints *tdpoly_list[prog->ndeps];
        int dep_loop_num_list[prog->ndeps];
        PlutoMatrix *trans[prog->ndeps];
        char **iters[prog->ndeps];
        char *indices[prog->ndeps];

        // for remote dependences only i.e., RAW dependences only
        generate_sigma_common(wacc_stmts, naccs, copy_level, prog, loop_num, pi_mappings, 1,
                &broadcast,
                &num_dest_loops,
                tdpoly_list,
                dep_loop_num_list,
                trans,
                iters,
                indices);

        if (!broadcast) {
            for (i=0; i<num_dest_loops; i++) {
                PlutoConstraints *tdpoly = tdpoly_list[i];
                int dep_loop_num = dep_loop_num_list[i];
                int dest_copy_level = copy_level[dep_loop_num];

                char *sigma_text = NULL;
                if (is_sigma_required) {
                    sigma_text = malloc(strlen("__p = pi_(,nprocs); if (__p != my_rank) { receiver_list[__p] = 1; }")
                            + 5 + strlen(indices[i])+1);
                    sprintf(sigma_text, "__p = pi_%d(%s,nprocs); if (__p != my_rank) { receiver_list[__p] = 1; }",
                            dep_loop_num, indices[i]);
                    pluto_add_stmt(sigma,tdpoly,trans[i],iters[i],sigma_text, IN_FUNCTION);
                }

                char *is_receiver_text = NULL;
                if (is_receiver_required) {
                    is_receiver_text = malloc(strlen("if (pi_(, nprocs) == my_rank) return 1;")
                            +5+strlen(indices[i])+1);
                    sprintf(is_receiver_text, "if (pi_%d(%s, nprocs) == my_rank) return 1;",
                            dep_loop_num, indices[i]);
                    pluto_add_stmt(is_receiver,tdpoly,trans[i],iters[i],is_receiver_text, IN_FUNCTION);
                }

                for (j=0; j<dest_copy_level; j++) {
                    free(iters[i][j]);
                }

                // if (is_sigma_required) pluto_stmt_print(stdout, sigma->stmts[0]);

                if (is_sigma_required) free(sigma_text);
                if (is_receiver_required) free(is_receiver_text);
                free(indices[i]);
                free(iters[i]);
                pluto_matrix_free(trans[i]);
                pluto_constraints_free(tdpoly);
            }

            if (is_sigma_required) {
                pluto_pad_stmt_transformations(sigma);
                // pluto_prog_print(sigma);
                if (sigma->nstmts >= 1) {
                    assert(sigma->stmts[0]->trans->nrows == sigma->num_hyperplanes);
                }

                pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts, 0, 0);
            }

            if (is_receiver_required) {
                pluto_pad_stmt_transformations(is_receiver);
                if (is_receiver->nstmts >= 1) {
                    assert(is_receiver->stmts[0]->trans->nrows == is_receiver->num_hyperplanes);
                }

                pluto_separate_stmts(is_receiver, is_receiver->stmts, is_receiver->nstmts, 0, 0);
            }
        }
    }

    FILE *cloogfp = NULL;

    if (is_sigma_required) {
		fprintf(outfp, "void sigma_%s_%d(", acc_name, loop_num);
		fprintf(headerfp, "void sigma_%s_%d(", acc_name, loop_num);
		for (i=0; i<src_copy_level+prog->npar; i++)    {
			if (i!=0) {
			    fprintf(outfp, ", ");
			    fprintf(headerfp, ", ");
			}
			if (i<=src_copy_level-1) {
			    fprintf(outfp, "int ts%d", i+1);
			    fprintf(headerfp, "int ts%d", i+1);
			}
			else {
			    fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
			    fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
			}
		}
		fprintf(outfp, ", int my_rank, int nprocs, int *receiver_list)\n{\n");
		fprintf(headerfp, ", int my_rank, int nprocs, int *receiver_list);\n");

		fprintf(outfp, "int __p;\n");

		if (broadcast) {
			fprintf(outfp, "\nfor (__p=0; __p<nprocs; __p++) if (__p != my_rank) receiver_list[__p] = 1;\n");
		}
		else {
			cloogfp = fopen("sigma.cloog", "w+");
			pluto_gen_cloog_file(cloogfp, sigma);
			rewind(cloogfp);
			generate_declarations(sigma, outfp);
			pluto_gen_cloog_code(sigma, -1, -1, cloogfp, outfp);
			fclose(cloogfp);

			for (i=0; i<sigma->nstmts; i++) {
				fprintf(outfp, "#undef S%d\n", i+1);
			}
		}

		fprintf(outfp, "}\n\n");
    }

    if (is_receiver_required) {
		fprintf(outfp, "int is_receiver_%s_%d(", acc_name, loop_num);
		fprintf(headerfp, "int is_receiver_%s_%d(", acc_name, loop_num);
		for (i=0; i<src_copy_level+prog->npar; i++)    {
			if (i!=0) {
			    fprintf(outfp, ", ");
			    fprintf(headerfp, ", ");
			}
			if (i<=src_copy_level-1) {
			    fprintf(outfp, "int ts%d", i+1);
			    fprintf(headerfp, "int ts%d", i+1);
			}
			else {
			    fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
			    fprintf(headerfp, "int %s", prog->params[i-src_copy_level]);
			}
		}
		fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
		fprintf(headerfp, ", int my_rank, int nprocs);\n");

		if (broadcast) {
			fprintf(outfp, "\nreturn 1;\n");
		}
		else {
			cloogfp = fopen("is_receiver.cloog", "w+");
			pluto_gen_cloog_file(cloogfp, is_receiver);
			rewind(cloogfp);
			generate_declarations(is_receiver, outfp);
			pluto_gen_cloog_code(is_receiver, -1, -1, cloogfp, outfp);
			fclose(cloogfp);

			for (i=0; i<is_receiver->nstmts; i++) {
				fprintf(outfp, "#undef S%d\n", i+1);
			}
            fprintf(outfp, "\nreturn 0;\n");
		}

		fprintf(outfp, "}\n\n");
    }
    
    if (is_sigma_required) pluto_prog_free(sigma);
    if (is_receiver_required) pluto_prog_free(is_receiver);
}

PlutoConstraints* get_receiver_tiles_of_dep(Dep *dep, 
        int src_copy_level, int dest_copy_level, PlutoProg *prog, int use_src_unique_dpolytope)
{
    int total_copy_level = src_copy_level + dest_copy_level;

    Stmt *src = prog->stmts[dep->src];
    Stmt *dest = prog->stmts[dep->dest];

    PlutoConstraints *tdpoly;
    if (!use_src_unique_dpolytope) tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
    else tdpoly = pluto_constraints_dup(dep->src_unique_dpolytope);
    assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

    pluto_constraints_project_out(tdpoly, src_copy_level, src->trans->nrows-src_copy_level);
    pluto_constraints_project_out(tdpoly, total_copy_level, dest->trans->nrows-dest_copy_level);
    assert(tdpoly->ncols == total_copy_level + prog->npar+1);
    return tdpoly;
}

/*
 * get_receiver_tiles: Returns receiver tiles for the given stmt_access_pairs and
 * for the given dependent loop
 * 'src_copy_level' outer schedule rows on the source statement and
 * 'dest_copy_level' outer schedule rows on the destination statement 
 * will be treated as * parameters along with global parameters 
 * pi_mappings: mapping from stmt id to parallel loop number or pi function #
 */
PlutoConstraints* get_receiver_tiles(struct stmt_access_pair **wacc_stmts, 
        int naccs, int src_copy_level, int dest_copy_level, PlutoProg *prog, int dep_loop_num, int *pi_mappings)
{
    int i, j;
    int total_copy_level = src_copy_level + dest_copy_level;

    assert(naccs != 0);

    PlutoConstraints *receiver_tiles = pluto_constraints_empty(total_copy_level+prog->npar+1);
    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

        Stmt *dest = prog->stmts[dep->dest];
        
        // find receiver tiles only for the given dependent loop
        if (pi_mappings[dest->id] != dep_loop_num) continue;

        /* If the dependence doesn't originate from this access */
        for (j=0; j<naccs; j++) {
            if (dep->src_acc == wacc_stmts[j]->acc) break;
        }
        if (j==naccs)   continue;

        PlutoConstraints *tdpoly = get_receiver_tiles_of_dep(dep, src_copy_level, dest_copy_level, prog, 0);
        pluto_constraints_unionize(receiver_tiles, tdpoly);
        pluto_constraints_free(tdpoly);
    }
    
    return receiver_tiles;
}


/*
 * Tau: generates senders for a given tile
 * 'copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 */
// not used currently
void generate_tau(struct stmt_access_pair **racc_stmts,
        int naccs, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings, FILE *outfp, FILE *headerfp)
{
    int i, j;
    char *acc_name;
    int dest_copy_level = copy_level[loop_num];

    assert(naccs != 0);

    acc_name = racc_stmts[0]->acc->name;

    PlutoProg *count_sending_tasks = NULL;
    count_sending_tasks = pluto_prog_alloc();

    // printf("Compute tau\n");

    for (i=0; i<dest_copy_level; i++) {
        char param[6];
        sprintf(param, "ts%d",i+1);
        pluto_prog_add_param(count_sending_tasks, param, count_sending_tasks->npar);
    }

    for (i=0; i<prog->npar; i++) {
        pluto_prog_add_param(count_sending_tasks, prog->params[i], count_sending_tasks->npar);
    }

#if 0
    for (i=0;i<src_copy_level; i++) {
        pluto_prog_add_hyperplane(count_sending_tasks,0,H_LOOP);
    }
#endif

    int num_src_loops;
    PlutoConstraints *tdpoly_list[prog->ndeps];
    int src_loop_num_list[prog->ndeps];
    PlutoMatrix *trans[prog->ndeps];
    char **iters[prog->ndeps];
    char *indices[prog->ndeps];

    generate_tau_common(racc_stmts, naccs, copy_level, prog, loop_num, pi_mappings,
            &num_src_loops,
            tdpoly_list,
            src_loop_num_list,
            trans,
            iters,
            indices);

    for (i=0; i<num_src_loops; i++) {
        PlutoConstraints *tdpoly = tdpoly_list[i];
        int src_loop_num = src_loop_num_list[i];
        int src_copy_level = copy_level[src_loop_num];

        char *count_sending_tasks_text = NULL;
        count_sending_tasks_text = malloc(strlen("if (pi_(,nprocs) != my_rank) count++;")
                + 5 + strlen(indices[i])+1);
        sprintf(count_sending_tasks_text, "if (pi_%d(%s,nprocs) != my_rank) count++;",
                src_loop_num, indices[i]);
        pluto_add_stmt(count_sending_tasks,tdpoly,trans[i],iters[i],count_sending_tasks_text, IN_FUNCTION);

        for (j=0; j<src_copy_level; j++) {
            free(iters[i][j]);
        }

        // pluto_stmt_print(stdout, count_sending_tasks->stmts[0]);

        free(count_sending_tasks_text);
        free(indices[i]);
        free(iters[i]);
        pluto_matrix_free(trans[i]);
        pluto_constraints_free(tdpoly);
    }

    pluto_pad_stmt_transformations(count_sending_tasks);
    // pluto_prog_print(count_sending_tasks);
    if (count_sending_tasks->nstmts >= 1) {
        assert(count_sending_tasks->stmts[0]->trans->nrows == count_sending_tasks->num_hyperplanes);
    }

    pluto_separate_stmts(count_sending_tasks, count_sending_tasks->stmts, count_sending_tasks->nstmts, 0, 0);

    FILE *cloogfp = NULL;

    fprintf(outfp, "int count_sending_tasks_%s_%d(", acc_name, loop_num);
    fprintf(headerfp, "int count_sending_tasks_%s_%d(", acc_name, loop_num);
    for (i=0; i<dest_copy_level+prog->npar; i++)    {
        if (i!=0) {
            fprintf(outfp, ", ");
            fprintf(headerfp, ", ");
        }
        if (i<=dest_copy_level-1) {
            fprintf(outfp, "int ts%d", i+1);
            fprintf(headerfp, "int ts%d", i+1);
        }
        else {
            fprintf(outfp, "int %s", prog->params[i-dest_copy_level]);
            fprintf(headerfp, "int %s", prog->params[i-dest_copy_level]);
        }
    }
    fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
    fprintf(headerfp, ", int my_rank, int nprocs);\n");

    fprintf(outfp, "int count = 0;\n");

    cloogfp = fopen("count_sending_tasks.cloog", "w+");
    pluto_gen_cloog_file(cloogfp, count_sending_tasks);
    rewind(cloogfp);
    generate_declarations(count_sending_tasks, outfp);
    pluto_gen_cloog_code(count_sending_tasks, -1, -1, cloogfp, outfp);
    fclose(cloogfp);

    for (i=0; i<count_sending_tasks->nstmts; i++) {
        fprintf(outfp, "#undef S%d\n", i+1);
    }

    fprintf(outfp, "return count;\n}\n\n");
    pluto_prog_free(count_sending_tasks);
}


void pluto_constraints_list_dep_add(PlutoConstraintsList *list, Dep *dep){

    assert(list != NULL);

    PlutoConstraintsList *new = pluto_constraints_list_alloc(pluto_constraints_dup(dep->dpolytope));
    new->next = list->next;
    list->next = new;

    new->deps = pluto_deps_list_dup(list->deps);
    pluto_deps_list_append(new->deps,dep);

    return;
}



/*adds a new element to constraints list
 * also adds the deps lsit
*/
void pluto_constraints_list_add_flowout(PlutoConstraintsList *list,const PlutoConstraints *cst, PlutoDepList *dep_list)
{
    assert(list != NULL);
    assert(cst != NULL);

    if(list->deps == NULL){
    	assert(list->next == NULL);
    	list->deps = dep_list;
    	list->constraints = pluto_constraints_dup(cst);
    	return;
    }

    PlutoConstraintsList *new = pluto_constraints_list_alloc(pluto_constraints_dup(cst));
    new->next = list->next;
    new->deps = dep_list;
    list->next = new;

    return;
}

Dep *pluto_dependence_dup(Dep *dep, PlutoConstraints *tdpoly, PlutoConstraints *src_unique_dpolytope){

	Dep *d = (Dep *)malloc(sizeof(Dep));

	d->depsat_poly = dep->depsat_poly;
	d->dest = dep->dest;
	d->dest_acc = dep->dest_acc;
	d->dirvec = dep->dirvec;
	d->dpolytope = pluto_constraints_dup(tdpoly);
	d->src_unique_dpolytope = pluto_constraints_dup(src_unique_dpolytope);
	d->id = dep->id;
	d->satisfaction_level = dep->satisfaction_level;
	d->satisfied = dep->satisfied;
	d->satvec = dep->satvec;
	d->src = dep->src;
	d->src_acc = dep->src_acc;
	d->type = dep->type;

	return d;

}

void set_dependence_tdpoly(Dep *d, PlutoConstraints *tdpoly){
	pluto_constraints_free(d->dpolytope);
	d->dpolytope = pluto_constraints_dup(tdpoly);
}

void pluto_dep_list_add(PlutoDepList* list, Dep *d){

    PlutoDepList *new = (PlutoDepList *)malloc(sizeof(PlutoDepList));
    new->dep = d;
    new->next = list->next;
    list->next = new;

    return;
}

void pluto_dep_list_list_add(PlutoDepListList* list, PlutoDepList *ls){

    PlutoDepListList *new = (PlutoDepListList *)malloc(sizeof(PlutoDepListList));
    new->dep_list = ls;
    new->next = list->next;
    list->next = new;

    return;
}

void compute_common_data_iterations(PlutoConstraints *td1, PlutoConstraints *td2,
		Dep *d1, Dep *d2, PlutoProg *prog, PlutoConstraints **intersect1, PlutoConstraints **intersect2)
{


	//if(d1->src == d2->src)
	if(0)
	{
		assert(td1->ncols == td2->ncols);

		*intersect1 = pluto_constraints_intersection(td1, td2);
		if(pluto_constraints_is_empty(*intersect1))
			*intersect1 = *intersect2 = NULL;
		else
			*intersect2 = pluto_constraints_dup(*intersect1);

		return;
	}
	else {

		//try finding the dependence between source iterators with their respective access functions
		PlutoConstraints *tdpoly = pluto_find_dependence(td1, td2, d1, d2, prog, NULL);

		//if there is no dependence try finding dep with next item in list
		if(tdpoly == NULL){
			*intersect1 = *intersect2 = NULL;
			return;
		}

		Stmt *src1 = prog->stmts[d1->src];
		Stmt *src2 = prog->stmts[d2->src];


		//get the both the source iterations for which the deps exits
		*intersect1 = pluto_constraints_dup(tdpoly);
		pluto_constraints_project_out(*intersect1,src1->dim ,src2->dim );

		*intersect2 = pluto_constraints_dup(tdpoly);
		pluto_constraints_project_out(*intersect2, 0, src1->dim);

	}

	return;
}

//Recursively do the dependence intersection and difference on all the elements of list with dep_to_split and put the
//intersected dependence in new list
void intersect_deps_in_list(PlutoDepList *list, PlutoDepList *intersect_list , Dep *dep_to_split, PlutoProg *prog){

	if(list == NULL)
		return;

    PlutoConstraints *intersect1, *intersect2, *diff2;

    Stmt *src = prog->stmts[dep_to_split->src];

	PlutoConstraints *td1, *td2;
	PlutoConstraints *tdpoly2;
	Dep *dep2, *curr_dep;
	int i = 0;

	curr_dep = list->dep;

	td1 = pluto_constraints_dup(dep_to_split->dpolytope);
	td2 = pluto_constraints_dup(curr_dep->dpolytope);

	//only source iterators are common for both the dependence
	//dest stmt could be different
	int td1_src_ncols = src->trans->nrows;
	int td1_dest_ncols = td1->ncols - src->trans->nrows - prog->npar-1;

	int td2_src_ncols = src->trans->nrows;
	int td2_dest_ncols = td2->ncols - src->trans->nrows - prog->npar-1;

	//project out dest iterators
	pluto_constraints_project_out(td1, td1_src_ncols, td1_dest_ncols);
	pluto_constraints_project_out(td2, td2_src_ncols, td2_dest_ncols);

	compute_common_data_iterations(td1, td2, dep_to_split, curr_dep, prog, &intersect1, &intersect2);

	//if there is no dependence return
	if(intersect1 == NULL && intersect2 == NULL){
		return;
	}
	else {

		tdpoly2 = pluto_constraints_dup(intersect2);

		for(i=0; i<td2_dest_ncols; i++)
			pluto_constraints_add_dim(tdpoly2, 0, NULL);

		pluto_constraints_intersect(tdpoly2,curr_dep->dpolytope);
		dep2 = pluto_dependence_dup(curr_dep, tdpoly2, tdpoly2);

		pluto_deps_list_append(intersect_list, dep2);

		//compute the remaining source iterators from curr_dep
		diff2 = pluto_constraints_difference(td2, intersect2);

		for(i=0; i<td2_dest_ncols; i++)
			pluto_constraints_add_dim(diff2, td2_src_ncols, NULL);

		pluto_constraints_intersect(curr_dep->dpolytope,diff2);

		intersect_deps_in_list(list->next, intersect_list, dep_to_split, prog);

		return;
	}
}

PlutoDepListList *pluto_dep_list_list_alloc(){
	PlutoDepListList *ls = (PlutoDepListList *) malloc(sizeof(PlutoDepListList));
	ls->dep_list = NULL;
	ls->next = NULL;
	return ls;
}

int check_for_disjoint_sets(PlutoConstraintsList *list){
	PlutoConstraintsList *outer, *inner;
	outer = inner = list;
	while(outer != NULL){
		inner = outer->next;
		while(inner != NULL){
			if(!pluto_constraints_is_empty(pluto_constraints_intersection(outer->constraints, inner->constraints)))
					return 0;

			inner = inner->next;
		}

		outer = outer->next;
	}

	return 1;
}

int check_for_disjoint_deps(PlutoDepListList *list, PlutoProg *prog){
	PlutoDepListList *outer, *inner;
	Dep *dep;
	Stmt *src, *dst;
	int td_src_ncols, td_dest_ncols;
	PlutoConstraints *tdpoly, *tdpoly_outer, *intersection;


	outer = inner = list;
	while(outer != NULL){
		dep = outer->dep_list->dep;
		tdpoly_outer = pluto_constraints_dup(dep->dpolytope);

		src = prog->stmts[dep->src];
		dst = prog->stmts[dep->dest];

		td_src_ncols = src->trans->nrows;
		td_dest_ncols = dst->trans->nrows;

		pluto_constraints_project_out(tdpoly_outer, td_src_ncols, td_dest_ncols);

		inner = outer->next;
		while(inner != NULL){
			dep = inner->dep_list->dep;
			tdpoly = pluto_constraints_dup(dep->dpolytope);

			src = prog->stmts[dep->src];
			dst = prog->stmts[dep->dest];

			td_src_ncols = src->trans->nrows;
			td_dest_ncols = tdpoly->ncols - src->trans->nrows - prog->npar-1;

			pluto_constraints_project_out(tdpoly, td_src_ncols, td_dest_ncols);

			intersection = pluto_constraints_intersection(tdpoly_outer, tdpoly);

			if(!pluto_constraints_is_empty(intersection))
					return 0;

			pluto_constraints_free(tdpoly);
			inner = inner->next;
		}

		pluto_constraints_free(tdpoly_outer);
		outer = outer->next;
	}

	return 1;
}

int check_for_valid_data(PlutoDepListList *list, int copy_level, PlutoProg *prog){

	PlutoDepListList *curr = list;
	Dep *dep;
	PlutoConstraints *tdpoly, *data, *datac;

	while(curr != NULL && curr->dep_list != NULL){

		//Do the compute data region to only one of Dep in list
		//As all the deps in list access same data
		dep = curr->dep_list->dep;
		tdpoly = pluto_constraints_dup(dep->dpolytope);

		Stmt *src = prog->stmts[dep->src];

		int td_src_ncols = src->trans->nrows;
		int td_dest_ncols = tdpoly->ncols - src->trans->nrows - prog->npar-1;

		pluto_constraints_project_out(tdpoly, td_src_ncols, td_dest_ncols);


		data = pluto_compute_region_data(src, tdpoly, dep->src_acc,
				copy_level, prog);

		IF_DEBUG(print_polylib_visual_sets("data1", data));
		PlutoDepList *d = curr->dep_list->next;

		while(d != NULL){
			dep = d->dep;
			tdpoly = pluto_constraints_dup(dep->dpolytope);

			src = prog->stmts[dep->src];

			td_src_ncols = src->trans->nrows;
			td_dest_ncols = tdpoly->ncols - src->trans->nrows - prog->npar-1;

			pluto_constraints_project_out(tdpoly, td_src_ncols, td_dest_ncols);

			datac = pluto_compute_region_data(src, tdpoly, dep->src_acc,
					copy_level, prog);

			IF_DEBUG(print_polylib_visual_sets("data2", datac));

			if(!pluto_constraints_is_empty(pluto_constraints_difference(datac, data))){
				return 0;
			}

			d = d->next;
		}

		curr = curr->next;
	}

	return 1;

}

void split_dependences(PlutoDepListList *dep_list_list,
        PlutoConstraints *curr_dep_poly, Dep *dep, PlutoProg *prog )
{
    assert(curr_dep_poly);

    /*for initial case when atomic dep is empty
    */
	Dep *dep_to_split = pluto_dependence_dup(dep, curr_dep_poly, curr_dep_poly);
    if(dep_list_list->dep_list == NULL) {
        dep_list_list->dep_list  = pluto_dep_list_alloc(dep_to_split);
        return;
    }

    PlutoDepListList *curr = dep_list_list, *next;
    PlutoConstraints *intersect1, *intersect2, *diff1, *diff2;
    int diff1_empty, diff2_empty;

    Stmt *src1 = prog->stmts[dep_to_split->src];

	PlutoConstraints *td1, *td2;
	PlutoConstraints *tdpoly1, *tdpoly2;
	Dep *dep1, *dep2, *curr_dep;
	dep1 = dep;
	int i = 0, k=0;

    while(curr != NULL) {

    	k++;

    	//n = check_for_valid_data(dep_list_list, 2, prog);
        next = curr->next;
        curr_dep = curr->dep_list->dep;

		Stmt *src2 = prog->stmts[curr_dep->src];

        td1 = pluto_constraints_dup(dep_to_split->dpolytope);
        td2 = pluto_constraints_dup(curr_dep->dpolytope);

		int td1_src_ncols = src1->trans->nrows;
		int td1_dest_ncols = td1->ncols - src1->trans->nrows - prog->npar-1;

		int td2_src_ncols = src2->trans->nrows;
		int td2_dest_ncols = td2->ncols - src2->trans->nrows - prog->npar-1;

		//project out dest iterators
		pluto_constraints_project_out(td1, td1_src_ncols, td1_dest_ncols);
		pluto_constraints_project_out(td2, td2_src_ncols, td2_dest_ncols);


		IF_DEBUG(print_polylib_visual_sets("td1", td1));
		IF_DEBUG(print_polylib_visual_sets("td2", td2));

		compute_common_data_iterations(td1, td2, dep_to_split, curr_dep, prog, &intersect1, &intersect2);

		//if there is no dependence try finding dep with next item in list
		if(intersect1 == NULL && intersect2 == NULL){
			curr = next;
			continue;
		}

		IF_DEBUG(print_polylib_visual_sets("intersect1", intersect1));
		IF_DEBUG(print_polylib_visual_sets("intersect2", intersect2));
		//compute the remaining source iterators from curr_dep
		diff1 = pluto_constraints_difference(td1, intersect1);
		diff2 = pluto_constraints_difference(td2, intersect2);

		IF_DEBUG(print_polylib_visual_sets("diff1", diff1));
		IF_DEBUG(print_polylib_visual_sets("diff2", diff2));

		diff1_empty = pluto_constraints_is_empty(diff1);
		diff2_empty = pluto_constraints_is_empty(diff2);


		if(diff1_empty && diff2_empty){
			//add dep_to_split to curr dep list
			pluto_dep_list_add(curr->dep_list, dep_to_split);
			//TODO: do the clean up
			return;
		}

		else if (!diff1_empty && diff2_empty){

			//add dep_to_split to curr dep list
			tdpoly1 = pluto_constraints_dup(intersect1);

			for(i=0; i<td1_dest_ncols; i++)
				pluto_constraints_add_dim(tdpoly1, td1_src_ncols, NULL);

			pluto_constraints_intersect(tdpoly1,dep_to_split->dpolytope);
			dep1 = pluto_dependence_dup(dep_to_split, tdpoly1, tdpoly1);

			pluto_dep_list_add(curr->dep_list, dep1);

			for(i=0; i<td2_src_ncols; i++)
				pluto_constraints_add_dim(diff1, td1_src_ncols, NULL);

			//diff1 has the source iterators for next iteration
			pluto_constraints_intersect(dep_to_split->dpolytope,diff1);

			//TODO: Do the clean up
			curr = next;
			continue;

		}
		else if(!diff2_empty) {


			tdpoly1 = pluto_constraints_dup(intersect1);

			for(i=0; i<td1_dest_ncols; i++)
				pluto_constraints_add_dim(tdpoly1, td1_src_ncols, NULL);

			pluto_constraints_intersect(tdpoly1,dep_to_split->dpolytope);
			dep1 = pluto_dependence_dup(dep_to_split, tdpoly1, tdpoly1);

			tdpoly2 = pluto_constraints_dup(intersect2);
			for(i=0; i<td2_dest_ncols; i++)
				pluto_constraints_add_dim(tdpoly2, td2_src_ncols, NULL);

			pluto_constraints_intersect(tdpoly2,curr_dep->dpolytope);
			dep2 = pluto_dependence_dup(curr_dep, tdpoly2, tdpoly2);


			//Add the new dep list for the new intersected deps which access same data
			PlutoDepList *list = pluto_dep_list_alloc(dep1);
			pluto_deps_list_append(list, dep2);

			pluto_dep_list_list_add(curr, list);


			for(i=0; i<td2_dest_ncols; i++)
				pluto_constraints_add_dim(diff2, td2_src_ncols, NULL);

			pluto_constraints_intersect(curr_dep->dpolytope,diff2);

			//Recursively split the remaining dependences in the curr dep list
			intersect_deps_in_list(curr->dep_list->next, list, dep_to_split, prog);

			for(i=0; i<td1_dest_ncols; i++)
				pluto_constraints_add_dim(diff1, td1_src_ncols, NULL);

			//diff1 has the source iterators for next iteration
			pluto_constraints_intersect(dep_to_split->dpolytope,diff1);

			if(diff1_empty)
				assert(pluto_constraints_is_empty(dep_to_split->dpolytope));

			//TODO: Add the cleanup code
			curr = next;

		}

    }

    //If it is dep that has some source iterators that do not intersect with any dep in list
    //add them to list
    if(!pluto_constraints_is_empty(dep_to_split->dpolytope)){

		PlutoDepList *list = pluto_dep_list_alloc(dep_to_split);

		pluto_dep_list_list_add(dep_list_list, list);
    }

    return;
}

void split_flow_out_set(PlutoProg *prog, PlutoConstraintsList *atomic_flowouts,
        PlutoConstraints *curr_flowout, Dep *dep)
{
    int i;
    assert(atomic_flowouts);
    assert(curr_flowout);

    /*for inital case when atomic flow out is empty
    */
    if(atomic_flowouts->constraints == NULL) {
        assert(atomic_flowouts->deps == NULL);
        atomic_flowouts->constraints = pluto_constraints_dup(curr_flowout);
        atomic_flowouts->deps = pluto_dep_list_alloc(dep);
        return;
    }

    PlutoConstraintsList *curr = atomic_flowouts, *next;
    PlutoConstraints *cst_intersect, *cst1, *cst2;
    PlutoConstraints *intersect, *intersect1;
    int curr_constraints_empty, curr_flowout_empty;

    curr_flowout_empty = pluto_constraints_is_empty(curr_flowout);
    assert(!curr_flowout_empty);
    while ((!curr_flowout_empty) && (curr != NULL)) {
        IF_DEBUG(printf("before intersection \n"););
        intersect = pluto_constraints_intersection(curr_flowout, curr->constraints);
        
        cst1 = pluto_constraints_dup(curr->constraints);
        pluto_constraints_project_out(cst1, cst1->ncols-prog->npar-1, prog->npar);
        cst2 = pluto_constraints_dup(curr_flowout);
        pluto_constraints_project_out(cst2, cst2->ncols-prog->npar-1, prog->npar);
        cst_intersect = pluto_constraints_intersection(cst1, cst2);

        IF_DEBUG(printf("after intersection \n"););
        if(pluto_constraints_is_empty(cst_intersect)) {
            //assert(pluto_constraints_is_empty(iterations_intersect));
            IF_DEBUG(printf("Intersection is empty\n"););
            curr = curr->next;
            pluto_constraints_free(cst1);
            pluto_constraints_free(cst2);
            pluto_constraints_free(cst_intersect);
            continue;
        }

        IF_DEBUG(printf("intersection \n"););
        IF_DEBUG(pluto_constraints_print(stdout, intersect););

        for (i=0; i<prog->npar; i++) {
            pluto_constraints_add_dim(cst_intersect, cst_intersect->ncols-1, NULL);
        }

        intersect = pluto_constraints_intersection(cst_intersect, curr->constraints);
        intersect1 = pluto_constraints_intersection(cst_intersect, curr_flowout);
        pluto_constraints_unionize(intersect, intersect1);
        assert(!pluto_constraints_is_empty(intersect));

        pluto_constraints_subtract(curr->constraints, intersect);

        pluto_constraints_subtract(curr_flowout, intersect);

        IF_DEBUG(printf("diff1 \n"););

        IF_DEBUG(pluto_constraints_print(stdout, curr->constraints););
        IF_DEBUG(printf("diff2 \n"););
        IF_DEBUG(pluto_constraints_print(stdout, curr_flowout););
        curr_constraints_empty = pluto_constraints_is_empty(curr->constraints);
        curr_flowout_empty = pluto_constraints_is_empty(curr_flowout);

        if(curr_constraints_empty){ 
            IF_DEBUG(printf("both intersection and curr flow out are equal\n"););
            pluto_constraints_unionize(curr->constraints, intersect);

            //add the dependency
            pluto_deps_list_append(curr->deps, dep);
            curr = curr->next;
        } else {
            next = curr->next;
            //Add insterction to the list and add parents deps
            pluto_constraints_list_add(curr, intersect, dep, 1);
            
            curr = next;
        }

        pluto_constraints_free(intersect);
        pluto_constraints_free(intersect1);
        pluto_constraints_free(cst_intersect);
        pluto_constraints_free(cst1);
        pluto_constraints_free(cst2);
    }

    //At the end append the remaining set, and current dependence
    if(!curr_flowout_empty) {
        pluto_constraints_list_add(atomic_flowouts, curr_flowout, dep, 0);
        pluto_constraints_free(curr_flowout);
    }

    return;
}

/*
 * 'copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 *
 * Output format:  copy_level, acc->nrows, prog->npar + 1
 */
PlutoConstraints *compute_flow_in_of_dep(Dep *dep, 
        int copy_level, PlutoProg *prog, int use_src_unique_dpolytope)
{
    int j;
    Stmt *src = prog->stmts[dep->src];
    Stmt *dest = prog->stmts[dep->dest];
    PlutoAccess *write_acc = dep->src_acc;

    // printf("Write access\n");
    // pluto_matrix_print(stdout, write_acc->mat);

    PlutoConstraints *tdpoly;
    if (!use_src_unique_dpolytope) tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
    else tdpoly = pluto_constraints_dup(dep->src_unique_dpolytope);
    assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

    //printf("tdpoly\n");
    //pluto_constraints_print(stdout, tdpoly);

    /* Dep source iterations corresponding to this target tile */
    /* tile_src format: [copy_level | src->trans->nrows + par + 1] */
    PlutoConstraints *tile_src = pluto_constraints_dup(tdpoly);

    PlutoConstraints *param_eq;

    /* No need to check whether source and target share all of the
     * copy_level loops; handled naturally since those scalar dimensions
     * would be part of the 'copy_level' dimensions. If the statements are 
     * separated somewhere, subtracting param_eq will not changed tdpoly
     */
    param_eq = get_context_equality(copy_level, src->trans->nrows, 
            tdpoly->ncols);

    /* Dep source iterations that fill within the same target tile */
    PlutoConstraints *tile_src_inside = 
        pluto_constraints_intersection(tile_src, param_eq);

    PlutoConstraints *scst;
    if (pluto_constraints_is_empty(tile_src_inside)) {
        /* All source iterations are outside the target tile */
        // target and source statements could still have the same set of surrounding loops
        // so incorrect to project out source copy_level loop iterators - retain them

        // project out target loop iterators inner to copy_level
        pluto_constraints_project_out(tile_src, src->trans->nrows+copy_level, 
                dest->trans->nrows-copy_level);
        
        // move target copy_level loop iterators to the beginning/top
        for (j=0; j<copy_level; j++) {
            pluto_constraints_add_dim(tile_src, 0, NULL);
        }
        for (j=0; j<copy_level; j++) {
            pluto_constraints_interchange_cols(tile_src, j, j+copy_level+src->trans->nrows);
        }
        for (j=0; j<copy_level; j++) {
            pluto_constraints_remove_dim(tile_src, copy_level+src->trans->nrows);
        }

        scst =  pluto_compute_region_data(src, tile_src, write_acc, 
                copy_level, prog);
    }else{ // some source iterations within the source tile
        // target and source statements definitely have outer copy_level surrounding loops common
        // so project out the source copy_level iterators too

        // project out source copy_level loop iterators and target loop iterators inner to copy_level
        pluto_constraints_project_out(tile_src_inside, 0, 
                copy_level);
        pluto_constraints_project_out(tile_src_inside, src->trans->nrows, 
                dest->trans->nrows-copy_level);
        //PlutoConstraints *data_inside = 
        //pluto_compute_region_data(src, tile_src_inside, write_acc, 
        //copy_level, prog);

        //IF_DEBUG(printf("Data inside\n"););
        //IF_DEBUG(pluto_constraints_print(stdout, data_inside););

        // project out source copy_level loop iterators and target loop iterators inner to copy_level
        pluto_constraints_project_out(tile_src, 0, 
                copy_level);
        pluto_constraints_project_out(tile_src, src->trans->nrows, 
                dest->trans->nrows-copy_level);

        // difference should be taken after the two sets are parameterized on the target copy_level iterators
        PlutoConstraints *tile_src_outside = pluto_constraints_difference(tile_src, 
                tile_src_inside);
        
        // move target copy_level loop iterators to the beginning/top
        for (j=0; j<copy_level; j++) {
            pluto_constraints_add_dim(tile_src_outside, 0, NULL);
        }
        for (j=0; j<copy_level; j++) {
            pluto_constraints_interchange_cols(tile_src_outside, j, j+src->trans->nrows);
        }
        for (j=0; j<copy_level; j++) {
            pluto_constraints_remove_dim(tile_src_outside, src->trans->nrows);
        }

        /* scst parameterized by the 
         * target tile copy_level params */
        scst = pluto_compute_region_data(src, tile_src_outside, write_acc, 
                copy_level, prog);
        pluto_constraints_free(tile_src_outside);
    }

    IF_DEBUG(printf("Data flowing in for dep %d\n", dep->id+1););
    IF_DEBUG(pluto_constraints_print(stdout, scst););

    pluto_constraints_free(tile_src);
    pluto_constraints_free(tile_src_inside);
    pluto_constraints_free(tdpoly);
    pluto_constraints_free(param_eq);
    return scst;
}

/*
 * 'copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 *
 * Output format:  copy_level, acc->nrows, prog->npar + 1
 */
PlutoConstraints *compute_flow_in(struct stmt_access_pair *racc_stmt, 
        int copy_level, PlutoProg *prog)
{
    int i;

    /* Stmt *stmt = racc_stmt->stmt; */
    PlutoAccess *racc = racc_stmt->acc;

    PlutoConstraints *uscst = pluto_constraints_empty(copy_level
            +racc->mat->nrows+prog->npar+1);

    IF_DEBUG(printf("Computing flow in set for %s\n", racc->name););

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

        // if (dep->dirvec[copy_level+1] == DEP_ZERO) continue;

        /* If the dependence doesn't originate from this access */
        if (dep->dest_acc != racc) continue;

        assert(dep->src_acc != NULL);

        PlutoConstraints *scst = compute_flow_in_of_dep(dep, copy_level, prog, 0);
        pluto_constraints_unionize(uscst, scst);
        pluto_constraints_free(scst);
    }
    //pluto_constraints_print(stdout, uscst);
    return uscst;
}

/*
 * 'src_copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 *
 * Output format:  src_copy_level, acc->nrows, prog->npar + 1
 */
PlutoConstraints *compute_flow_out_of_dep(Dep *dep, 
        int src_copy_level, int *copy_level, PlutoProg *prog, int split, PlutoConstraints **dcst1, int *pi_mappings)
{
    int j;
    Stmt *src = prog->stmts[dep->src];
    Stmt *dest = prog->stmts[dep->dest];
    PlutoAccess *read_acc = dep->dest_acc;

    int dep_loop_num = pi_mappings[dest->id];
    int dest_copy_level = copy_level[dep_loop_num];
    int total_copy_level = src_copy_level + dest_copy_level;

    // printf("Read access\n");
    // pluto_matrix_print(stdout, read_acc->mat);

    PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
    assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

    //printf("tdpoly\n");
    //pluto_constraints_print(stdout, tdpoly);

    /* Dep target iterations corresponding to this source tile */
    /* tile_dest format: [src_copy_level | dest->trans->nrows + par + 1] */
    PlutoConstraints *tile_dest = pluto_constraints_dup(tdpoly);

    PlutoConstraints *param_eq;

    /* No need to check whether source and target share all of the
     * src_copy_level loops; handled naturally since those scalar dimensions
     * would be part of the 'src_copy_level' dimensions. If the statements are 
     * separated somewhere, subtracting param_eq will not changed tiledest
     */
    param_eq = get_context_equality(src_copy_level, src->trans->nrows, 
            tdpoly->ncols);

    /* Dep target iterations that fill within the same source tile */
    PlutoConstraints *tile_dest_inside = 
        pluto_constraints_intersection(tile_dest, param_eq);

    PlutoConstraints *dcst;
    if (pluto_constraints_is_empty(tile_dest_inside)) {
        /* All target iterations are outside the source tile */
        // target and source statements could still have the same set of surrounding loops
        // so incorrect to project out target src_copy_level loop iterators - retain them
        // project out only source loop iterators inner to src_copy_level
        pluto_constraints_project_out(tile_dest, src_copy_level,
                src->trans->nrows-src_copy_level);

        dep->src_unique_dpolytope = pluto_constraints_dup(tdpoly);
        if (split) {
            assert(dcst1 != NULL);
            
            PlutoConstraints *copyleveleq, *tile_dest_same_copy_level;
            copyleveleq = pluto_constraints_alloc(1, tile_dest->ncols);
            pluto_constraints_add_equality(copyleveleq);
            copyleveleq->val[0][src_copy_level-1] = 1;
            copyleveleq->val[0][total_copy_level-1] = -1;

            tile_dest_same_copy_level = pluto_constraints_intersection(tile_dest, copyleveleq);

            if (!pluto_constraints_is_empty(tile_dest_same_copy_level)) {
                /* dcst1 is parameterized by the source src_copy_level iterators */
                *dcst1 = pluto_compute_region_data(dest, tile_dest_same_copy_level, read_acc,
                        src_copy_level, prog);
                pluto_constraints_subtract(tile_dest, tile_dest_same_copy_level);
            }
            pluto_constraints_free(copyleveleq);
            pluto_constraints_free(tile_dest_same_copy_level);
        }
    }else{ // some target iterations within the source tile
        // target and source statements definitely have outer src_copy_level surrounding loops common
        // so the target src_copy_level iterators can be projected out
        // !!!roshan not sure if projecting out target src_copy_level iterators is necessary
        // project out target src_copy_level loop iterators and source loop iterators inner to src_copy_level
        pluto_constraints_project_out(tile_dest_inside, src_copy_level,
                src->trans->nrows);
        pluto_constraints_project_out(tile_dest, src_copy_level,
                src->trans->nrows);

        // difference should be taken after the two sets are parameterized on the source src_copy_level iterators
        pluto_constraints_subtract(tile_dest, tile_dest_inside);

        dep->src_unique_dpolytope = pluto_constraints_dup(tile_dest);
        for(j=0; j<src->trans->nrows; j++)
            pluto_constraints_add_dim(dep->src_unique_dpolytope , src_copy_level, NULL);
        pluto_constraints_intersect(dep->src_unique_dpolytope, tdpoly);
    }

    /* dcst is parameterized by the source src_copy_level iterators */
    dcst = pluto_compute_region_data(dest, tile_dest, read_acc,
            src_copy_level, prog);

    IF_DEBUG(printf("Data flowing out for dep %d\n", dep->id+1););
    IF_DEBUG(pluto_constraints_print(stdout, dcst));
    IF_DEBUG(print_polylib_visual_sets("d0", dcst));

    pluto_constraints_free(tile_dest);
    pluto_constraints_free(tile_dest_inside);
    pluto_constraints_free(tdpoly);
    pluto_constraints_free(param_eq);
    return dcst;
}

/*
 * 'src_copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 *
 * Output format:  src_copy_level, acc->nrows, prog->npar + 1
 */
PlutoConstraints *compute_flow_out(struct stmt_access_pair *wacc_stmt, 
        int src_copy_level, int *copy_level, PlutoProg *prog, int *pi_mappings)
{
    int i;

    /* Stmt *stmt = wacc_stmt->stmt; */
    PlutoAccess *wacc = wacc_stmt->acc;

    PlutoConstraints *udcst = pluto_constraints_empty(src_copy_level
            +wacc->mat->nrows+prog->npar+1);

    IF_DEBUG(printf("Computing flow out set for %s\n", wacc->name););

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];

		if (options->data_dist) {
			/* Only RAW and RAR deps matter when data is distributed*/
			if (dep->type != OSL_DEPENDENCE_RAW && dep->type != OSL_DEPENDENCE_RAR) continue;
		}
		else {
			/* Only RAW deps matter */
			if (dep->type != OSL_DEPENDENCE_RAW) continue;
		}

        // if (dep->dirvec[src_copy_level+1] == DEP_ZERO) continue;

        /* If the dependence doesn't originate from this access */
        if (dep->src_acc != wacc) continue;

        assert(dep->dest_acc != NULL);

        PlutoConstraints *dcst = compute_flow_out_of_dep(dep, src_copy_level, copy_level, prog, 0, NULL, pi_mappings);

        IF_DEBUG(print_polylib_visual_sets("d1", dcst));

        pluto_constraints_unionize(udcst, dcst);
        pluto_constraints_free(dcst);
    }
    //pluto_constraints_print(stdout, udcst);
    return udcst;
}

void split_deps_acc_flowout(PlutoConstraintsList *atomic_flowouts, int copy_level, int access_nrows, PlutoProg *prog) {

	PlutoConstraintsList *curr = atomic_flowouts;

	//int n = check_for_disjoint_sets(atomic_flowouts);

	if(curr->constraints != NULL)
		return;

	while(curr != NULL && curr->constraints != NULL) {

		PlutoDepList *dep_list = curr->deps;
        PlutoDepList *prev_dep_list = NULL;
		PlutoConstraints *data = pluto_constraints_dup(curr->constraints);
		IF_DEBUG(print_polylib_visual_sets("dat", data));

		assert(data->ncols == copy_level + access_nrows + prog->npar +  1);

		//Create a identity access function for data with outer copy level parameters
		PlutoMatrix *data_access = pluto_matrix_alloc(access_nrows, data->ncols);
		pluto_matrix_set(data_access, 0);
		int i;
		for(i=0; i<access_nrows; i++){
			data_access->val[i][copy_level + i] = 1;
		}

		//Dep *d = pluto_dependence_dup(dep_list->dep, data);

		while(dep_list != NULL) {
			assert(dep_list->dep != NULL);
			Dep *dep = pluto_dependence_dup(dep_list->dep, dep_list->dep->dpolytope, dep_list->dep->src_unique_dpolytope);

			Stmt *dest = prog->stmts[dep->dest];
			Stmt *src = prog->stmts[dep->src];

			int src_ncols = src->trans->nrows;
			int dest_ncols = dest->trans->nrows;

			PlutoConstraints *src_iterators = pluto_constraints_dup(dep->src_unique_dpolytope);

			assert(src_iterators != NULL);

			if(pluto_constraints_is_empty(src_iterators)) {
				dep_list = dep_list->next;
				continue;
			}

			pluto_constraints_project_out(src_iterators, src_ncols, dest_ncols);
			IF_DEBUG(print_polylib_visual_sets("src", src_iterators));

			PlutoConstraints *tdpoly = pluto_find_dependence(src_iterators, data, dep, dep, prog, data_access);

            if (tdpoly!=NULL) {
                PlutoConstraints *parm_eq = get_context_equality(copy_level, src_ncols, tdpoly->ncols);

                pluto_constraints_intersect(tdpoly, parm_eq);

                assert(tdpoly != NULL);

    /*
                PlutoConstraints* intersect = pluto_constraints_dup(tdpoly);
                //pluto_constraints_project_out(intersect, copy_level+access_nrows, src_ncols);
                pluto_constraints_project_out(intersect, 0, src_ncols);

                IF_DEBUG(print_polylib_visual_sets("int", intersect));

                PlutoConstraints *diff = pluto_constraints_difference(intersect, data);
                assert(pluto_constraints_is_empty(diff));

                //IF_DEBUG(print_polylib_visual_sets("tdpoly", tdpoly));

                //pluto_constraints_project_out(tdpoly, 0, copy_level + access_nrows);

    */
                pluto_constraints_project_out(tdpoly, src_ncols, copy_level + access_nrows);


                IF_DEBUG(print_polylib_visual_sets("intersect", tdpoly));

                for(i=0; i<dest_ncols; i++)
                    pluto_constraints_add_dim(tdpoly, src_ncols, NULL);

                pluto_constraints_intersect(dep->src_unique_dpolytope, tdpoly);

                if(pluto_constraints_is_empty(tdpoly))
                    tdpoly = NULL;

            }

            if ((tdpoly==NULL) || pluto_constraints_is_empty(dep->src_unique_dpolytope)) {
                if (prev_dep_list == NULL) {
                    curr->deps = dep_list->next;
                    dep_list->next = NULL;
                    pluto_deps_list_free(dep_list);
                    dep_list = curr->deps;
                }
                else {
                    prev_dep_list->next = dep_list->next;
                    dep_list->next = NULL;
                    pluto_deps_list_free(dep_list);
                    dep_list = prev_dep_list->next;
                }
            } else {

                dep_list->dep = dep;

                prev_dep_list = dep_list;
                dep_list = dep_list->next;
            }
		}

		curr = curr->next;
	}
}

/*
 * Computes the flow out sets with dependence based spliting
 * 'src_copy_level' outer schedule rows will be treated as
 * parameters along with global parameters 
 *
 * Output format:  src_copy_level, acc->nrows, prog->npar + 1
 */
void compute_flow_out_partitions(struct stmt_access_pair *wacc_stmt,
        int src_copy_level, int *copy_level, PlutoProg *prog, PlutoConstraintsList *atomic_flowouts, int *pi_mappings)
{
    int i;

    PlutoAccess *wacc = wacc_stmt->acc;

    //PlutoConstraintsList *atomic_flowouts = pluto_constraints_list_alloc(NULL);
    // pluto_constraints_empty(src_copy_level +wacc->mat->nrows+prog->npar+1);

    IF_DEBUG(printf("Computing flow out set with dependence spliting for %s\n", wacc->name););

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

        // if (dep->dirvec[src_copy_level+1] == DEP_ZERO) continue;

        /* If the dependence doesn't originate from this access */
        if (dep->src_acc != wacc) continue;

        assert(dep->dest_acc != NULL);

        PlutoConstraints *dcst1 = NULL;
        PlutoConstraints *dcst = compute_flow_out_of_dep(dep, src_copy_level, copy_level, prog, 1, &dcst1, pi_mappings);

        if(!pluto_constraints_is_empty(dcst)){
            assert(!pluto_constraints_is_empty(dep->src_unique_dpolytope));
            /* split the dependence into atomic sections*/
            split_flow_out_set(prog, atomic_flowouts, dcst, dep);
        }

        if((dcst1 != NULL) && !pluto_constraints_is_empty(dcst1)){
            assert(!pluto_constraints_is_empty(dep->src_unique_dpolytope));
            /* split the dependence into atomic sections*/
            split_flow_out_set(prog, atomic_flowouts, dcst1, dep);
        }
    }
    return;
}


/*
 * Computes read-in set
 * racc_stmt: access pair for which first read needs to be computed
 * copy_level: number of (outer) schedule rows to be treated as
 * parameters along with global parameters
 *
 * Output format:  [copy_level, acc->nrows, prog->npar + 1]
 */
PlutoConstraints *compute_read_in(struct stmt_access_pair *racc_stmt,
        int copy_level, PlutoProg *prog)
{
    int i;

    Stmt *rstmt = racc_stmt->stmt;
    PlutoAccess *racc = racc_stmt->acc;

    PlutoConstraints *srcdomain = pluto_get_new_domain(rstmt);

    PlutoConstraints *urcst =  pluto_compute_region_data(rstmt,
            srcdomain, racc, copy_level, prog);

    return urcst;

    IF_DEBUG(printf("Computing read in set for %s\n", racc->name););

    IF_DEBUG(printf("Data read by the tile\n"););
    IF_DEBUG(pluto_constraints_print(stdout, urcst););

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];

        if (dep->type != OSL_DEPENDENCE_RAR) continue;

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* Only dependences with this one as dest access */
        if (racc != dep->dest_acc) continue;

        IF_DEBUG(printf("For dep %d\n", dep->id+1));

        const PlutoAccess *racc_src = dep->src_acc;

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //printf("tdpoly\n");
        //pluto_constraints_print(stdout, tdpoly);

        PlutoConstraints *param_eq;

        /* No need to check whether source and target share all of the
         * copy_level loops; handled naturally since those scalar dimensions
         * would be part of the 'copy_level' dimensions. If the statements are
         * separated somewhere, subtracting param_eq will not changed tdpoly
         */
        param_eq = get_context_equality(copy_level, src->trans->nrows,
                tdpoly->ncols);

        PlutoConstraints *tile_dest_outside = pluto_constraints_difference(tdpoly, param_eq);

        /* Parameteric in the outer source 'copy_level' dimensions (since
         * urcst is constructed that way) */
        pluto_constraints_project_out(tile_dest_outside, src->trans->nrows, dest->trans->nrows);

        IF_DEBUG(printf("Dep target iters outside of tile that write to same variable subseq to tile exec\n"
                    ););
        IF_DEBUG(pluto_constraints_print(stdout, tile_dest_outside););

        /* Values that'll be read outside the tile with source inside;
         * tile_dest_outside format [ src copy_level params | dest | par | 1]
         */
        PlutoConstraints *rcst =
            pluto_compute_region_data(src, tile_dest_outside, racc_src,
                    copy_level, prog);

        IF_DEBUG(printf("Values written outside for Dep %d with source inside tile\n", i+1););
        IF_DEBUG(pluto_constraints_print(stdout, rcst););

        urcst = pluto_constraints_subtract(urcst, rcst);

        pluto_constraints_free(tdpoly);
        pluto_constraints_free(rcst);
        pluto_constraints_free(tile_dest_outside);
        pluto_constraints_free(param_eq);
    }
    pluto_constraints_free(srcdomain);
    // IF_DEBUG(printf("Last write out\n"););
    // IF_DEBUG(pluto_constraints_print(stdout, uwcst););
    return urcst;
}

/*
 * Computes write-out set
 * wacc_stmt: access pair for which last writes need to be computed
 * copy_level: number of (outer) schedule rows to be treated as
 * parameters along with global parameters 
 *
 * NOTE: requires transitive WAW deps
 *
 * Output format:  [copy_level, acc->nrows, prog->npar + 1]
 */
PlutoConstraints *compute_write_out(struct stmt_access_pair *wacc_stmt, 
        int copy_level, PlutoProg *prog)
{
    int i;

    Stmt *wstmt = wacc_stmt->stmt;
    PlutoAccess *wacc = wacc_stmt->acc;
    Dep **deps;
    int ndeps;
    if (options->lastwriter) {
        // contains transitive dependences only due to WAR dependences
        deps = prog->transdeps;
        ndeps = prog->ntransdeps;
    }else{
        // contains transitive dependences
        deps = prog->deps;
        ndeps = prog->ndeps;
    }

    PlutoConstraints *srcdomain = pluto_get_new_domain(wstmt);

    /* Locations written to inside this tile */
    PlutoConstraints *uwcst =  pluto_compute_region_data(wstmt, 
            srcdomain, wacc, copy_level, prog);

    IF_DEBUG(printf("Computing write out set for %s\n", wacc->name););

    IF_DEBUG(printf("Data written to in tile\n"););
    IF_DEBUG(pluto_constraints_print(stdout, uwcst););

    // requires transitive WAR dependences if they exist
    for (i=0; i<ndeps; i++)   {
        Dep *dep = deps[i];
        // if (dep->type != OSL_DEPENDENCE_WAW && dep->type != OSL_DEPENDENCE_WAR) continue;
        if (dep->type != OSL_DEPENDENCE_WAW) continue;

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* Only dependences with this one as source access */
        if (dep->src_acc != wacc) continue;

        IF_DEBUG(printf("For dep %d\n", dep->id+1));

        const PlutoAccess *wacc_src = dep->src_acc;

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);

        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //printf("tdpoly\n");
        //pluto_constraints_print(stdout, tdpoly);

        PlutoConstraints *param_eq;

        /* No need to check whether source and target share all of the
         * copy_level loops; handled naturally since those scalar dimensions
         * would be part of the 'copy_level' dimensions. If the statements are 
         * separated somewhere, subtracting param_eq will not changed tdpoly
         */
        param_eq = get_context_equality(copy_level, src->trans->nrows, 
                tdpoly->ncols);

        /* Source/target iteration pairs with target lying outside the tile */
        PlutoConstraints *tile_dest_outside = pluto_constraints_difference(tdpoly, 
                param_eq);

        /* Parameteric in the outer source 'copy_level' dimensions (since
         * uwcst is constructed that way) */
        /* Yields iterations inside the tile whose write locations are again 
         * written to outside the tile */
        pluto_constraints_project_out(tile_dest_outside, src->trans->nrows,
                dest->trans->nrows);

        IF_DEBUG(printf("Dep target iters outside of tile that write to same variable subseq to tile exec\n"
                    ););
        IF_DEBUG(pluto_constraints_print(stdout, tile_dest_outside););

        /* Values that'll be written outside the tile with source inside; 
         * tile_dest_outside format [ src copy_level params | dest | par | 1]
         */
        /* Locations that will again be written to outside */
        PlutoConstraints *wcst = 
            pluto_compute_region_data(src, tile_dest_outside, wacc_src, 
                    copy_level, prog);

        IF_DEBUG(printf("Values written outside for Dep %d with source inside tile\n", i+1););
        IF_DEBUG(pluto_constraints_print(stdout, wcst););

        /* Subtracting out locations that will be written to outside */
        uwcst = pluto_constraints_subtract(uwcst, wcst);

        pluto_constraints_free(tdpoly);
        pluto_constraints_free(wcst);
        pluto_constraints_free(tile_dest_outside);
        pluto_constraints_free(param_eq);
    }
    pluto_constraints_free(srcdomain);
    // IF_DEBUG(printf("Last write out\n"););
    // IF_DEBUG(pluto_constraints_print(stdout, uwcst););
    return uwcst;
}
