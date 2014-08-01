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

#include "scoplib/statement.h"
#include "scoplib/access.h"

#include <isl/map.h>
#include <isl/mat.h>
#include <isl/set.h>
#include <isl/flow.h>
#include <isl/union_map.h>


/*
 * 'copy_level' outer schedule rows will be treated as
 * parameters along with global parameters; \pi(those copy_level rows,
 * global params, nprocs) is what is constructed
 *
 * loop_num: parallel loop number
 */
void generate_pi(FILE *outfp, int copy_level, const PlutoProg *prog, int loop_num,
        int *dimensions_to_skip, int num_dimensions_to_skip)
{
    int i;

    fprintf(outfp, "\
#include <math.h>\n\
#include <assert.h>\n\
#include \"pi_defs.h\"\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

    fprintf(outfp, "int pi_%d(", loop_num);
    assert(copy_level>0);
    int index_to_skip = 0;
    int index_offset = 1;
    index_offset++; // due to the scalar dimension that will be added to separate write-out statements from the rest
    for (i=0; i<copy_level; i++) {
        if ((index_to_skip >= num_dimensions_to_skip) || (i < dimensions_to_skip[index_to_skip])) {
            fprintf(outfp, " int t%d", i+index_offset);
            i++;
            break;
        } else {
            assert(i == dimensions_to_skip[index_to_skip]);
            index_to_skip++;
        }
    }
    for (; i<copy_level; i++) {
        if ((index_to_skip >= num_dimensions_to_skip) || (i < dimensions_to_skip[index_to_skip])) {
            fprintf(outfp, ", int t%d", i+index_offset);
        } else {
            assert(i == dimensions_to_skip[index_to_skip]);
            index_to_skip++;
        }
    }
    for (i=0; i<prog->npar; i++) {
        fprintf(outfp, ", int %s", prog->params[i]);
    }
    fprintf(outfp, ", int nprocs)\n{\n");

    copy_level++; // due to the scalar dimension that will be added to separate write-out statements from the rest
    fprintf(outfp, "\
            int __p, __lb, __ub;\n\
            long compute_start, compute_end;\n\
            __lb = %s%d;\n\
            __ub = %s%d;\n\
            long __n = __ub - __lb + 1;\n\
            for (__p=0; __p<nprocs; __p++)    {\n\
            if (__p < __n%%nprocs)  {\n\
            compute_start =  __lb + (__n/nprocs)*__p + __p;\n\
            /* procs with id < __n%%nprocs get an extra iteration to distributed\n\
             * the remainder */\n\
            compute_end = compute_start + (__n/nprocs)-1 + 1;\n\
            }else{\n\
            compute_start =  __lb + (__n/nprocs)*__p + __n%%nprocs;\n\
            compute_end = compute_start + (__n/nprocs) - 1;\n\
            }\n\
            if (t%d >= compute_start && t%d <= compute_end) return __p;\n\
            }\n\tassert(0);\n\treturn -1;\n}\n", "_LB_REPLACE_ME_DISTLOOG", loop_num,
            "_UB_REPLACE_ME_DISTLOOG", loop_num,
            copy_level, copy_level);
}

void pluto_gen_sigma_cloog_code(const PlutoProg *sigma, FILE *cloogfp, FILE *outfp, int partition_id) {
    int i;
    for(i = 0; i < sigma->nstmts; i++){
        sigma->stmts[i]->id += partition_id+1;
    }

    pluto_gen_cloog_file(cloogfp, sigma);

    for(i = 0; i < sigma->nstmts; i++){
        sigma->stmts[i]->id -= partition_id+1;
    }
    rewind(cloogfp);

    // pluto_gen_cloog_code(sigma, cloogfp, stdout);
    // rewind(cloogfp);

    generate_declarations(sigma, outfp);
    pluto_gen_cloog_code(sigma, 1, sigma->num_hyperplanes, cloogfp, outfp);
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
        int naccs, int *copy_level, PlutoProg *prog, PlutoConstraintsList *partitions, int loop_num, int *pi_mappings)
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
            PlutoConstraints *tdpoly = pluto_constraints_dup(dep->src_unique_dpolytype);
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

            for (i=0;i<src_copy_level; i++) {
                pluto_prog_add_hyperplane(sigma_check,0,H_LOOP);
            }
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
                        pluto_constraints_add_dim(tdpoly, 0);
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
                sigma_text = malloc(strlen("add_proc_to_receiver_list(pi_(, nprocs),my_rank)")
                        +5+strlen(indices)+1);
                sprintf(sigma_text, "add_proc_to_receiver_list(pi_%d(%s, nprocs),my_rank)",
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

                pluto_add_stmt(sigma,tdpoly,trans,iters,sigma_text, SIGMA);
                pluto_add_stmt(is_receiver,tdpoly,trans,iters,is_receiver_text, SIGMA);
                if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) pluto_add_stmt(sigma_check,tdpoly,trans,iters,sigma_check_text, SIGMA);

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
            pluto_detect_hyperplane_types(prog);
            //pluto_prog_print(sigma);
            if (sigma->nstmts >= 1) {
                assert(sigma->stmts[0]->trans->nrows == sigma->num_hyperplanes);
            }

            pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts, 0);
            pluto_separate_stmts(is_receiver, is_receiver->stmts, is_receiver->nstmts, 0);
            if (options->fop_unicast_runtime && !broadcast_of_partition[partition_id]) pluto_separate_stmts(sigma_check, sigma_check->stmts, sigma_check->nstmts, 0);
        }

        FILE *outfp, *cloogfp=NULL;
        outfp = fopen("sigma_fop.c", "a");
        assert(outfp != NULL);

        fprintf(outfp, "void sigma_%s_%d_%d(", acc_name, partition_id+1, loop_num);
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) fprintf(outfp, ", ");
            if (i<=src_copy_level-1) fprintf(outfp, "int ts%d", i+1);
            else fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");

        if (broadcast_of_partition[partition_id]) {
            fprintf(outfp, "int __p;\nfor (__p=0; __p<nprocs; __p++) add_proc_to_receiver_list(__p, my_rank);\n");
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
        for (i=0; i<src_copy_level+prog->npar; i++)    {
            if (i!=0) fprintf(outfp, ", ");
            if (i<=src_copy_level-1) fprintf(outfp, "int ts%d", i+1);
            else fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
        }
        fprintf(outfp, ", int my_rank, int nprocs)\n{\n");

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
            for (i=0; i<src_copy_level+prog->npar; i++)    {
                if (i!=0) fprintf(outfp, ", ");
                if (i<=src_copy_level-1) fprintf(outfp, "int ts%d", i+1);
                else fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
            }
            fprintf(outfp, ", int my_rank, int nprocs)\n{\n");
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

        fclose(outfp);
        curr = curr->next;
        partition_id++;
    }
    assert(partition_id == num_partitions);

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
        int naccs, int *copy_level, PlutoProg *prog, int loop_num, int *pi_mappings)
{
    int i, j;
    char *acc_name;
    int src_copy_level = copy_level[loop_num];

    assert(naccs != 0);

    acc_name = wacc_stmts[0]->acc->name;

    PlutoProg *sigma = pluto_prog_alloc();

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

    PlutoConstraints *tdpoly_list[prog->ndeps];
    int dep_loop_num_list[prog->ndeps];
    int num_dest_loops = 0;
    int broadcast = !options->commopt; // nocommopt => no precise determination of receivers => broadcast

    if (!broadcast) {
        for (i=0; i<prog->ndeps; i++)   {
            Dep *dep = prog->deps[i];
            /* Only RAW deps matter */
            if (dep->type != OSL_DEPENDENCE_RAW) continue;

            Stmt *src = prog->stmts[dep->src];
            Stmt *dest = prog->stmts[dep->dest];

            /* If the dependence doesn't originate from this access */
            for (j=0; j<naccs; j++) {
                if (src->writes[0] == wacc_stmts[j]->acc) break;
            }
            if (j==naccs)   continue;

            /* Add a statement for this dependence */
            int dep_loop_num = pi_mappings[dest->id];
            if (dep_loop_num == -1) { 
                // destination statement will be executed by all processors
                broadcast = 1;
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
        }
    }

    if (!broadcast) {
        for (i=0; i<num_dest_loops; i++) {
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
                    pluto_constraints_add_dim(tdpoly, 0);
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
            sigma_text = malloc(strlen("add_proc_to_receiver_list(pi_(, nprocs),my_rank)") 
                    + 5 + strlen(indices)+1);
            sprintf(sigma_text, "add_proc_to_receiver_list(pi_%d(%s, nprocs),my_rank)", 
                    dep_loop_num, indices);
            pluto_add_stmt(sigma,tdpoly,trans,iters,sigma_text, SIGMA);

            for (j=0; j<dest_copy_level; j++) {
                free(iters[j]);
            }

            // pluto_stmt_print(stdout, sigma->stmts[0]);

            free(sigma_text);
            free(indices);
            free(iters);
            pluto_matrix_free(trans);
            pluto_constraints_free(tdpoly);
        }

        //pluto_constraints_print(stdout, udcst);

        pluto_pad_stmt_transformations(sigma);
        pluto_detect_hyperplane_types(prog);
        //pluto_prog_print(sigma);
        if (sigma->nstmts >= 1) {
            assert(sigma->stmts[0]->trans->nrows == sigma->num_hyperplanes);
        }

        pluto_separate_stmts(sigma, sigma->stmts, sigma->nstmts, 0);
    }

    FILE *outfp, *cloogfp = NULL;

    outfp = fopen("sigma.c", "a");

    assert(outfp != NULL);

    fprintf(outfp, "void sigma_%s_%d(", acc_name, loop_num);
    for (i=0; i<src_copy_level+prog->npar; i++)    {
        if (i!=0) fprintf(outfp, ", ");
        if (i<=src_copy_level-1) fprintf(outfp, "int ts%d", i+1);
        else fprintf(outfp, "int %s", prog->params[i-src_copy_level]);
    }
    fprintf(outfp, ", int my_rank, int nprocs)\n{\n");

    if (broadcast) {
        fprintf(outfp, "int __p;\nfor (__p=0; __p<nprocs; __p++) add_proc_to_receiver_list(__p, my_rank);\n");
    }
    else {
        cloogfp = fopen("sigma.cloog", "w+");
        pluto_gen_cloog_file(cloogfp, sigma);
        rewind(cloogfp);
        // pluto_gen_cloog_code(sigma, cloogfp, stdout);
        // rewind(cloogfp);
        generate_declarations(sigma, outfp);
        pluto_gen_cloog_code(sigma, 1, sigma->num_hyperplanes, cloogfp, outfp);
        fclose(cloogfp);

        for (i=0; i<sigma->nstmts; i++) {
            fprintf(outfp, "#undef S%d\n", i+1);
        }
    }

    fprintf(outfp, "}\n\n");

    fclose(outfp);

    pluto_prog_free(sigma);
}

PlutoConstraints* get_receiver_tiles_of_dep(Dep *dep, 
        int src_copy_level, int dest_copy_level, PlutoProg *prog, int use_tile_dest_outside)
{
    int total_copy_level = src_copy_level + dest_copy_level;

    Stmt *src = prog->stmts[dep->src];
    Stmt *dest = prog->stmts[dep->dest];

    PlutoConstraints *tdpoly;
    if (!use_tile_dest_outside) tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
    else tdpoly = pluto_constraints_dup(dep->src_unique_dpolytype);
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
void generate_tau(struct stmt_access_pair *wacc_stmt, 
        int copy_level, PlutoProg *prog)
{
    int i, j;

    Stmt *stmt = wacc_stmt->stmt;
    PlutoAccess *wacc = wacc_stmt->acc;

    printf("Compute tau\n");
    PlutoProg *tau = pluto_prog_alloc();

    for (i=0; i<copy_level; i++) {
        char param[6];
        sprintf(param, "td%d",i+1);
        pluto_prog_add_param(tau, param, tau->npar);
    }

    for (i=0; i<prog->npar; i++) {
        pluto_prog_add_param(tau, prog->params[i], tau->npar);
    }

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

        Stmt *src = prog->stmts[dep->src];
        Stmt *dest = prog->stmts[dep->dest];

        /* If the dependence doesn't originate from this access */
        if (src->writes[0] != wacc) continue;

        assert(src->id == stmt->id);

        PlutoConstraints *tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
        assert(tdpoly->ncols == src->trans->nrows+dest->trans->nrows+prog->npar+1);

        //printf("tdpoly\n");
        //pluto_constraints_print(stdout, tdpoly);

        pluto_constraints_project_out(tdpoly, copy_level, src->trans->nrows-copy_level);
        pluto_constraints_project_out(tdpoly, 2*copy_level, dest->trans->nrows-copy_level);
        assert(tdpoly->ncols == 2*copy_level + prog->npar+1);

        IF_DEBUG(printf("Tile space dep poly (src, dest, param, 1)\n"););
        IF_DEBUG(pluto_constraints_print(stdout, tdpoly););

        char **iters = malloc(copy_level*sizeof(char *));
        char *indices = malloc(copy_level*8+1);
        strcpy(indices, "");
        for (j=0; j<copy_level; j++) {
            iters[j] = malloc(5);
            sprintf(iters[j], "d%d", j+1);
            sprintf(indices+strlen(indices), "[%s]", iters[j]);
        }
        // printf("%s\n", indices);

        PlutoMatrix *trans = pluto_matrix_identity(copy_level);
        for (j=0; j<copy_level+prog->npar; j++) {
            pluto_matrix_add_col(trans, trans->ncols);
        }
        /* constant part */
        pluto_matrix_add_col(trans, trans->ncols);
        char *tau_text = malloc(strlen("add_proc_to_sender_list(pi_table,my_rank)") + strlen(indices)+1);
        sprintf(tau_text, "add_proc_to_sender_list(pi_table%s,my_rank)",indices);
        pluto_add_stmt(tau,tdpoly,trans,iters,tau_text,TAU);

        // pluto_stmt_print(stdout, tau->stmts[0]);

        free(tau_text);
        pluto_matrix_free(trans);
        pluto_constraints_free(tdpoly);
    }
    //pluto_constraints_print(stdout, udcst);

    /* Separate code for different dependences */
    pluto_separate_stmts(tau, tau->stmts, tau->nstmts, 0);

    FILE *cloogfp = fopen("tau.cloog", "w+");
    pluto_gen_cloog_file(cloogfp, tau);
    rewind(cloogfp);
    FILE *outfp = fopen("tau.c", "w");

    fprintf(outfp, "\
#include <math.h>\n\
#include <assert.h>\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");


        fprintf(outfp, "extern int pi_table");
    for (i=0; i<copy_level; i++) {
        fprintf(outfp, "[%d]", PI_TABLE_SIZE);
    }
    fprintf(outfp, ";\n");

    fprintf(outfp, "void tau(");
    for (i=0; i<copy_level+prog->npar; i++)    {
        if (i!=0) fprintf(outfp, ", ");
        if (i<=copy_level-1) fprintf(outfp, "int td%d", i+1);
        else fprintf(outfp, "int %s", prog->params[i-copy_level]);
    }
    fprintf(outfp, ", int my_rank, int nprocs)\n{\n");

    // pluto_gen_cloog_code(tau, cloogfp, stdout);
    // rewind(cloogfp);
    generate_declarations(tau, outfp);
    pluto_gen_cloog_code(tau, 1, tau->num_hyperplanes, cloogfp, outfp);

    fprintf(outfp, "}\n\n");

    fclose(outfp);
    fclose(cloogfp);

    pluto_prog_free(tau);
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

Dep *pluto_dependence_dup(Dep *dep, PlutoConstraints *tdpoly, PlutoConstraints *tile_dest_outside){

	Dep *d = (Dep *)malloc(sizeof(Dep));

	d->depsat_poly = dep->depsat_poly;
	d->dest = dep->dest;
	d->dest_acc = dep->dest_acc;
	d->dirvec = dep->dirvec;
	d->dpolytope = pluto_constraints_dup(tdpoly);
	d->src_unique_dpolytype = pluto_constraints_dup(tile_dest_outside);
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
			pluto_constraints_add_dim(tdpoly2, 0);

		pluto_constraints_intersect(tdpoly2,curr_dep->dpolytope);
		dep2 = pluto_dependence_dup(curr_dep, tdpoly2, tdpoly2);

		pluto_deps_list_append(intersect_list, dep2);

		//compute the remaining source iterators from curr_dep
		diff2 = pluto_constraints_difference(td2, intersect2);

		for(i=0; i<td2_dest_ncols; i++)
			pluto_constraints_add_dim(diff2, td2_src_ncols);

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
				pluto_constraints_add_dim(tdpoly1, td1_src_ncols);

			pluto_constraints_intersect(tdpoly1,dep_to_split->dpolytope);
			dep1 = pluto_dependence_dup(dep_to_split, tdpoly1, tdpoly1);

			pluto_dep_list_add(curr->dep_list, dep1);

			for(i=0; i<td2_src_ncols; i++)
				pluto_constraints_add_dim(diff1, td1_src_ncols);

			//diff1 has the source iterators for next iteration
			pluto_constraints_intersect(dep_to_split->dpolytope,diff1);

			//TODO: Do the clean up
			curr = next;
			continue;

		}
		else if(!diff2_empty) {


			tdpoly1 = pluto_constraints_dup(intersect1);

			for(i=0; i<td1_dest_ncols; i++)
				pluto_constraints_add_dim(tdpoly1, td1_src_ncols);

			pluto_constraints_intersect(tdpoly1,dep_to_split->dpolytope);
			dep1 = pluto_dependence_dup(dep_to_split, tdpoly1, tdpoly1);

			tdpoly2 = pluto_constraints_dup(intersect2);
			for(i=0; i<td2_dest_ncols; i++)
				pluto_constraints_add_dim(tdpoly2, td2_src_ncols);

			pluto_constraints_intersect(tdpoly2,curr_dep->dpolytope);
			dep2 = pluto_dependence_dup(curr_dep, tdpoly2, tdpoly2);


			//Add the new dep list for the new intersected deps which access same data
			PlutoDepList *list = pluto_dep_list_alloc(dep1);
			pluto_deps_list_append(list, dep2);

			pluto_dep_list_list_add(curr, list);


			for(i=0; i<td2_dest_ncols; i++)
				pluto_constraints_add_dim(diff2, td2_src_ncols);

			pluto_constraints_intersect(curr_dep->dpolytope,diff2);

			//Recursively split the remaining dependences in the curr dep list
			intersect_deps_in_list(curr->dep_list->next, list, dep_to_split, prog);

			for(i=0; i<td1_dest_ncols; i++)
				pluto_constraints_add_dim(diff1, td1_src_ncols);

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
            pluto_constraints_add_dim(cst_intersect, cst_intersect->ncols-1);
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
        int copy_level, PlutoProg *prog, int use_tile_dest_outside)
{
    int j;
    Stmt *src = prog->stmts[dep->src];
    Stmt *dest = prog->stmts[dep->dest];
    PlutoAccess *write_acc = dep->src_acc;

    // printf("Write access\n");
    // pluto_matrix_print(stdout, write_acc->mat);

    PlutoConstraints *tdpoly;
    if (!use_tile_dest_outside) tdpoly = pluto_get_transformed_dpoly(dep, src, dest);
    else tdpoly = pluto_constraints_dup(dep->src_unique_dpolytype);
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
            pluto_constraints_add_dim(tile_src, 0);
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
            pluto_constraints_add_dim(tile_src_outside, 0);
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

        dep->src_unique_dpolytype = pluto_constraints_dup(tdpoly);
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

        dep->src_unique_dpolytype = pluto_constraints_dup(tile_dest);
        for(j=0; j<src->trans->nrows; j++)
            pluto_constraints_add_dim(dep->src_unique_dpolytype , src_copy_level);
        pluto_constraints_intersect(dep->src_unique_dpolytype, tdpoly);
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
        /* Only RAW deps matter */
        if (dep->type != OSL_DEPENDENCE_RAW) continue;

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

	if(curr->constraints == NULL)
		return;

	while(curr != NULL) {

		PlutoDepList *dep_list = curr->deps;
        PlutoDepList *prev_dep_list = NULL;
		PlutoConstraints *data = pluto_constraints_dup(curr->constraints);
		IF_DEBUG(print_polylib_visual_sets("dat", data));

		assert(data->ncols == copy_level + access_nrows + prog->npar +  1);

		//Create a identity access function for data with outer copy level parameters
		PlutoMatrix *data_access = pluto_matrix_alloc(access_nrows, data->ncols);
		pluto_matrix_initialize(data_access, 0);
		int i;
		for(i=0; i<access_nrows; i++){
			data_access->val[i][copy_level + i] = 1;
		}

		//Dep *d = pluto_dependence_dup(dep_list->dep, data);

		while(dep_list != NULL) {
			assert(dep_list->dep != NULL);
			Dep *dep = pluto_dependence_dup(dep_list->dep, dep_list->dep->dpolytope, dep_list->dep->src_unique_dpolytype);

			Stmt *dest = prog->stmts[dep->dest];
			Stmt *src = prog->stmts[dep->src];

			int src_ncols = src->trans->nrows;
			int dest_ncols = dest->trans->nrows;

			PlutoConstraints *src_iterators = pluto_constraints_dup(dep->src_unique_dpolytype);

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
                    pluto_constraints_add_dim(tdpoly, src_ncols);

                pluto_constraints_intersect(dep->src_unique_dpolytype, tdpoly);
            }

            if ((tdpoly==NULL) || pluto_constraints_is_empty(dep->src_unique_dpolytype)) {
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
            assert(!pluto_constraints_is_empty(dep->src_unique_dpolytype));
            /* split the dependence into atomic sections*/
            split_flow_out_set(prog, atomic_flowouts, dcst, dep);
        }

        if((dcst1 != NULL) && !pluto_constraints_is_empty(dcst1)){
            assert(!pluto_constraints_is_empty(dep->src_unique_dpolytype));
            /* split the dependence into atomic sections*/
            split_flow_out_set(prog, atomic_flowouts, dcst1, dep);
        }
    }
    return;
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
        if (wacc != src->writes[0]) continue;

        IF_DEBUG(printf("For dep %d\n", dep->id+1));

        const PlutoAccess *wacc_src = src->writes[0];

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
