/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
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
#include <assert.h>
#include <string.h>

#include "pluto.h"
#include "program.h"
#include "transforms.h"

#include "scoplib/access.h"

/* Get list of <Stmt, read acc> lists - each list corresponds to 
 * statements with read accesses to the same data; *num_data gives the number of 
 * data elements; nstmts_per_acc[i] is the number of stmts in list[i] 
 * Consider statements passed via 'stmts'
 * Return: num_data - number of variables, i.e., number of lists
 * nstmts_per_acc[]: number of stmt, acc pairs in each list
 * */
// not being used currently
struct stmt_access_pair ***get_read_access_with_stmts(Stmt **stmts, 
        int nstmts, int *num_data,
        int **nstmts_per_acc)
{
    int i, j, k, curr_num;
    int *num_stmts_per_acc;
    struct stmt_access_pair ***racc_stmts;

    curr_num = 0;
    num_stmts_per_acc = NULL;
    racc_stmts = NULL;

    for (i=0; i<nstmts; i++)  {
        Stmt *stmt = stmts[i];
        for (k=0; k<stmt->nreads; k++){
            struct stmt_access_pair *new = malloc(sizeof(struct stmt_access_pair));
            new->stmt = stmt;
            new->acc = stmt->reads[k];

            for (j=0; j<curr_num; j++)  {
                if (!strcmp(stmt->reads[k]->name, racc_stmts[j][0]->acc->name)) {
                    /* Add to end of array */
                    racc_stmts[j] = (struct stmt_access_pair **) realloc(racc_stmts[j], 
                        (num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
                    racc_stmts[j][num_stmts_per_acc[j]] = new;

                    num_stmts_per_acc[j]++;
                    break;
                }
            }
            if (j==curr_num)    {
                /* New data variable */
                racc_stmts = (struct stmt_access_pair ***)realloc(racc_stmts, 
                    (curr_num+1)*sizeof(struct stmt_access_pair **));
                racc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
                    sizeof(struct stmt_access_pair*));
                racc_stmts[curr_num][0] = new;

                num_stmts_per_acc = (int *)realloc(num_stmts_per_acc, 
                    (curr_num+1)*sizeof(int));
                num_stmts_per_acc[curr_num] = 1;
                curr_num++;
            }
        }
    }

    *num_data = curr_num;
    *nstmts_per_acc = num_stmts_per_acc;

    return racc_stmts;
}

/* Get list of <Stmt, write acc> lists - each list corresponds to 
 * statements with write accesses to the same data; *num_data gives the number of 
 * data elements; nstmts_per_acc[i] is the number of stmts in list[i] 
 * Consider statements passed via 'stmts'
 * Return: num_data - number of variables, i.e., number of lists
 * nstmts_per_acc[]: number of stmt, acc pairs in each list
 * */
struct stmt_access_pair ***get_write_access_with_stmts(Stmt **stmts, 
        int nstmts, int *num_data,
        int **nstmts_per_acc)
{
    int i, j, curr_num;
    int *num_stmts_per_acc;
    struct stmt_access_pair ***wacc_stmts;

    curr_num = 0;
    num_stmts_per_acc = NULL;
    wacc_stmts = NULL;

    for (i=0; i<nstmts; i++)  {
        Stmt *stmt = stmts[i];
        struct stmt_access_pair *new = malloc(sizeof(struct stmt_access_pair));
        new->stmt = stmt;
        new->acc = stmt->writes[0];

        for (j=0; j<curr_num; j++)  {
            if (!strcmp(stmt->writes[0]->name, wacc_stmts[j][0]->acc->name)) {
                /* Add to end of array */
                wacc_stmts[j] = (struct stmt_access_pair **) realloc(wacc_stmts[j], 
                    (num_stmts_per_acc[j]+1)*sizeof(struct stmt_access_pair*));
                wacc_stmts[j][num_stmts_per_acc[j]] = new;

                num_stmts_per_acc[j]++;
                break;
            }
        }
        if (j==curr_num)    {
            /* New data variable */
            wacc_stmts = (struct stmt_access_pair ***)realloc(wacc_stmts, 
                (curr_num+1)*sizeof(struct stmt_access_pair **));
            wacc_stmts[curr_num] = (struct stmt_access_pair **) malloc(
                sizeof(struct stmt_access_pair*));
            wacc_stmts[curr_num][0] = new;

            num_stmts_per_acc = (int *)realloc(num_stmts_per_acc, 
                (curr_num+1)*sizeof(int));
            num_stmts_per_acc[curr_num] = 1;
            curr_num++;
        }
    }

    *num_data = curr_num;
    *nstmts_per_acc = num_stmts_per_acc;

    return wacc_stmts;
}


static void pluto_mark_statements(PlutoProg *prog)
{
    int i;

    pluto_prog_add_param(prog, "my_rank", prog->npar);
    pluto_constraints_add_lb(prog->context, prog->npar-1, 0);

    for (i=0; i<prog->nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        if (stmt->type == LW_COPY_IN) {
            /* Add my_rank == 0; write-out set is copied-back only by master proc */
            pluto_constraints_set_var(stmt->domain, stmt->domain->ncols-2, 0);
        }
    }
}


/* Get access string */
char *reconstruct_access(PlutoAccess *acc)
{
    int ndims, i;
    ndims = acc->mat->nrows;

    char *access;

    access = (char *) malloc(strlen(acc->name)+ndims*(2+4)+1);

    strcpy(access, acc->name);

    for (i=0; i<ndims; i++) {
        strcat(access, "[");
        char tmp[5];
        sprintf(tmp, "d%d", i+1);
        strcat(access, tmp);
        strcat(access, "]");
    }

    return access;
}


void free_char_array_buffers(char** buffer, int size) {

    int i = 0;
    for(i = 0; i < size; ++i) {
        free(buffer[i]);
    }

    free(buffer);

    return;
}

/* 
 * Optimized communication code generation for the write-out set
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 * This function can be called by any flow-based communication scheme
 * caller is responsible for freeing write_out_stmts returned
 */  
Stmt **gen_write_out_code(struct stmt_access_pair **wacc_stmts, int num_accs,
        PlutoProg *prog, int *copy_level, int loop_num)
{
    int i, src_copy_level, acc_nrows;
    src_copy_level = copy_level[loop_num];

    assert(num_accs >= 1);
    char *access = reconstruct_access(wacc_stmts[0]->acc);
    acc_nrows = wacc_stmts[0]->acc->mat->nrows;
    Stmt *anchor_stmt = wacc_stmts[0]->stmt;
    char *acc_name = anchor_stmt->writes[0]->name;

    IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

    /* To be inside a loop: can't foresee other use */
    if (src_copy_level >= 1)    {
        assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));
    }

    PlutoConstraints *write_out = NULL;
    for (i=0; i<num_accs; i++)  {
        PlutoConstraints *write_out_one;
        //write_out_one = compute_last_writes(wacc_stmts[k], src_copy_level, prog);
        write_out_one = compute_write_out(wacc_stmts[i], src_copy_level, prog);
        if (write_out == NULL) write_out = pluto_constraints_dup(write_out_one);
        else write_out = pluto_constraints_unionize(write_out, write_out_one);
        pluto_constraints_free(write_out_one);
    }

    IF_DEBUG(printf("Write out set\n"););
    IF_DEBUG(pluto_constraints_print(stdout, write_out));

    // if opencl is not enabled, generate the MPI code for distributed cluster
    //
#ifdef PLUTO_OPENCL
    if(!options->opencl) {
#endif

// ===== from here on starts the code generation =====================================================

    /***************************************************************************************************/

    char *lw_buf_size = get_parametric_bounding_box(write_out, src_copy_level, acc_nrows, prog->npar, (const char **)prog->params);
    char *lw_recv_buf_size = malloc(1280);
    strcpy(lw_recv_buf_size, lw_buf_size);

    if (src_copy_level >= 1) {
        PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
        char *extent;
        /* Just the first one in flow_out is enough (rest all should give the
         * same since they are all under the same parallel loop and each
         * iteration of the parallel loop writes to distinct data) */
        get_parametric_extent_const(anchor_stmt_new_dom, src_copy_level-1, prog->npar,
                (const char **)prog->params, &extent);
        sprintf(lw_buf_size+strlen(lw_buf_size), "*ceilf((%s)/(float)nprocs)",
            extent);
        sprintf(lw_recv_buf_size+strlen(lw_recv_buf_size), "*(%s + nprocs)",
            extent);
        free(extent);
        pluto_constraints_free(anchor_stmt_new_dom);
    }
    IF_DEBUG(printf("Last writer buffer size for %s: %s\n", acc_name, lw_buf_size););

    if (src_copy_level >= 1)    {
        assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
                prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);
    }

    char *lwbufname = concat("lw_buf_", acc_name);
    char *lw_count_name = concat("lw_count_", acc_name);

    sprintf(prog->decls+strlen(prog->decls), "%s\
            = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &lw_buf_size_%s);\n", 
            lwbufname, lwbufname, lw_buf_size, acc_name);

    char *lw_holder_text = strdup("/* lw_holder */");

    char *lw_stmt_text = malloc(strlen(lwbufname) + strlen("[")+
        strlen(lw_count_name) + strlen("++] = ") + strlen(access) + 1);
    sprintf(lw_stmt_text, "%s[%s++] = %s", lwbufname, lw_count_name, access);

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *lw_holder = create_helper_stmt(anchor_stmt, src_copy_level, lw_holder_text, LW_COPY_OUT);
    Stmt *lw_copy_stmt = create_helper_stmt(anchor_stmt, src_copy_level, lw_stmt_text, LW_COPY_OUT);
    //IF_DEBUG(pluto_stmt_print(stdout, flow_copy_stmt););
    //IF_DEBUG(pluto_stmt_print(stdout, lw_copy_stmt););

    /* Add dimensions that actually scan the data space */
    for (i=0; i < acc_nrows; i++) {
        char iter[5];
        sprintf(iter, "d%d", i+1);
        pluto_stmt_add_dim(lw_copy_stmt, lw_copy_stmt->dim, lw_copy_stmt->trans->nrows, iter, H_LOOP, prog);
    }

    /* Add constraints on copy loops - completes domain of copy stmt */
    pluto_constraints_intersect(lw_copy_stmt->domain, write_out);

    char *lwrecvbufname = concat("lw_recv_buf_", acc_name);
    char *lw_displsname = concat("displs_lw_", acc_name);
    char *lw_recv_counts_name = concat("lw_recv_counts_", acc_name);

    char *comm_text_w = malloc(1024);

    if (options->commreport)
        sprintf(comm_text_w,  "IF_TIME(t_writeout_start = rtclock());");
    else
        strcpy(comm_text_w, "");
    sprintf(comm_text_w+strlen(comm_text_w), "\
        MPI_Gather(&%s, 1, MPI_INT,\
        %s, 1, MPI_INT, 0, MPI_COMM_WORLD);", lw_count_name, lw_recv_counts_name);
    sprintf(comm_text_w+strlen(comm_text_w), "MPI_Gatherv(%s, %s, MPI_DOUBLE,\
        %s, %s, %s, MPI_DOUBLE, 0, MPI_COMM_WORLD); %s = 0; lw_prev_proc=-1;",
            lwbufname, lw_count_name, lwrecvbufname, lw_recv_counts_name,
            lw_displsname, lw_count_name);
    if (options->commreport)
        sprintf(comm_text_w+strlen(comm_text_w),  "IF_TIME(t_writeout += rtclock() - t_writeout_start);");

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *comm_stmt_w = create_helper_stmt(anchor_stmt, src_copy_level-1, comm_text_w, LW_COMM_CALL);

    sprintf(prog->decls+strlen(prog->decls), 
            "%s = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &lw_recv_buf_size_%s);\n",
            lwrecvbufname, lwrecvbufname, lw_recv_buf_size, acc_name);

    sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
        %s[__p] = __p*(int)ceilf((%s)/(float)nprocs);}\n\n", lw_displsname,
        lw_recv_buf_size);
    char *lw_copyback_text = malloc(strlen(access) + strlen(" = lw_recv_buf_")
        + strlen(acc_name) + strlen("[") + strlen(lw_displsname) +
        strlen("[proc]+ count++];") + 1);
    sprintf(lw_copyback_text, "%s = lw_recv_buf_%s[%s[proc]+ count++];",
        access, acc_name, lw_displsname);

    char lw_proc_stmt_text[1024];
    sprintf(lw_proc_stmt_text, "proc = pi_%d(", loop_num);
    for (i=0; i<src_copy_level+prog->npar; i++)    {
        if (i>=1)  sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), ", ");
        if (i<=src_copy_level-1) sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), "t%d", i+1);
        else sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), "%s", prog->params[i-src_copy_level]);
    }
    sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text), ", nprocs)");
    sprintf(lw_proc_stmt_text+strlen(lw_proc_stmt_text),
        "; if (proc != lw_prev_proc) {lw_prev_proc = proc; count=0;}\
        if (lw_recv_counts_%s[proc] == 0) continue", acc_name);

    Stmt *lw_copyback_stmt = create_helper_stmt(anchor_stmt, src_copy_level, lw_copyback_text, LW_COPY_IN);
    Stmt *lw_copyback_proc_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
        lw_proc_stmt_text, LW_COPY_IN);

    /* Add dimensions that actually scan the data space */
    for (i=0; i < acc_nrows; i++) {
        char iter[5];
        sprintf(iter, "d%d", i+1);
        pluto_stmt_add_dim(lw_copyback_stmt, lw_copyback_stmt->dim,	lw_copyback_stmt->trans->nrows, iter, H_LOOP, prog);
    }

    assert(lw_copyback_stmt->domain->ncols == write_out->ncols);

    /* Add constraints for write-out set copy loops */
    lw_copyback_stmt->domain = pluto_constraints_intersect(lw_copyback_stmt->domain, write_out);

    // printf("lw_copyback_stmt before var\n");
    // IF_DEBUG(pluto_stmt_print(stdout, lw_copyback_stmt););

    // printf("lw_copyback_stmt\n");
    // IF_DEBUG(pluto_stmt_print(stdout, lw_copyback_stmt););
    /***************************************************************************************************/

    // write out stmts organized as follows (used later during code generation)
    // pack-guard stmt, pack stmt, comm stmt, unpack-guard stmt, unpack stmt
    Stmt **write_out_stmts = (Stmt **) malloc(5*sizeof(Stmt *));
    write_out_stmts[0] = lw_holder;
    write_out_stmts[1] = lw_copy_stmt;
    write_out_stmts[2] = comm_stmt_w;
    write_out_stmts[3] = lw_copyback_proc_stmt;
    write_out_stmts[4] = lw_copyback_stmt;

    pluto_constraints_free(write_out);
    free(lw_stmt_text);
    free(lw_holder_text);
    free(lw_copyback_text);
    free(lw_buf_size);
    free(lw_recv_buf_size);
    free(lw_count_name);
    free(lwbufname);
    free(lwrecvbufname);
    free(lw_displsname);
    free(lw_recv_counts_name);
    free(comm_text_w);
    free(access);

    return write_out_stmts;
#ifdef PLUTO_OPENCL
    }
    else{
        pluto_opencl_codegen(); // just dummy for now
        return NULL;
    }
#endif

}

void generate_pack_or_unpack(FILE *packfp, PlutoProg *prog, 
        PlutoConstraints *constraints, char *stmttext, PlutoStmtType stmttype, char **iters, int src_copy_level, int acc_nrows)
{
    FILE *packcloogfp;
    PlutoProg *pack;
    PlutoMatrix *packtrans;
    int i, total_level = src_copy_level + acc_nrows;
    PlutoConstraints *tdpoly;
    
    pack = pluto_prog_alloc();

    for (i=0; i<src_copy_level; i++) {
        char param[6];
        sprintf(param, "ts%d",i+1);
        pluto_prog_add_param(pack, param, pack->npar);
    }

    for (i=0; i<prog->npar; i++) {
        pluto_prog_add_param(pack, prog->params[i], pack->npar);
    }

    for (i=0;i<src_copy_level; i++) {
        pluto_prog_add_hyperplane(pack,0,H_LOOP);
    }

    tdpoly = pluto_constraints_dup(constraints);
    // move acc_nrows loop iterators to the beginning/top
    for (i=0; i<acc_nrows; i++) {
        pluto_constraints_add_dim(tdpoly, 0);
    }
    for (i=0; i<acc_nrows; i++) {
        pluto_constraints_interchange_cols(tdpoly, i, i+total_level);
    }
    for (i=0; i<acc_nrows; i++) {
        pluto_constraints_remove_dim(tdpoly, total_level);
    }

    packtrans = pluto_matrix_identity(acc_nrows);
    for (i=0; i<src_copy_level+prog->npar+1; i++) {
        pluto_matrix_add_col(packtrans, packtrans->ncols);
    }

    pluto_add_stmt(pack, tdpoly, packtrans, iters, stmttext, stmttype);

    packcloogfp = fopen("packunpack.cloog", "w+");
    assert(packcloogfp != NULL);
    pluto_gen_cloog_file(packcloogfp, pack);
    rewind(packcloogfp);

    generate_declarations(pack, packfp);
    pluto_gen_cloog_code(pack, 1, pack->num_hyperplanes, packcloogfp, packfp);

    pluto_matrix_free(packtrans);
    pluto_constraints_free(tdpoly);
    fclose(packcloogfp);
    pluto_prog_free(pack);
}

/* 
 * Optimized communication code generation wih dependence spliting
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop should be the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 */  
Stmt **gen_comm_code_opt_fop(struct stmt_access_pair **wacc_stmts, int num_accs, int nloops,
        PlutoProg *prog, int *copy_level, int loop_num, int *pi_mappings, int* num_comm_stmts)
{
    int i, k, src_copy_level, acc_nrows;
    src_copy_level = copy_level[loop_num];
    assert(src_copy_level>=1);

    assert(num_accs >= 1);
    char *access = reconstruct_access(wacc_stmts[0]->acc);
    acc_nrows = wacc_stmts[0]->acc->mat->nrows;
    Stmt *anchor_stmt = wacc_stmts[0]->stmt;
    char *acc_name = anchor_stmt->writes[0]->name;

    IF_DEBUG(printf("%d accesses: %s\n", num_accs, acc_name););

    /* To be inside a loop: can't foresee other use */
    assert((prog->hProps[src_copy_level-1].type == H_LOOP) || (prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP));

    /* Sender-side copying */
    PlutoConstraintsList *atomic_flowouts = pluto_constraints_list_alloc((PlutoConstraints *) NULL);
    PlutoConstraints *flow_out = NULL; // used only for get parametric extent/bounding box
    for (k=0; k<num_accs; k++)  {
        PlutoConstraints *flow_out_one;
		flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
		if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
		else{
			pluto_constraints_unionize(flow_out, flow_out_one);
		}
		pluto_constraints_free(flow_out_one);

        compute_flow_out_partitions(wacc_stmts[k], src_copy_level, copy_level, prog, atomic_flowouts, pi_mappings);
    }

    PlutoConstraintsList *prev = NULL, *curr = atomic_flowouts;
    PlutoConstraintsList *const_cst = NULL;
    while(curr != NULL) {
        PlutoConstraints *cst = curr->constraints;
        int const_bounds = 1;
        while (cst != NULL) {
            for (i=src_copy_level; i<src_copy_level+acc_nrows; i++) {
                if (get_const_bound_difference(cst, i) == -1) {
                    const_bounds = 0;
                    break;
                }
            }
            if (const_bounds == 0) break;
            cst = cst->next;
        }
        if (const_bounds == 1) {
            if (const_cst == NULL) {
                const_cst = curr;
            } 
            else {
                pluto_constraints_unionize(const_cst->constraints, curr->constraints);
                PlutoDepList *dep_list = curr->deps;
                while (dep_list != NULL) {
                    pluto_deps_list_append(const_cst->deps, dep_list->dep);
                    dep_list = dep_list->next;
                }

                prev->next = curr->next;
                curr->next = NULL;
                pluto_constraints_list_free(curr);
                curr = prev->next;
                continue;
            }
        }
        prev = curr;
        curr = curr->next;
    }

    split_deps_acc_flowout(atomic_flowouts, src_copy_level , acc_nrows, prog);

    generate_sigma_dep_split(wacc_stmts, num_accs, copy_level, prog, atomic_flowouts, loop_num, pi_mappings);

    // if opencl is not enabled, generate the MPI code for distributed cluster
    //
#ifdef PLUTO_OPENCL
    if(!options->opencl) {
#endif

// ===== from here on starts the code generation =====================================================

    /***************************************************************************************************/

    char *send_buf_size = 
        get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows, 
                prog->npar, (const char **)prog->params);
    char *recv_buf_size = malloc(1024 * 8);
    char *displs_size = malloc(1024 * 8);
    strcpy(recv_buf_size, send_buf_size);
    strcpy(displs_size, send_buf_size);

    PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
    char *extent;
    /* Just the first one in flow_out is enough (rest all should give the
     * same since they are all under the same parallel loop and each
     * iteration of the parallel loop writes to distinct data) */
    get_parametric_extent_const(anchor_stmt_new_dom, src_copy_level-1, prog->npar,
            (const char **)prog->params, &extent);
    sprintf(send_buf_size+strlen(send_buf_size), 
            "*ceilf((%s)/(float)nprocs)", extent);
    /* The + nprocs is needed since the displacement has to be set larger
     * when some processors have more iterations - worst case when the
     * first processor has an extra iteration; one then needs
     * (nprocs-1)*num_values_per_iteration additional space) */
    sprintf(recv_buf_size+strlen(recv_buf_size), "*(%s + nprocs)", 
            extent);
    sprintf(displs_size+strlen(displs_size), "*ceilf((%s)/(float)nprocs)",
            extent);
    free(extent);
    pluto_constraints_free(anchor_stmt_new_dom);
    IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

    char *sendbufname = concat("send_buf_", acc_name);
    char *send_counts_name = concat("send_counts_", acc_name);

    sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
            %s[__p] = 0;\n}\n", send_counts_name);
    sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
            %s[__p]\
            = (double *) polyrt_max_alloc(%s[__p], sizeof(double)*(%s), &send_buf_size_%s[__p]);\n}", 
            sendbufname, sendbufname, send_buf_size, acc_name);

    char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
            strlen(send_counts_name) + strlen("++] = ") + strlen(access) + 1);
    sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_counts_name, access);

    char *displsname = concat("displs_", acc_name);
    char *currdisplsname = concat("curr_displs_", acc_name);
    char *recvbufname = concat("recv_buf_", acc_name);
    char *recv_counts_name = concat("recv_counts_", acc_name);

    IF_DEBUG(printf("Recv buffer size for %s: %s\n", acc_name, recv_buf_size););

    sprintf(prog->decls+strlen(prog->decls), 
            "%s = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &recv_buf_size_%s);\n",
            recvbufname, recvbufname, recv_buf_size, acc_name);
    sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
            %s[__p] = __p*(%s);}\n\n", displsname, displs_size);

    char *comm_text = malloc(2048);

    if (options->commreport)
        sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
    else
        strcpy(comm_text, "");
    sprintf(comm_text+strlen(comm_text), "\
            MPI_Alltoall(%s, 1, MPI_INT,\
                %s, 1, MPI_INT, MPI_COMM_WORLD);", 
            send_counts_name, recv_counts_name);
    sprintf(comm_text+strlen(comm_text),  "\
            req_count=0;\
            for (__p=0; __p<nprocs; __p++) {\
            if (%s[__p] >= 1) {",
            send_counts_name);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
    sprintf(comm_text+strlen(comm_text),  "\
            MPI_Isend(%s[__p], %s[__p], MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}", 
            sendbufname, send_counts_name);
    sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
            if(%s[__p] >= 1) {\
            MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
            MPI_Waitall(req_count, reqs, stats);\
            for (__p=0; __p<nprocs; __p++) {\
            %s[__p] = 0;}",
            recv_counts_name, recvbufname, displsname, recv_counts_name,
            send_counts_name);
    sprintf(comm_text+strlen(comm_text), "for (__p=0; __p<nprocs; __p++) {\
            %s[__p] = %s[__p]; }", currdisplsname, displsname);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");

    /* Receiver-side copy */
    /* With send/recv -based more exact communication scheme, recv side copy
     * does not include all the sender-side copies; processes that do not send
     * data are to be skipped */
    assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
            prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

    char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
            + strlen("[") + strlen(currdisplsname) + strlen("++]") + 1);
    sprintf(flow_copyback_text, "%s = %s[%s++]",
            access, recvbufname, currdisplsname);

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *comm_stmt = create_helper_stmt(anchor_stmt, src_copy_level-1, comm_text, COMM_CALL);

    FILE *packfp = fopen("packunpack.c", "a");
    assert(packfp != NULL);

    char args[1024];
    strcpy(args,"");
    sprintf(args+strlen(args), "t%d", 1);
    /* make it src_copy_level+1 since we use an extra dimension to separate
     * statements */
    for (i=1; i<src_copy_level; i++) {
        sprintf(args+strlen(args), ",t%d", i+1);
    }

    char passed_args[1024];
    strcpy(passed_args,"");
    sprintf(passed_args+strlen(passed_args), "ts%d", 1);
    for (i=1; i<src_copy_level; i++) {
        sprintf(passed_args+strlen(passed_args), ",ts%d", i+1);
    }

    char decl_args[1024];
    strcpy(decl_args,"");
    sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
    for (i=1; i<src_copy_level; i++) {
        sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
    }

    char params[1024];
    strcpy(params, "");
    if (prog->npar>=1) {
        sprintf(params+strlen(params), "%s", prog->params[0]);
        for (i=1; i<prog->npar; i++) {
            sprintf(params+strlen(params), ",%s", prog->params[i]);
        }
    }

    char **iters;
    iters = malloc(acc_nrows * sizeof(char *));
    for (i=0; i < acc_nrows; i++) {
        iters[i] = malloc(5);
        sprintf(iters[i], "d%d", i+1);
    }

    char *pack_stmt_text = malloc(8192);
    strcpy(pack_stmt_text, "");
    char *unpack_stmt_text = malloc(8192);
    sprintf(unpack_stmt_text, "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "if ((my_rank != proc) && (%s[proc] > 0)) { ", recv_counts_name);
    k = 0;
    curr = atomic_flowouts;
    while(curr != NULL) {
        if ((curr->constraints == NULL) || pluto_constraints_is_empty(curr->constraints)) {
            curr = curr->next;
            continue;
        }
        
        int j, l;
        int broadcast = 0;

        char acc_name_k[512];
        sprintf(acc_name_k, "%s_%d", acc_name, k+1);

        if (options->fop_unicast_runtime) {
            PlutoConstraints *foifi_sets[nloops];
            PlutoConstraints *receivers[nloops];
            PlutoDepList *dep_list = curr->deps;
            for (l=0; l<nloops; l++) {
                foifi_sets[l] = NULL;
                receivers[l] = NULL;
            }
            while (dep_list != NULL) {
                Dep *dep = dep_list->dep;
                Stmt *dest = prog->stmts[dep->dest];
                int dependent_loop = pi_mappings[dest->id];
                if (dependent_loop == -1) { 
                    // destination statement will be executed by all processors
                    broadcast = 1;
                    break;
                }

                int dest_copy_level = copy_level[dependent_loop];

                PlutoConstraints *fi = compute_flow_in_of_dep(dep, dest_copy_level, prog, 1);
                PlutoConstraints *rt = get_receiver_tiles_of_dep(dep, src_copy_level, dest_copy_level, prog, 1);

                if (foifi_sets[dependent_loop] == NULL) 
                    foifi_sets[dependent_loop] = pluto_constraints_dup(fi);
                else
                    pluto_constraints_unionize(foifi_sets[dependent_loop], fi);

                if (receivers[dependent_loop] == NULL) 
                    receivers[dependent_loop] = pluto_constraints_dup(rt);
                else
                    pluto_constraints_unionize(receivers[dependent_loop], rt);

                pluto_constraints_free(fi);
                pluto_constraints_free(rt);
                dep_list = dep_list->next;
            }
            if (!broadcast) {
                char sigma_check_text[1024];
                sprintf(sigma_check_text, "distinct_recv = sigma_check_%s_%d(%s,%s", acc_name_k, loop_num, args, params);

                sprintf(pack_stmt_text+strlen(pack_stmt_text), "%s, my_rank, nprocs); if (distinct_recv == 1) {", sigma_check_text);
                sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "%s, proc, nprocs); if (distinct_recv == 1) {", sigma_check_text);

                for (l=0; l<nloops; l++) {
                    PlutoConstraints *foifi = foifi_sets[l];
                    if (foifi == NULL) continue;
                    
                    int dest_copy_level = copy_level[l];
                    assert(dest_copy_level>=1);
                    int total_copy_level = src_copy_level + dest_copy_level;
                    assert(foifi->ncols == dest_copy_level + acc_nrows + prog->npar + 1);

                    PlutoConstraints *fo = pluto_constraints_dup(curr->constraints);
                    // add target tile iterators for flow-out set
                    for (j=0; j<dest_copy_level; j++) {
                        pluto_constraints_add_dim(fo,src_copy_level);
                    }

                    // add source tile iterators for flow-in set
                    for (j=0; j<src_copy_level; j++) {
                        pluto_constraints_add_dim(foifi,0);
                    }

                    pluto_constraints_intersect(foifi, fo);

                    if (pluto_constraints_is_empty(foifi)) {
                        pluto_constraints_free(fo);
                        pluto_constraints_free(foifi);
                        pluto_constraints_free(receivers[l]);
                        continue;
                    }

                    PlutoConstraints *receiver_tiles = receivers[l];
                    assert(receiver_tiles->ncols == total_copy_level + prog->npar + 1);
                    assert(!pluto_constraints_is_empty(receiver_tiles));

                    char dest_args[1024];
                    strcpy(dest_args,"");
                    sprintf(dest_args+strlen(dest_args), "t%d", src_copy_level+1);
                    for (i=src_copy_level+1; i<total_copy_level; i++) {
                        sprintf(dest_args+strlen(dest_args), ",t%d", i+1);
                    }

                    char decl_dest_args[1024];
                    strcpy(decl_dest_args,"");
                    sprintf(decl_dest_args+strlen(decl_dest_args), "int ts%d", src_copy_level+1);
                    for (i=src_copy_level+1; i<total_copy_level; i++) {
                        sprintf(decl_dest_args+strlen(decl_dest_args), ",int ts%d", i+1);
                    }

                    char **dest_iters;
                    dest_iters = malloc(dest_copy_level * sizeof(char *));
                    for (i=0; i<dest_copy_level; i++) {
                        dest_iters[i] = malloc(5);
                        sprintf(dest_iters[i], "t%d", i+src_copy_level+1);
                    }

                    char pi_text[1024];
                    sprintf(pi_text, "recv_proc = pi_%d(%s,%s, nprocs)", l, dest_args, params);

                    sprintf(pack_stmt_text+strlen(pack_stmt_text), "pack_foifi_%s_%d_%d(%s,%s,%s, my_rank, nprocs", 
                            acc_name_k, loop_num, l, args, sendbufname, send_counts_name);
                    if (options->variables_not_global) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
                    sprintf(pack_stmt_text+strlen(pack_stmt_text), ");");

                    char flow_cg_text[1024];
                    sprintf(flow_cg_text, "%s; if (recv_proc != my_rank) \
                            { %s[recv_proc] = pack_foifi_recv_%s_%d_%d(%s,%s,%s[recv_proc],%s[recv_proc]", 
                            pi_text, send_counts_name, acc_name_k, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
                    if (options->variables_not_global) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);
                    sprintf(flow_cg_text+strlen(flow_cg_text), "); }");
                    
                    fprintf(packfp, "int pack_foifi_%s_%d_%d(%s,double **%s,int *%s, int my_rank, int nprocs", 
                            acc_name_k, loop_num, l, decl_args, sendbufname, send_counts_name);
                    if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
                    fprintf(packfp,"){\n");
                    fprintf(packfp, "\nint recv_proc;\n");
                    generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_cg_text, COPY_OUT, dest_iters, src_copy_level, dest_copy_level);
                    fprintf(packfp, "\nreturn 0;\n#undef S1;\n}\n\n");
                    
                    fprintf(packfp, "int pack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s", 
                            acc_name_k, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
                    if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
                    fprintf(packfp,"){\n");
                    generate_pack_or_unpack(packfp, prog, foifi, flow_copy_text, COPY_OUT, iters, total_copy_level, acc_nrows);
                    fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", send_counts_name);

                    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "%s[proc] = unpack_foifi_%s_%d_%d(%s,%s,%s[proc], my_rank, nprocs", 
                            currdisplsname, acc_name_k, loop_num, l, args, recvbufname, currdisplsname);
                    if (options->variables_not_global) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);
                    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

                    char flow_copyback_guard_text[1024];
                    sprintf(flow_copyback_guard_text, "%s; if (recv_proc == my_rank) \
                            { %s = unpack_foifi_recv_%s_%d_%d(%s,%s,%s,%s", 
                            pi_text, currdisplsname, acc_name_k, loop_num, l, passed_args, dest_args, recvbufname, currdisplsname);
                    if (options->variables_not_global) sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), ",%s", acc_name);
                    sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), "); }");
                    
                    fprintf(packfp, "int unpack_foifi_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs", 
                            acc_name_k, loop_num, l, decl_args, recvbufname, currdisplsname);
                    if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
                    fprintf(packfp,"){\n");
                    fprintf(packfp, "\nint recv_proc;\n");
                    generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_copyback_guard_text, FOIFI_COPY_IN, dest_iters, src_copy_level, dest_copy_level);
                    fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);
                    
                    fprintf(packfp, "int unpack_foifi_recv_%s_%d_%d(%s,%s,double *%s,int %s", 
                            acc_name_k, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
                    if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
                    fprintf(packfp,"){\n");
                    generate_pack_or_unpack(packfp, prog, foifi, flow_copyback_text, FOIFI_COPY_IN, iters, total_copy_level, acc_nrows);
                    fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);

                    for (i=0; i<dest_copy_level; i++) {
                        free(dest_iters[i]);
                    }
                    free(dest_iters);

                    pluto_constraints_free(fo);
                    pluto_constraints_free(foifi);
                    pluto_constraints_free(receiver_tiles);
                }

                sprintf(pack_stmt_text+strlen(pack_stmt_text), "} else {");
                sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "} else {");
            }
        }

        char sigma_text[1024];
        sprintf(sigma_text, "sigma_%s_%d(%s,%s", acc_name_k, loop_num, args, params);

        // !!!roshan should add the below keyword to parallelize per-receiver copy loop
        // this currently leads to slow down; should be investigated
            /*__omp_par_for_guided \*/
        sprintf(pack_stmt_text+strlen(pack_stmt_text), "clear_sender_receiver_lists(nprocs); %s, my_rank, nprocs); \
            for (__p=0; __p<nprocs; __p++) { if (receiver_list[__p] != 0) \
            { %s[__p] = pack_%s_%d(%s, %s[__p], %s[__p]", 
            sigma_text, send_counts_name, acc_name_k, loop_num, args, sendbufname, send_counts_name);
        if (options->variables_not_global) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
        sprintf(pack_stmt_text+strlen(pack_stmt_text), "); } }");

        if (options->fop_unicast_runtime && !broadcast) sprintf(pack_stmt_text+strlen(pack_stmt_text), "}");

        if(curr!= NULL && curr->deps == NULL){
        	curr->constraints = pluto_constraints_empty(src_copy_level + acc_nrows + prog->npar +1);
        }
        assert(curr->constraints->ncols == src_copy_level + acc_nrows + prog->npar + 1);
        
        fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, sendbufname, send_counts_name);
        if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
        fprintf(packfp,"){\n");
        generate_pack_or_unpack(packfp, prog, curr->constraints, flow_copy_text, COPY_OUT, iters, src_copy_level, acc_nrows);
        fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", send_counts_name);

        sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
            "if (is_receiver_%s_%d(%s,%s,my_rank,nprocs) != 0) \
            { %s[proc] = unpack_%s_%d(%s,%s,%s[proc]", 
            acc_name_k, loop_num, args, params, currdisplsname, acc_name_k, loop_num, args, recvbufname, currdisplsname);
        if (options->variables_not_global) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);
        sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "); }");

        if (options->fop_unicast_runtime && !broadcast) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");
        
        fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name_k, loop_num, decl_args, recvbufname, currdisplsname);
        if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
        fprintf(packfp,"){\n");
        generate_pack_or_unpack(packfp, prog, curr->constraints, flow_copyback_text, COPY_IN, iters, src_copy_level, acc_nrows);
        fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);

        k++;
        curr = curr->next;
    }
    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *pack_stmt = create_helper_stmt(anchor_stmt, src_copy_level, pack_stmt_text, COPY_OUT);
    Stmt *unpack_stmt = create_helper_stmt(anchor_stmt, src_copy_level,
        unpack_stmt_text, COPY_IN);

    // comm stmts organized as follows
    // pack stmt, comm stmt, unpack stmt
    *num_comm_stmts = 3;
    Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
    copy_comm_stmts[0] = pack_stmt;
    copy_comm_stmts[1] = comm_stmt;
    copy_comm_stmts[2] = unpack_stmt;

    for (i=0; i<acc_nrows; i++) {
        free(iters[i]);
    }
    free(iters);

    free(displsname);
    free(currdisplsname);
    free(sendbufname);
    free(send_counts_name);
    free(recv_counts_name);
    free(pack_stmt_text);
    free(unpack_stmt_text);
    free(flow_copy_text);
    free(flow_copyback_text);
    free(comm_text);
    free(displs_size);
    free(recv_buf_size);
    free(send_buf_size);
    free(access);
    fclose(packfp);

    pluto_constraints_list_free(atomic_flowouts);

    return copy_comm_stmts;
#ifdef PLUTO_OPENCL
    }
    else{
        pluto_opencl_codegen(); // just dummy for now
        return NULL;
    }
#endif

}

void free_stmt_array_buffers(Stmt** buffer, int size) {

    int i = 0;
    for(i = 0; i < size; ++i) {
        free(buffer[i]);
    }

    free(buffer);

    return;
}

/* 
 * Optimized communication code generation using FOIFI scheme
 * copy_level: number of outer loops to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop is normally the parallel loop)
 * wacc_stmts: all <statements, write access> of this data variable in this loop
 * num_accs: number of write accesses of this data variable in this loop
 * flow_in: flow_in set of this data variable in each loop - indexed by loop
 * This function is called per data variable
 */  
Stmt **gen_comm_code_opt_foifi(struct stmt_access_pair **wacc_stmts, int num_accs,
				int nloops,
        PlutoProg *prog, int *copy_level, int loop_num, int *pi_mappings, int *num_comm_stmts)
{
    int i, j, k, l, src_copy_level, acc_nrows;
    src_copy_level = copy_level[loop_num];
    assert(src_copy_level>=1);

    assert(num_accs >= 1);
    char *access = reconstruct_access(wacc_stmts[0]->acc);
    acc_nrows = wacc_stmts[0]->acc->mat->nrows;
    Stmt *anchor_stmt = wacc_stmts[0]->stmt;
    char *acc_name = anchor_stmt->writes[0]->name;

    /* To be inside a loop: can't foresee other use */
    assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
            prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

    PlutoConstraints *foifi_sets[nloops];
    for (l=0; l<nloops; l++) {
        foifi_sets[l] = NULL;
    }

    int broadcast = 0;
    /* Sender-side copying */
    PlutoConstraints *flow_out = NULL; // used only for get parametric extent/bounding box
    for (k=0; k<num_accs; k++)  {
        PlutoConstraints *flow_out_one;
        flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
        if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
        else{
            flow_out = pluto_constraints_unionize(flow_out, flow_out_one);
        }
        pluto_constraints_free(flow_out_one);

        if (broadcast) continue;

        for (i=0; i<prog->ndeps; i++)   {
            Dep *dep = prog->deps[i];
            /* Only RAW deps matter */
            if (dep->type != OSL_DEPENDENCE_RAW) continue;

            // if (dep->dirvec[copy_level+1] == DEP_ZERO) continue;

            /* If the dependence doesn't originate from this access */
            if (dep->src_acc != wacc_stmts[k]->acc) continue;

            assert(dep->dest_acc != NULL);

            Stmt *dest = prog->stmts[dep->dest];
            int dependent_loop = pi_mappings[dest->id];
            if (dependent_loop == -1) { 
                // destination statement will be executed by all processors
                broadcast = 1;
                break;
            }

            int dest_copy_level = copy_level[dependent_loop];

            PlutoConstraints *fo = compute_flow_out_of_dep(dep, src_copy_level, copy_level, prog, 0, NULL, pi_mappings);
            // add target tile iterators for flow-out set
            for (j=0; j<dest_copy_level; j++) {
                pluto_constraints_add_dim(fo,src_copy_level);
            }

            PlutoConstraints *fi = compute_flow_in_of_dep(dep, dest_copy_level, prog, 0);
            // add source tile iterators for flow-in set
            for (j=0; j<src_copy_level; j++) {
                pluto_constraints_add_dim(fi,0);
            }

            PlutoConstraints *foifi = pluto_constraints_intersection(fo, fi);

            if (foifi_sets[dependent_loop] == NULL) 
                foifi_sets[dependent_loop] = pluto_constraints_dup(foifi);
            else
                pluto_constraints_unionize(foifi_sets[dependent_loop], foifi);

            pluto_constraints_free(fo);
            pluto_constraints_free(fi);
            pluto_constraints_free(foifi);
        }
    }

    // === from here on we begin the code generation common for all dependent loops===================================
    //

    // !!!roshan may need to increase send buffer size to an arbitrarily large number
    // due to the possible amount of duplication in the FOIFI scheme
    char *send_buf_size = 
        get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows, 
                prog->npar, (const char **)prog->params);
    sprintf(send_buf_size+strlen(send_buf_size), "*__FOIFI_MAX_DUP_FACTOR");
    char *recv_buf_size = malloc(1024);
    char *displs_size = malloc(1024);
    strcpy(recv_buf_size, send_buf_size);
    strcpy(displs_size, send_buf_size);

    PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
    char *extent;
    /* Just the first one in flow_out is enough (rest all should give the
     * same since they are all under the same parallel loop and each
     * iteration of the parallel loop writes to distinct data) */
    get_parametric_extent_const(anchor_stmt_new_dom, src_copy_level-1, prog->npar,
            (const char **)prog->params, &extent);
    sprintf(send_buf_size+strlen(send_buf_size), 
            "*ceilf((%s)/(float)nprocs)", extent);
    /* The + nprocs is needed since the displacement has to be set larger
     * when some processors have more iterations - worst case when the
     * first processor has an extra iteration; one then needs
     * (nprocs-1)*num_values_per_iteration additional space) */
    sprintf(recv_buf_size+strlen(recv_buf_size), "*(%s + nprocs)", 
            extent);
    sprintf(displs_size+strlen(displs_size), "*ceilf((%s)/(float)nprocs)",
            extent);
    /* Assumes load-balanced distribution */
    free(extent);
    pluto_constraints_free(anchor_stmt_new_dom);
    IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

    char *sendbufname = concat("send_buf_", acc_name);
    char *send_counts_name = concat("send_counts_", acc_name);

    sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
            %s[__p] = 0;\n}\n", send_counts_name);
    sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
            %s[__p]\
            = (double *) polyrt_max_alloc(%s[__p], sizeof(double)*(%s), &send_buf_size_%s[__p]);\n}", 
            sendbufname, sendbufname, send_buf_size, acc_name);

    char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
            strlen(send_counts_name) + strlen("++] = ") + strlen(access) + 1);
    sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_counts_name, access);

    char *displsname = concat("displs_", acc_name);
    char *currdisplsname = concat("curr_displs_", acc_name);
    char *recvbufname = concat("recv_buf_", acc_name);
    char *recv_counts_name = concat("recv_counts_", acc_name);

    char *comm_text = malloc(2048);

    /* Message passing code (MPI calls) */
    if (options->commreport)
        sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
    else
        strcpy(comm_text, "");
    sprintf(comm_text+strlen(comm_text), "\
            MPI_Alltoall(%s, 1, MPI_INT,\
                %s, 1, MPI_INT, MPI_COMM_WORLD);", 
            send_counts_name, recv_counts_name);
    sprintf(comm_text+strlen(comm_text),  "\
            req_count=0;\
            for (__p=0; __p<nprocs; __p++) {\
            if (%s[__p] >= 1) {\
            assert(\"increase __FOIFI_MAX_DUP_FACTOR\" && (%s[__p]*8 <= send_buf_size_%s[__p]));",
            send_counts_name, send_counts_name, acc_name);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s[__p]);", send_counts_name);
    sprintf(comm_text+strlen(comm_text),  "\
            MPI_Isend(%s[__p], %s[__p], MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}", 
            sendbufname, send_counts_name);
    sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
            if(%s[__p] >= 1) {\
            MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
            MPI_Waitall(req_count, reqs, stats);\
            for (__p=0; __p<nprocs; __p++) {\
            %s[__p] = 0;}",
            recv_counts_name, recvbufname, displsname, recv_counts_name,
            send_counts_name);
    sprintf(comm_text+strlen(comm_text), "for (__p=0; __p<nprocs; __p++) {\
            %s[__p] = %s[__p]; }", currdisplsname, displsname);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *comm_stmt = create_helper_stmt(anchor_stmt, src_copy_level-1, comm_text, COMM_CALL);

    /* Receiver-side copy */
    /* With send/recv -based more exact communication scheme, recv side copy
     * does not include all the sender-side copies; processes that do not send
     * data are to be skipped */
    IF_DEBUG(printf("Recv buffer size for %s: %s\n", acc_name, recv_buf_size););

    sprintf(prog->decls+strlen(prog->decls), 
            "%s = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &recv_buf_size_%s);\n",
            recvbufname, recvbufname, recv_buf_size, acc_name);
    sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
            %s[__p] = __p*(%s);}\n\n", displsname, displs_size);

    char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
            + strlen("[") + strlen(currdisplsname) + strlen("++]") + 1);
    sprintf(flow_copyback_text, "%s = %s[%s++]",
            access, recvbufname, currdisplsname);

    char args[1024];
    strcpy(args,"");
    sprintf(args+strlen(args), "t%d", 1);
    /* make it src_copy_level+1 since we use an extra dimension to separate
     * statements */
    for (i=1; i<src_copy_level; i++) {
        sprintf(args+strlen(args), ",t%d", i+1);
    }

    char passed_args[1024];
    strcpy(passed_args,"");
    sprintf(passed_args+strlen(passed_args), "ts%d", 1);
    for (i=1; i<src_copy_level; i++) {
        sprintf(passed_args+strlen(passed_args), ",ts%d", i+1);
    }

    char decl_args[1024];
    strcpy(decl_args,"");
    sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
    for (i=1; i<src_copy_level; i++) {
        sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
    }

    char params[1024];
    strcpy(params, "");
    if (prog->npar>=1) {
        sprintf(params+strlen(params), "%s", prog->params[0]);
        for (i=1; i<prog->npar; i++) {
            sprintf(params+strlen(params), ",%s", prog->params[i]);
        }
    }

    char **iters;
    iters = malloc(acc_nrows * sizeof(char *));
    for (i=0; i < acc_nrows; i++) {
        iters[i] = malloc(5);
        sprintf(iters[i], "d%d", i+1);
    }

    FILE *packfp = fopen("packunpack.c", "a");
    assert(packfp != NULL);

    char *pack_stmt_text = malloc(8192);
    strcpy(pack_stmt_text,"");
    char *unpack_stmt_text = malloc(8192);
    sprintf(unpack_stmt_text, "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "if ((my_rank != proc) && (%s[proc] > 0)) { ", recv_counts_name);
    if (broadcast) {
        if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
            sprintf(pack_stmt_text+strlen(pack_stmt_text), 
                    "for (__p=0; __p<nprocs; __p++) \
                    { if (__p != my_rank) \
                    { %s[__p] = pack_%s_%d(%s,%s[__p],%s[__p]", 
                    send_counts_name, acc_name, loop_num, args, sendbufname, send_counts_name);
            if (options->variables_not_global) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
            sprintf(pack_stmt_text+strlen(pack_stmt_text), "); } }");

            fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_counts_name);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            generate_pack_or_unpack(packfp, prog, flow_out, flow_copy_text, COPY_OUT, iters, src_copy_level, acc_nrows);
            fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", send_counts_name);

            sprintf(unpack_stmt_text+strlen(unpack_stmt_text),
                "%s[proc] = unpack_%s_%d(%s,%s,%s[proc]", 
                currdisplsname, acc_name, loop_num, args, recvbufname, currdisplsname);
            if (options->variables_not_global) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);
            sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

            fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, recvbufname, currdisplsname);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            generate_pack_or_unpack(packfp, prog, flow_out, flow_copyback_text, COPY_IN, iters, src_copy_level, acc_nrows);
            fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);
        }
    }
    else {
        for (l=0; l<nloops; l++) {
            PlutoConstraints *foifi = foifi_sets[l];
            if (foifi == NULL) continue;
            if (pluto_constraints_is_empty(foifi)) {
                pluto_constraints_free(foifi);
                continue;
            }
            
            int dest_copy_level = copy_level[l];
            assert(dest_copy_level>=1);
            int total_copy_level = src_copy_level + dest_copy_level;
            assert(foifi->ncols == total_copy_level + acc_nrows + prog->npar + 1);

            PlutoConstraints *receiver_tiles = get_receiver_tiles(wacc_stmts, num_accs, src_copy_level, dest_copy_level, prog, l, pi_mappings);
            assert(receiver_tiles->ncols == total_copy_level + prog->npar + 1);
            assert(!pluto_constraints_is_empty(receiver_tiles));

            IF_DEBUG(printf("Data flow out intersection flow in set for %s and dependent loop %d\n", acc_name, l););
            IF_DEBUG(pluto_constraints_print(stdout, foifi));

        // === from here on we begin the code generation for each dependent loop==========================================
        //
#ifdef PLUTO_OPENCL
            if(!options->opencl) {
#endif

            char dest_args[1024];
            strcpy(dest_args,"");
            sprintf(dest_args+strlen(dest_args), "t%d", src_copy_level+1);
            for (i=src_copy_level+1; i<total_copy_level; i++) {
                sprintf(dest_args+strlen(dest_args), ",t%d", i+1);
            }

            char decl_dest_args[1024];
            strcpy(decl_dest_args,"");
            sprintf(decl_dest_args+strlen(decl_dest_args), "int ts%d", src_copy_level+1);
            for (i=src_copy_level+1; i<total_copy_level; i++) {
                sprintf(decl_dest_args+strlen(decl_dest_args), ",int ts%d", i+1);
            }

            char **dest_iters;
            dest_iters = malloc(dest_copy_level * sizeof(char *));
            for (i=0; i<dest_copy_level; i++) {
                dest_iters[i] = malloc(5);
                sprintf(dest_iters[i], "t%d", i+src_copy_level+1);
            }

            char pi_text[1024];
            sprintf(pi_text, "recv_proc = pi_%d(%s,%s, nprocs)", l, dest_args, params);

            sprintf(pack_stmt_text+strlen(pack_stmt_text), "pack_%s_%d_%d(%s,%s,%s, my_rank, nprocs", 
                    acc_name, loop_num, l, args, sendbufname, send_counts_name);
            if (options->variables_not_global) sprintf(pack_stmt_text+strlen(pack_stmt_text), ",%s", acc_name);
            sprintf(pack_stmt_text+strlen(pack_stmt_text), ");");

            char flow_cg_text[1024];
            sprintf(flow_cg_text, "%s; if (recv_proc != my_rank) \
                    { %s[recv_proc] = pack_recv_%s_%d_%d(%s,%s,%s[recv_proc],%s[recv_proc]", 
                    pi_text, send_counts_name, acc_name, loop_num, l, passed_args, dest_args, sendbufname, send_counts_name);
            if (options->variables_not_global) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);
            sprintf(flow_cg_text+strlen(flow_cg_text), "); }");
            
            fprintf(packfp, "int pack_%s_%d_%d(%s,double **%s,int *%s, int my_rank, int nprocs", 
                    acc_name, loop_num, l, decl_args, sendbufname, send_counts_name);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            fprintf(packfp, "\nint recv_proc;\n");
            generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_cg_text, COPY_OUT, dest_iters, src_copy_level, dest_copy_level);
            fprintf(packfp, "\nreturn 0;\n#undef S1;\n}\n\n");
            
            fprintf(packfp, "int pack_recv_%s_%d_%d(%s,%s,double *%s,int %s", 
                    acc_name, loop_num, l, decl_args, decl_dest_args, sendbufname, send_counts_name);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            generate_pack_or_unpack(packfp, prog, foifi, flow_copy_text, COPY_OUT, iters, total_copy_level, acc_nrows);
            fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", send_counts_name);

            sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "\
                    %s[proc] = unpack_%s_%d_%d(%s,%s,%s[proc], my_rank, nprocs", 
                    currdisplsname, acc_name, loop_num, l, args, recvbufname, currdisplsname);
            if (options->variables_not_global) sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ",%s", acc_name);
            sprintf(unpack_stmt_text+strlen(unpack_stmt_text), ");");

            char flow_copyback_guard_text[1024];
            sprintf(flow_copyback_guard_text, "%s; if (recv_proc == my_rank) \
                    { %s = unpack_recv_%s_%d_%d(%s,%s,%s,%s", 
                    pi_text, currdisplsname, acc_name, loop_num, l, passed_args, dest_args, recvbufname, currdisplsname);
            if (options->variables_not_global) sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), ",%s", acc_name);
            sprintf(flow_copyback_guard_text+strlen(flow_copyback_guard_text), "); }");
            
            fprintf(packfp, "int unpack_%s_%d_%d(%s,double *%s,int %s, int my_rank, int nprocs", acc_name, loop_num, l, decl_args, recvbufname, currdisplsname);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            fprintf(packfp, "\nint recv_proc;\n");
            generate_pack_or_unpack(packfp, prog, receiver_tiles, flow_copyback_guard_text, FOIFI_COPY_IN, dest_iters, src_copy_level, dest_copy_level);
            fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);
            
            fprintf(packfp, "int unpack_recv_%s_%d_%d(%s,%s,double *%s,int %s", acc_name, loop_num, l, decl_args, decl_dest_args, recvbufname, currdisplsname);
            if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
            fprintf(packfp,"){\n");
            generate_pack_or_unpack(packfp, prog, foifi, flow_copyback_text, FOIFI_COPY_IN, iters, total_copy_level, acc_nrows);
            fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", currdisplsname);

            for (i=0; i<dest_copy_level; i++) {
                free(dest_iters[i]);
            }
            free(dest_iters);

            pluto_constraints_free(foifi);
            pluto_constraints_free(receiver_tiles);

#ifdef PLUTO_OPENCL
            } 
            else {
                pluto_opencl_codegen(); // dummy for now
            }
#endif
        }
    }
    sprintf(unpack_stmt_text+strlen(unpack_stmt_text), "}");

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *pack_stmt = create_helper_stmt(anchor_stmt, src_copy_level, pack_stmt_text, FOIFI_COPY_OUT);
    Stmt *unpack_stmt = create_helper_stmt(anchor_stmt, src_copy_level, 
            unpack_stmt_text, COPY_IN);

    // comm stmts organized as follows
    // pack stmt, comm stmt, unpack stmt
    *num_comm_stmts = 3;
    Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
    copy_comm_stmts[0] = pack_stmt;
    copy_comm_stmts[1] = comm_stmt;
    copy_comm_stmts[2] = unpack_stmt;

    for (i=0; i<acc_nrows; i++) {
        free(iters[i]);
    }
    free(iters);

    free(pack_stmt_text);
    free(unpack_stmt_text);
    free(flow_copy_text);
    free(flow_copyback_text);
    free(comm_text);
    free(displs_size);

    free(displsname);
    free(sendbufname);
    free(send_counts_name);
    free(recv_counts_name);
    free(recv_buf_size);
    free(send_buf_size);
    fclose(packfp);

    pluto_constraints_free(flow_out);
    free(access);
    return copy_comm_stmts;
}

/* 
 * Optimized communication code generation
 * copy_level: copy level of all loops
 * loop_num: source loop
 * copy_level[loop_num]: number of outer loops of source loop to be treated as parameters
 * (copy_level-1)^th (0-indexed) loop is normally the parallel loop)
 * wacc_stmts:  all <statements, wacc> writing to this data variable
 * This function is called per data variable
 */  
Stmt **gen_comm_code_opt(struct stmt_access_pair **wacc_stmts, int num_accs,
        PlutoProg *prog, int *copy_level, int loop_num, int *pi_mappings, int *num_comm_stmts)
{
    int i, k, src_copy_level, acc_nrows;
    src_copy_level = copy_level[loop_num];
    assert(src_copy_level>=1);

    assert(num_accs >= 1);
    char *access = reconstruct_access(wacc_stmts[0]->acc);
    acc_nrows = wacc_stmts[0]->acc->mat->nrows;
    Stmt *anchor_stmt = wacc_stmts[0]->stmt;
    char *acc_name = anchor_stmt->writes[0]->name;

    // printf("%d accesses: %s\n", num_accs, acc_name);

    /* To be inside a loop: can't foresee other use */
    assert(prog->hProps[src_copy_level-1].type == H_LOOP ||
            prog->hProps[src_copy_level-1].type == H_TILE_SPACE_LOOP);

    /* Sender-side copying */
    PlutoConstraints *flow_out = NULL;
    for (k=0; k<num_accs; k++)  {
        PlutoConstraints *flow_out_one;
        flow_out_one = compute_flow_out(wacc_stmts[k], src_copy_level, copy_level, prog, pi_mappings);
        if (flow_out == NULL) flow_out = pluto_constraints_dup(flow_out_one);
        else{
            flow_out = pluto_constraints_unionize(flow_out, flow_out_one);
        }
        pluto_constraints_free(flow_out_one);
    }
	IF_DEBUG(print_polylib_visual_sets("data", flow_out));

    generate_sigma(wacc_stmts, num_accs, copy_level, prog, loop_num, pi_mappings);
    // generate_tau(wacc_stmts, num_accs, copy_level, prog);

    IF_DEBUG(printf("Data flow out set for %s\n", acc_name););
    IF_DEBUG(pluto_constraints_print(stdout, flow_out));
    assert(flow_out->ncols == src_copy_level + acc_nrows + prog->npar + 1);

// ==== from here we start the code generation =================================
//
#ifdef PLUTO_OPENCL
    if(!options->opencl) {
#endif

    char *send_buf_size = 
        get_parametric_bounding_box(flow_out, src_copy_level, acc_nrows, 
                prog->npar, (const char **)prog->params);
    char *recv_buf_size = malloc(1024);
    char *displs_size = malloc(1024);
    strcpy(recv_buf_size, send_buf_size);
    strcpy(displs_size, send_buf_size);

    PlutoConstraints *anchor_stmt_new_dom = pluto_get_new_domain(anchor_stmt);
    char *extent;
    /* Just the first one in flow_out is enough (rest all should give the
     * same since they are all under the same parallel loop and each
     * iteration of the parallel loop writes to distinct data) */
    get_parametric_extent_const(anchor_stmt_new_dom, src_copy_level-1, prog->npar,
            (const char **)prog->params, &extent);
    // printf("Extent is %s\n", extent);
    sprintf(send_buf_size+strlen(send_buf_size), 
            "*ceilf((%s)/(float)nprocs)", extent);
    /* The + nprocs is needed since the displacement has to be set larger
     * when some processors have more iterations - worst case when the
     * first processor has an extra iteration; one then needs
     * (nprocs-1)*num_values_per_iteration additional space) */
    sprintf(recv_buf_size+strlen(recv_buf_size), "*(%s + nprocs)", 
            extent);
    sprintf(displs_size+strlen(displs_size), "*ceilf((%s)/(float)nprocs)",
            extent);
    /* Assumes load-balanced distribution */
    free(extent);
    pluto_constraints_free(anchor_stmt_new_dom);
    IF_DEBUG(printf("Send buffer size for %s: %s\n", acc_name, send_buf_size););

    FILE *packfp = fopen("packunpack.c", "a");
    assert(packfp != NULL);

    char args[1024];
    strcpy(args,"");
    sprintf(args+strlen(args), "t%d", 1);
    /* make it src_copy_level+1 since we use an extra dimension to separate
     * statements */
    for (i=1; i<src_copy_level; i++) {
        sprintf(args+strlen(args), ",t%d", i+1);
    }

    char decl_args[1024];
    strcpy(decl_args,"");
    sprintf(decl_args+strlen(decl_args), "int ts%d", 1);
    for (i=1; i<src_copy_level; i++) {
        sprintf(decl_args+strlen(decl_args), ",int ts%d", i+1);
    }

    char params[1024];
    strcpy(params, "");
    if (prog->npar>=1) {
        sprintf(params+strlen(params), "%s", prog->params[0]);
        for (i=1; i<prog->npar; i++) {
            sprintf(params+strlen(params), ",%s", prog->params[i]);
        }
    }

    char **iters;
    iters = malloc(acc_nrows * sizeof(char *));
    for (i=0; i < acc_nrows; i++) {
        iters[i] = malloc(5);
        sprintf(iters[i], "d%d", i+1);
    }

    /* Sender-side copy */
    char *sendbufname = concat("send_buf_", acc_name);
    char *send_counts_name = concat("send_counts_", acc_name);
    char *send_count_name = concat("send_count_", acc_name);

    sprintf(prog->decls+strlen(prog->decls), "%s\
            = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &send_buf_size_%s);\n", 
            sendbufname, sendbufname, send_buf_size, acc_name);

    char *flow_copy_text = malloc(strlen(sendbufname) + strlen("[")+
            strlen(send_count_name) + strlen("++] = ") + strlen(access) + 1);
    sprintf(flow_copy_text, "%s[%s++] = %s", sendbufname, send_count_name, access);

    char *sigma_text = malloc(1024);
    sprintf(sigma_text, "sigma_%s_%d(%s,%s", acc_name, loop_num, args, params);

    char *sigma_sender_text = malloc(1024);
    strcpy(sigma_sender_text, sigma_text);
    sprintf(sigma_sender_text+strlen(sigma_sender_text), ", my_rank, nprocs)");

    char *flow_cg_text = malloc(1024);
    strcpy(flow_cg_text, "");

    if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
        sprintf(flow_cg_text+strlen(flow_cg_text), "clear_sender_receiver_lists(nprocs); %s; \
                if (need_to_send(nprocs)) \
                { %s = pack_%s_%d(%s,%s,%s", 
                sigma_sender_text, send_count_name, acc_name, loop_num, args, sendbufname, send_count_name);
        if (options->variables_not_global) sprintf(flow_cg_text+strlen(flow_cg_text), ",%s", acc_name);
        sprintf(flow_cg_text+strlen(flow_cg_text), "); }");

        fprintf(packfp, "int pack_%s_%d(%s,double *%s,int %s", acc_name, loop_num, decl_args, sendbufname, send_count_name);
        if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
        fprintf(packfp,"){\n");
        generate_pack_or_unpack(packfp, prog, flow_out, flow_copy_text, COPY_OUT, iters, src_copy_level, acc_nrows);
        fprintf(packfp, "\nreturn %s;\n#undef S1;\n}\n\n", send_count_name);
    }

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *flow_copy_guard = create_helper_stmt(anchor_stmt, src_copy_level, flow_cg_text, COPY_OUT);

    /* sigma, clear statements */
    char *clear_text = "clear_sender_receiver_lists(nprocs)";

    Stmt *send_recv_list_clear = create_helper_stmt(anchor_stmt, src_copy_level-1, clear_text, COPY_OUT);
    Stmt *sigma_stmt = create_helper_stmt(anchor_stmt, src_copy_level, sigma_sender_text, SIGMA);

    char *displsname = concat("displs_", acc_name);
    char *recvbufname = concat("recv_buf_", acc_name);
    char *recv_counts_name = concat("recv_counts_", acc_name);

    char *comm_text = malloc(2048);

    /* Message passing code (MPI calls) */
    if (options->commreport)
        sprintf(comm_text,  "IF_TIME(t_comm_start = rtclock());");
    else
        strcpy(comm_text, "");
    sprintf(comm_text+strlen(comm_text), "\
            for (__p=0; __p<nprocs; __p++) {\
            %s[__p] = receiver_list[__p]? %s: 0;\
            }\
            MPI_Alltoall(%s, 1, MPI_INT,\
                %s, 1, MPI_INT, MPI_COMM_WORLD);", 
            send_counts_name, send_count_name, send_counts_name, recv_counts_name);
    sprintf(comm_text+strlen(comm_text),  "\
            req_count=0;\
            for (__p=0; __p<nprocs; __p++) {\
            if(%s[__p] >= 1) {",
            send_counts_name);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(__total_count += %s);", send_count_name);
    sprintf(comm_text+strlen(comm_text),  "\
            MPI_Isend(%s, %s, MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}", 
            sendbufname, send_count_name);
    sprintf(comm_text+strlen(comm_text),  "for (__p=0; __p<nprocs; __p++) {\
            if(%s[__p] >= 1) {\
            MPI_Irecv(%s+%s[__p], %s[__p], MPI_DOUBLE,\
                __p, 123, MPI_COMM_WORLD, &reqs[req_count++]);}}\
            MPI_Waitall(req_count, reqs, stats);\
            %s = 0; prev_proc = -1; ", 
            recv_counts_name, recvbufname, displsname, recv_counts_name,
            send_count_name);
    if (options->commreport)
        sprintf(comm_text+strlen(comm_text),  "IF_TIME(t_comm += rtclock() - t_comm_start);");

    /*sprintf(text, "MPI_Allgather(send_buf, send_count, MPI_DOUBLE,\
      recv_buf, max_recv_count, MPI_DOUBLE, MPI_COMM_WORLD)");*/

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *comm_stmt = create_helper_stmt(anchor_stmt, src_copy_level-1, comm_text, COMM_CALL);
    //IF_DEBUG(pluto_stmt_print(stdout, comm_stmt););

    /* Receiver-side copy */
    /* With send/recv -based more exact communication scheme, recv side copy
     * does not include all the sender-side copies; processes that do not send
     * data are to be skipped */
    IF_DEBUG(printf("Recv buffer size for %s: %s\n", acc_name, recv_buf_size););

    sprintf(prog->decls+strlen(prog->decls), 
            "%s = (double *) polyrt_max_alloc(%s, sizeof(double)*(%s), &recv_buf_size_%s);\n",
            recvbufname, recvbufname, recv_buf_size, acc_name);
    sprintf(prog->decls+strlen(prog->decls), "\tfor (__p=0; __p<nprocs; __p++) {\
            %s[__p] = __p*(%s);}\n\n", displsname, displs_size);

    char *flow_copyback_text = malloc(strlen(access) + strlen(" = ") + strlen(recvbufname)
            + strlen("[") + strlen(displsname) + strlen("[proc]+ count++];") + 1);
    sprintf(flow_copyback_text, "%s = %s[%s + count++];",
            access, recvbufname, displsname);

    char proc_stmt_text[512];
    strcpy(proc_stmt_text, "");

    if ((flow_out != NULL) && !pluto_constraints_is_empty(flow_out)) {
        sprintf(proc_stmt_text+strlen(proc_stmt_text), "proc = pi_%d(%s,%s, nprocs); ", loop_num, args, params);
        sprintf(proc_stmt_text+strlen(proc_stmt_text), "clear_sender_receiver_lists(nprocs); \
                %s, proc, nprocs); ", sigma_text);
        sprintf(proc_stmt_text+strlen(proc_stmt_text), 
                "; if (proc != prev_proc) {prev_proc = proc; count=0;}\
                if ((%s[proc] > 0) && need_to_send(nprocs)) \
                { count = unpack_%s_%d(%s,%s,%s[proc],count", 
                recv_counts_name, acc_name, loop_num, args, recvbufname, displsname);
        if (options->variables_not_global) sprintf(proc_stmt_text+strlen(proc_stmt_text), ",%s", acc_name);
        sprintf(proc_stmt_text+strlen(proc_stmt_text), "); }");

        fprintf(packfp, "int unpack_%s_%d(%s,double *%s,int %s,int count", acc_name, loop_num, decl_args, recvbufname, displsname);
        if (options->variables_not_global) fprintf(packfp, ", __DECLARATION_OF_%s", acc_name);
        fprintf(packfp,"){\n");
        generate_pack_or_unpack(packfp, prog, flow_out, flow_copyback_text, COPY_IN, iters, src_copy_level, acc_nrows);
        fprintf(packfp, "\nreturn count;\n#undef S1;\n}\n\n");
    }

    /* Create statement stub to get hold of the outer loops for the copy
     * stmt; rest of the loops to be added are the actual copy loops */
    Stmt *copyback_proc_stmt = create_helper_stmt(anchor_stmt, src_copy_level, 
            proc_stmt_text, COPY_IN);

    // comm stmts organized as follows
    // pack stmt, clear stmt, sigma stmt, comm stmt, unpack stmt
    *num_comm_stmts = 5;
    Stmt **copy_comm_stmts = (Stmt **) malloc((*num_comm_stmts)*sizeof(Stmt *));
    copy_comm_stmts[0] = flow_copy_guard;
    copy_comm_stmts[1] = send_recv_list_clear;
    copy_comm_stmts[2] = sigma_stmt;
    copy_comm_stmts[3] = comm_stmt;
    copy_comm_stmts[4] = copyback_proc_stmt;

    for (i=0; i<acc_nrows; i++) {
        free(iters[i]);
    }
    free(iters);

    pluto_constraints_free(flow_out);
    free(sigma_text);
    free(flow_copy_text);
    free(flow_cg_text);
    free(flow_copyback_text);
    free(comm_text);
    free(displs_size);

    free(displsname);
    free(sendbufname);
    free(send_count_name);
    free(send_counts_name);
    free(recv_counts_name);
    free(recv_buf_size);
    free(send_buf_size);
    free(access);
    fclose(packfp);

    return copy_comm_stmts;

#ifdef PLUTO_OPENCL
    }
    else {
        pluto_opencl_codegen(); // dummy for now
        return NULL;
    }
#endif

}



int pluto_distmem_codegen(PlutoProg *prog, FILE *cloogfp, FILE *outfp)
{
    if (options->commopt_foifi) {
        fprintf(outfp, "#include <assert.h>\n\n");
        fprintf(outfp, "##ifndef __FOIFI_MAX_DUP_FACTOR\n");
        fprintf(outfp, "##define __FOIFI_MAX_DUP_FACTOR 1\n");
        fprintf(outfp, "##endif\n\n");
    }
    fprintf(outfp, "#include <mpi.h>\n\n");
    fprintf(outfp, "\tvoid *polyrt_max_alloc(void *buf, size_t size, size_t *curr_size);\n\n");
    fprintf(outfp, "#define MPI \n\n");
    fprintf(outfp, "#include <limits.h>\n\n");


    generate_declarations(prog, outfp);

    if (options->commreport) {
        fprintf(outfp, "\tdouble t_comm_start, t_comm = 0.0, t_globalcomm = 0.0;\n");
        fprintf(outfp, "\tdouble t_comp_start, t_comp = 0.0, t_globalcomp = 0.0;\n");
        fprintf(outfp, "\tdouble t_pack_start, t_pack = 0.0, t_globalpack = 0.0;\n");
        fprintf(outfp, "\tdouble t_writeout_start, t_writeout = 0.0;\n");
        fprintf(outfp, "\tdouble t_unpack_start, t_unpack = 0.0, t_globalunpack = 0.0;\n");
        fprintf(outfp, "\tdouble t_local = 0.0, t_global = 0.0;\n");
        fprintf(outfp, "\tdouble __total_count = 0, __total_count_all = 0;\n");
    }

    fprintf(outfp, "\n##ifndef GLOBAL_MY_RANK\n\tint my_rank;\n##endif\n");
    fprintf(outfp, "\tint nprocs, my_start, my_end, _i,\
            __p, proc, recv_proc, prev_proc, lw_prev_proc, distinct_recv;\n");
    fprintf(outfp, "\tint count;\n");
    fprintf(outfp, "int req_count;\n");
    fprintf(outfp, "\tint _lb_dist, _ub_dist;\n");
    fprintf(outfp, "\tMPI_Init(NULL, NULL);\n");
    fprintf(outfp, "\tMPI_Comm_rank(MPI_COMM_WORLD, &my_rank);\n");
    fprintf(outfp, "\tMPI_Comm_size(MPI_COMM_WORLD, &nprocs);\n");
    fprintf(outfp, "%s", prog->decls);
    fprintf(outfp, "\textern int *sender_list, *receiver_list; \n");
    fprintf(outfp, "MPI_Request %s[2*nprocs];\n", "reqs");
    fprintf(outfp, "MPI_Status %s[2*nprocs];\n\n", "stats");

    fprintf(outfp, "\tpolyrt_init(nprocs);\n\n");
    if (options->commreport) {
        fprintf(outfp, "\tt_local = rtclock();\n");
    }

    pluto_gen_cloog_code(prog, -1, -1, cloogfp, outfp);

    if (options->commreport) {
        fprintf(outfp, "##ifdef TIME\n");
        fprintf(outfp, "\tt_local = rtclock() - t_local - t_writeout;\n");
        fprintf(outfp, "\tchar buffer[4096];\n");
        fprintf(outfp, "\tstrcpy(buffer, \"\");\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Node %%d\\n\", my_rank);\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Write-out time: %%lf\\n\", t_writeout);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Values communicated: %%0.0lf\\n\", __total_count);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Computation time: %%lf\\n\", t_comp);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Communication time: %%lf\\n\", t_comm);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Packing time: %%lf\\n\", t_pack);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Unpacking time: %%lf\\n\", t_unpack);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Sum of split times: %%lf\\n\", t_comp + t_comm + t_pack + t_unpack);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"Total time minus write-out time: %%lf\\n\", t_local);\n\n");
        fprintf(outfp, "\tsprintf(buffer+strlen(buffer), \"-------\");\n");
        fprintf(outfp, "\tfprintf(stdout, \"%%s\\n\", buffer);\n");
        fprintf(outfp, "\tMPI_Reduce(&__total_count, &__total_count_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tMPI_Reduce(&t_comp, &t_globalcomp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tMPI_Reduce(&t_comm, &t_globalcomm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tMPI_Reduce(&t_pack, &t_globalpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tMPI_Reduce(&t_unpack, &t_globalunpack, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tMPI_Reduce(&t_local, &t_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);\n");
        fprintf(outfp, "\tif (my_rank==0) {\n");
        fprintf(outfp, "\t\tstrcpy(buffer, \"SUMMARY\\n\");\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Write-out time spent in master node: %%0.6lf s\\n\", t_writeout);\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Total communication volume across all nodes: %%0.6lf GB\\n\", __total_count_all*8/(1024*1024*1024));\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum computation time spent across all nodes: %%0.6lf s\\n\", t_globalcomp);\n\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum communication time spent across all nodes: %%0.6lf s\\n\", t_globalcomm);\n\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum packing time spent across all nodes: %%0.6lf s\\n\", t_globalpack);\n\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum unpacking time spent across all nodes: %%0.6lf s\\n\", t_globalunpack);\n\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"Maximum total time spent across all nodes: %%0.6lf s\\n\", t_global);\n");
        fprintf(outfp, "\t\tsprintf(buffer+strlen(buffer), \"-------\");\n");
        fprintf(outfp, "\t\tfprintf(stdout, \"%%s\\n\", buffer);\n");
        fprintf(outfp, "\t}\n");
        fprintf(outfp, "##endif\n");
    }
    fprintf(outfp, "\tMPI_Finalize();\n");

    return 0; 
}

void pluto_add_distmem_decls(PlutoProg *prog)
{
    int i, num=0;

    PlutoAccess **waccs;
    waccs = pluto_get_all_waccs(prog, &num);

    for (i=0; i<num; i++) {
        char *name = waccs[i]->name;

        /* lw count name */
        sprintf(prog->decls+strlen(prog->decls), 
                "int lw_count_%s = 0;\n", name );
        sprintf(prog->decls+strlen(prog->decls), 
                "int send_counts_%s[nprocs];\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "int recv_counts_%s[nprocs];\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "int lw_recv_counts_%s[nprocs];\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "int displs_%s[nprocs];\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "int displs_lw_%s[nprocs];\n", name);
        if (options->commopt_foifi || options->commopt_fop) {
            sprintf(prog->decls+strlen(prog->decls), 
                    "double *send_buf_%s[nprocs];\n", name);
            sprintf(prog->decls+strlen(prog->decls), 
                    "int curr_displs_%s[nprocs];\n", name);
            sprintf(prog->decls+strlen(prog->decls), 
                    "size_t send_buf_size_%s[nprocs];\n", name);
            sprintf(prog->decls+strlen(prog->decls), "for (__p=0; __p<nprocs; __p++) {\n\
                    send_buf_%s[__p] = NULL;\n\
                    send_buf_size_%s[__p] = 0;\n}", name, name);
        }
        else {
            /* send buf count name */
            sprintf(prog->decls+strlen(prog->decls), 
                    "int send_count_%s = 0;\n", name);
            sprintf(prog->decls+strlen(prog->decls), 
                    "double *send_buf_%s = NULL;\n", name);
            sprintf(prog->decls+strlen(prog->decls), 
                    "size_t send_buf_size_%s = 0;\n", name);
        }
        /* To store previous size */
        sprintf(prog->decls+strlen(prog->decls), 
                "double *recv_buf_%s = NULL;\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "size_t recv_buf_size_%s = 0;\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "double *lw_buf_%s = NULL;\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "size_t lw_buf_size_%s = 0;\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "double *lw_recv_buf_%s = NULL;\n", name);
        sprintf(prog->decls+strlen(prog->decls), 
                "size_t lw_recv_buf_size_%s = 0;\n", name);
    }

}


/* Distributed memory parallelization */
int pluto_distmem_parallelize(PlutoProg *prog)
{
    int i, j, l, k, dist_parallel_loop=0;
    int nstmts;
    Stmt ****copy_comm_stmts, ****write_out_stmts;
    Stmt ***sep_stmts;
    int *pi_mappings, *num_data, *copy_level;

    nstmts = prog->nstmts;

    // print_hyperplane_properties(prog);

    /* Create a new pi.c file */

    FILE *pifp = fopen("pi.c", "w");
    assert(pifp != NULL);

    fprintf(pifp, "\
#include <math.h>\n\
#include <assert.h>\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

    /* Create a new sigma.c file */

    FILE *sigmafp = NULL;
    if (options->commopt_fop) {
        sigmafp = fopen("sigma_fop.c", "w");
        if (options->fop_unicast_runtime) {
            fprintf(sigmafp, "#ifndef __FOP_UNICAST_RECV_LIMIT\n");
            fprintf(sigmafp, "#define __FOP_UNICAST_RECV_LIMIT 1\n");
            fprintf(sigmafp, "#endif\n\n");
        }
    }else if (!options->commopt_foifi) {
        sigmafp = fopen("sigma.c", "w");
    }
    if (sigmafp != NULL) {
        fprintf(sigmafp, "\
#include <math.h>\n\
#include <assert.h>\n\
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))\n\
#define floord(n,d) floor(((double)(n))/((double)(d)))\n\
#define max(x,y)    ((x) > (y)? (x) : (y))\n\
#define min(x,y)    ((x) < (y)? (x) : (y))\n\n");

        fclose(sigmafp);
    }

    pluto_add_distmem_decls(prog);

    int nloops=0;
    Ploop **loops = pluto_get_dom_parallel_loops(prog, &nloops);
    /* Loops should be in the increasing order of their depths */
    qsort(loops, nloops, sizeof(Ploop *), pluto_loop_compar);

    IF_DEBUG(printf("distmem: parallelizing loops\n"););
    IF_DEBUG(pluto_loops_print(loops,nloops););

    pi_mappings = malloc(nstmts*sizeof(int));
    for (i=0; i<prog->nstmts; i++) {
        pi_mappings[i] = -1;
    }

    for (l=0; l<nloops; l++) {
        for (i=0; i<loops[l]->nstmts; i++) {
            pi_mappings[loops[l]->stmts[i]->id] = l;
        }
    }
    // if pi_mappings of a statement is -1, then that statement is not distributed

    /* Stripmine parallel loop to achieve block-cyclic partitioning */
    if (options->blockcyclic) {
        for (l=0; l<nloops; l++) {
            int tsizes[1];
            Ploop *bloop = pluto_loop_dup(loops[l]);
            bloop->depth += l;
            Band band = {bloop, 1};
            tsizes[0]=options->cyclesize;
            pluto_tile_band(prog, &band, tsizes);
        }
    }

    copy_comm_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
    write_out_stmts = (Stmt ****) malloc(nloops*sizeof(Stmt ***));
    sep_stmts = (Stmt ***) malloc(nloops*sizeof(Stmt **));

    num_data = (int *) malloc(nloops*sizeof(int));
    copy_level = (int *) malloc(nloops*sizeof(int));

    int num_comm_stmts = 0;

    char packfilename[15];
    strcpy(packfilename,"packunpack.c");
    FILE *appendfp = fopen(".appendfilename", "w");
    assert(appendfp != NULL);
    fprintf(appendfp, "%s\n", packfilename);
    fclose(appendfp);
    FILE *packfp = fopen(packfilename, "w");
    fclose(packfp);

    for (l=0; l<nloops; l++) {
        if (options->blockcyclic) {
            /* FIXME: copy_level[l] in the presence of blockcyclic
             * will depend on the order in which these parallel loops were tiled 
             * for block cylic partitioning (will work if they were in the
             * increasing order of depths: i.e., loops in ploops array were in
             * increasing order of their depths
             */
            copy_level[l] = loops[l]->depth + l + 1;
        }else{
            copy_level[l] = loops[l]->depth + 1;
        }
    }

    for (l=0; l<nloops; l++) {
        Ploop *loop = loops[l];

        int *num_stmts_per_wacc; // indexed by data variable
        struct stmt_access_pair ***wacc_stmts; // indexed by data variable
        wacc_stmts = get_write_access_with_stmts(loop->stmts, 
            loop->nstmts, &num_data[l], &num_stmts_per_wacc);

        copy_comm_stmts[l] = (Stmt ***) malloc(num_data[l]*sizeof(Stmt **));
        write_out_stmts[l] = (Stmt ***) malloc(num_data[l]*sizeof(Stmt **));
        for (i=0; i<num_data[l]; i++)    {
            if (options->commopt_fop) {
                copy_comm_stmts[l][i] = 
                    gen_comm_code_opt_fop(wacc_stmts[i],num_stmts_per_wacc[i],nloops,prog,
                        copy_level,l,pi_mappings,&num_comm_stmts);
            }else if (options->commopt_foifi) {
                copy_comm_stmts[l][i] = 
                    gen_comm_code_opt_foifi(wacc_stmts[i],num_stmts_per_wacc[i],nloops,prog,
                        copy_level,l,pi_mappings,&num_comm_stmts);
            }else {
                copy_comm_stmts[l][i] = 
                    gen_comm_code_opt(wacc_stmts[i],num_stmts_per_wacc[i],prog,
                        copy_level,l,pi_mappings,&num_comm_stmts);
            }
            write_out_stmts[l][i] = 
                gen_write_out_code(wacc_stmts[i],num_stmts_per_wacc[i],prog,copy_level,l);
        }
        for (i=0; i<num_data[l]; i++) {
            free(wacc_stmts[i]);
        }
        free(wacc_stmts);
        free(num_stmts_per_wacc);
    }

    int num_write_out_stmts = 5; // as defined in gen_write_out_code
    int count;
    for (l=0; l<nloops; l++) {
        Ploop *loop = loops[l];
        IF_DEBUG(printf("separating for loop: \n"););
        IF_DEBUG(pluto_loop_print(loop););

        /* Add statements to PlutoProg */
        sep_stmts[l] = (Stmt **)malloc(num_comm_stmts*num_data[l]*sizeof(Stmt *));
        count = 0;
        for (i=0; i<num_data[l]; i++)    {
            for (k=0; k<num_comm_stmts; k++) {
                pluto_add_given_stmt(prog, copy_comm_stmts[l][i][k]);
                sep_stmts[l][count++] = copy_comm_stmts[l][i][k];
            }
        }

        for (i=0; i<num_data[l]; i++) {
            for (k=0; k<num_write_out_stmts; k++) {
                pluto_add_given_stmt(prog, write_out_stmts[l][i][k]);
            }
        }
    }

    int scalar_dimensions_added[nloops*2];
    int num_scalar_dimensions_added = 0;

    for (l=0; l<nloops; l++) {
        // depth of the loop before any scalar dimensions were added
        if (options->blockcyclic) {
            dist_parallel_loop = loops[l]->depth + l;
        }else{
            dist_parallel_loop = loops[l]->depth;
        }

        // maintain the scalar dimensions added in sorted order
        for (i=0; i<num_scalar_dimensions_added; i++) {
            int orig_dimension = scalar_dimensions_added[i]-(i+1);
            // compare original dimensions before any scalar dimensinos were added
            if (orig_dimension>=dist_parallel_loop) {
                // update the depth of the loops after scalar dimensions were added by the previous loops
                // i represents the number of dimensions added before original depth of the loop
                dist_parallel_loop += i; 
                for (j=num_scalar_dimensions_added; j>i; j--) {
                    scalar_dimensions_added[j] = scalar_dimensions_added[j-1] + 1;
                }
                scalar_dimensions_added[i++] = dist_parallel_loop;
                num_scalar_dimensions_added++;
                break;
            }
        }
        if (i==num_scalar_dimensions_added) {
            // update the depth of the loops after scalar dimensions were added by the previous loops
            // i represents the number of dimensions added before original depth of the loop
            dist_parallel_loop += i;
            scalar_dimensions_added[i++] = dist_parallel_loop;
            num_scalar_dimensions_added++;
        }
        for (; i<num_scalar_dimensions_added; i++) {
            // compare current dimensions after scalar dimensions were added by the previous loops
            if (scalar_dimensions_added[i]>=(dist_parallel_loop+2)) {
                for (j=num_scalar_dimensions_added; j>i; j--) {
                    scalar_dimensions_added[j] = scalar_dimensions_added[j-1] + 1;
                }
                scalar_dimensions_added[i++] = dist_parallel_loop+2;
                num_scalar_dimensions_added++;
                break;
            }
        }
        if (i==num_scalar_dimensions_added) {
            scalar_dimensions_added[i++] = dist_parallel_loop+2;
            num_scalar_dimensions_added++;
        }

        /* Parallel loop is distributed (loop distribution) around these
         * statements - except for write-out statements */
        pluto_separate_stmts(prog, sep_stmts[l], num_comm_stmts*num_data[l], dist_parallel_loop);
        /* Add scalar dimension to separate the fused statements out inside (immediately inside); since they are
         * to be fused for dist_parallel_loop - used only for write-out pack and unpack */
        pluto_separate_stmts(prog, NULL, 0, dist_parallel_loop+2);

        /* dist_parallel_loop's position now changes */
        dist_parallel_loop++;

        FILE *outfp = fopen(".distmem", "w");
        if (outfp)  {
            fprintf(outfp, "t%d", dist_parallel_loop+1);
            fclose(outfp);
        }

        for (i=0; i<num_data[l]; i++)    {
            free(copy_comm_stmts[l][i]);
        }
        free(copy_comm_stmts[l]);
        free(sep_stmts[l]);
    } 

    for (l=0; l<nloops; l++) {
        // depth of the loop before any scalar dimensions were added
        if (options->blockcyclic) {
            dist_parallel_loop = loops[l]->depth + l;
        }else{
            dist_parallel_loop = loops[l]->depth;
        }

        // find the depth of the loop after scalar dimensions were added
        int num_scalar_dimensions_added_before = 0;
        for (i=0; i<num_scalar_dimensions_added; i++) {
            int orig_dimension = scalar_dimensions_added[i]-(i+1);
            // compare original dimensions before any scalar dimensinos were added
            if (orig_dimension>=dist_parallel_loop) {
                num_scalar_dimensions_added_before = i;
                dist_parallel_loop += i; 
                break;
            }
        }
        if (i==num_scalar_dimensions_added) {
            num_scalar_dimensions_added_before = i;
            dist_parallel_loop += i; 
        }

        /* Generate pi */
        generate_pi(pifp, dist_parallel_loop+1, prog, l, scalar_dimensions_added, num_scalar_dimensions_added_before);
    }
    fclose(pifp);

    // scalar dimensions have already been added before and after dist_parallel_loop
    /* Add outermost scalar dimension to separate the write-out statements from the rest */
    pluto_separate_stmts(prog, NULL, 0, 0);

    for (l=0; l<nloops; l++) {
        // depth of the loop before any scalar dimensions were added
        if (options->blockcyclic) {
            dist_parallel_loop = loops[l]->depth + l;
        }else{
            dist_parallel_loop = loops[l]->depth;
        }

        // find the depth of the loop after scalar dimensions were added
        for (i=0; i<num_scalar_dimensions_added; i++) {
            int orig_dimension = scalar_dimensions_added[i]-(i+1);
            // compare original dimensions before any scalar dimensinos were added
            if (orig_dimension>=dist_parallel_loop) {
                dist_parallel_loop += i; 
                break;
            }
        }
        if (i==num_scalar_dimensions_added) {
            dist_parallel_loop += i; 
        }

        dist_parallel_loop++; // due to the scalar dimension added to separate write-out statements from the rest

        int count = 0;
        for (i=0; i<num_data[l]; i++) {
            Stmt *stmt1, *stmt2;
            // write out stmts organized as follows (as defined gen_write_out_code)
            // pack-guard stmt, pack stmt, comm stmt, unpack-guard stmt, unpack stmt
            // separte write-out statements from the rest and
            // order them accordingly in the scalar dimension before dist_parallel_loop

            /* fuse write-out copy holder and text */
            stmt1 = write_out_stmts[l][i][0];
            stmt2 = write_out_stmts[l][i][1];
            stmt1->trans->val[0][stmt1->trans->ncols-1] = 1;
            stmt2->trans->val[0][stmt2->trans->ncols-1] = 1;
            stmt1->trans->val[dist_parallel_loop-1][stmt1->trans->ncols-1] = count;
            stmt2->trans->val[dist_parallel_loop-1][stmt2->trans->ncols-1] = count;
            /* separate them out immediately inside */
            stmt1->trans->val[dist_parallel_loop+1][stmt1->trans->ncols-1] = 0;
            stmt2->trans->val[dist_parallel_loop+1][stmt2->trans->ncols-1] = 1;
            count++;

            stmt1 = write_out_stmts[l][i][2];
            stmt1->trans->val[0][stmt1->trans->ncols-1] = 1;
            stmt1->trans->val[dist_parallel_loop-1][stmt1->trans->ncols-1] = count;
            count++;

            /* fuse write-out copy back holder and text */
            stmt1 = write_out_stmts[l][i][3];
            stmt2 = write_out_stmts[l][i][4];
            stmt1->trans->val[0][stmt1->trans->ncols-1] = 1;
            stmt2->trans->val[0][stmt2->trans->ncols-1] = 1;
            stmt1->trans->val[dist_parallel_loop-1][stmt1->trans->ncols-1] = count;
            stmt2->trans->val[dist_parallel_loop-1][stmt2->trans->ncols-1] = count;
            /* separate them out immediately inside */
            stmt1->trans->val[dist_parallel_loop+1][stmt1->trans->ncols-1] = 0;
            stmt2->trans->val[dist_parallel_loop+1][stmt2->trans->ncols-1] = 1;
            count++;
        }

        for (i=0; i<num_data[l]; i++)    {
            free(write_out_stmts[l][i]);
        }
        free(write_out_stmts[l]);
    }

    free(copy_comm_stmts);
    free(write_out_stmts);
    free(sep_stmts);

    free(num_data);
    free(copy_level);
    free(pi_mappings);

    pluto_mark_statements(prog);

    pluto_loops_free(loops,nloops);

    return (nloops==0);
}
