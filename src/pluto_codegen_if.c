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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution.
 *
 * Pluto codegen interface
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include <cloog/cloog.h>

#include "version.h"

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"
#include "program.h"
#include "ast_transform.h"

int num_common_tiles;

static int splitLoops(struct clast_stmt *s, int loop_level, int n, struct clast_stmt **split, 
        int stmt_num, CloogState *state, CloogDomain *domain);

static int get_first_point_loop(Stmt *stmt, const PlutoProg *prog)
{
    int i, first_point_loop;

    if (stmt->type != ORIG) {
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (!pluto_is_hyperplane_scalar(stmt, i)) {
                return i;
            }
        }
        /* No non-scalar hyperplanes */
        return 0;
    }

    for (i=stmt->last_tile_dim+1; i<stmt->trans->nrows; i++)   {
        if (stmt->hyp_types[i] == H_LOOP)  break;
    }

    if (i < prog->num_hyperplanes) {
        first_point_loop = i;
    }else{
        /* Should come here only if
         * it's a 0-d statement */
        first_point_loop = 0;
    }

    return first_point_loop;
}


/* Generate and print .cloog file from the transformations computed */
void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog)
{
    int i;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int npar = prog->npar;

    IF_DEBUG(printf("[Pluto] generating Cloog file...\n"));
    fprintf(fp, "# CLooG script generated automatically by PLUTO %s\n", PLUTO_VERSION);
    fprintf(fp, "# language: C\n");
    fprintf(fp, "c\n\n");

    /* Context: setting conditions on parameters */
    pluto_constraints_print_polylib(fp, prog->context);


    /* Setting parameter names */
    fprintf(fp, "\n1\n");
    for (i=0; i<npar; i++)  {
        fprintf(fp, "%s ", prog->params[i]);
    }
    fprintf(fp, "\n\n");

    fprintf(fp, "# Number of statements\n");
    fprintf(fp, "%d\n\n", nstmts);

    /* Print statement domains */
    for (i=0; i<nstmts; i++)    {
        fprintf(fp, "# S%d (%s)\n", stmts[i]->id+1, stmts[i]->text);
        fprintf(fp, "%d # of domains\n", 
                pluto_constraints_num_in_list(stmts[i]->domain));
        pluto_constraints_print_polylib(fp, stmts[i]->domain);
        fprintf(fp, "0 0 0\n\n");
    }

    fprintf(fp, "# we want cloog to set the iterator names\n");
    fprintf(fp, "0\n\n");

    fprintf(fp, "# of scattering functions\n");
    if (nstmts >= 1 && stmts[0]->trans != NULL) {
        fprintf(fp, "%d\n\n", nstmts);

        /* Print scattering functions */
        for (i=0; i<nstmts; i++) {
            fprintf(fp, "# T(S%d)\n", i+1);
            PlutoConstraints *sched = pluto_stmt_get_schedule(stmts[i]);
            pluto_constraints_print_polylib(fp, sched);
            fprintf(fp, "\n");
            pluto_constraints_free(sched);
        }

        /* Setting target loop names (all stmts have same number of hyperplanes */
        fprintf(fp, "# we will set the scattering dimension names\n");
        fprintf(fp, "%d\n", stmts[0]->trans->nrows);
        for (i=0; i<stmts[0]->trans->nrows; i++) {
            fprintf(fp, "t%d ", i+1);
        }
        fprintf(fp, "\n");
    }else{
        fprintf(fp, "0\n\n");
    }
}

static void gen_stmt_macro(const Stmt *stmt, FILE *outfp)
{
    int j;

    for (j=0; j<stmt->dim; j++) {
        if (stmt->iterators[j] == NULL) {
            printf("Iterator name not set for S%d; required \
                    for generating declarations\n", stmt->id+1);
            assert(0);
        }
    }
    fprintf(outfp, "#define S%d", stmt->id+1);
    fprintf(outfp, "(");
    for (j=0; j<stmt->dim; j++)  {
        if (j!=0)   fprintf(outfp, ",");
        fprintf(outfp, "%s", stmt->iterators[j]);
    }
    fprintf(outfp, ")\t");

    /* Generate pragmas for Bee/Cl@k */
    if (options->bee)   {
        fprintf(outfp, " __bee_schedule");
        for (j=0; j<stmt->trans->nrows; j++)    {
            fprintf(outfp, "[");
            pretty_print_affine_function(outfp, stmt->trans->val[j], 
                    stmt->dim, stmt->iterators);
            fprintf(outfp, "]");
        }
        fprintf(outfp, " _NL_DELIMIT_ ");
    }
    fprintf(outfp, "%s\n", stmt->text);
}


/* Generate variable declarations and macros */
int generate_declarations(const PlutoProg *prog, FILE *outfp)
{
    int i;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    /* Generate statement macros */
    for (i=0; i<nstmts; i++)    {
        gen_stmt_macro(stmts[i], outfp);
    }
    fprintf(outfp, "\n");

    /* Insert "dynsched" specific code
     * Number of parameters to task_computation is fixed as 4 in the static code - scheduler.h
     * So hardcode the parameters to the function as t1, t2, t3, t4
     */
    if (options->dynschedule) {

        fprintf(outfp,"\n\nvoid task_computation(int t1, int t2, int t3, int t4)\n{\n\n");
        fprintf(outfp, "\tint ");
        for (i = 4; i<stmts[0]->trans->nrows; i++)  {
            if (i!=4) fprintf(outfp, ", ");
            fprintf(outfp, "t%d", i+1);
        }
        fprintf(outfp, ";\n\n");

    }
    else {
        /* Scattering iterators. */
        if (prog->num_hyperplanes >= 1)    {
            fprintf(outfp, "\t\tint ");
            for (i=0; i<prog->num_hyperplanes; i++)  {
                if (i!=0) fprintf(outfp, ", ");
                fprintf(outfp, "t%d", i+1);
                if (prog->hProps[i].unroll)   {
                    fprintf(outfp, ", t%dt, newlb_t%d, newub_t%d", i+1, i+1, i+1);
                }
            }
            fprintf(outfp, ";\n\n");
        }
    }

    if (options->parallel || options->distmem)   {
        fprintf(outfp, "\tint lb, ub, lbp, ubp, lbd, ubd, lb2, ub2;\n");
    }
    /* For vectorizable loop bound replacement */
    fprintf(outfp, "\tregister int lbv, ubv;\n\n");

    return 0;
}


/* Call cloog and generate code for the transformed program
 *
 * cloogf, cloogl: set to -1 if you want the function to decide
 *
 * --cloogf, --cloogl overrides everything; next cloogf, cloogl if != -1,
 *  then the function takes care of the rest
 */
int pluto_gen_cloog_code(const PlutoProg *prog, int cloogf, int cloogl,
        FILE *cloogfp, FILE *outfp)
{
    CloogInput *input ;
    CloogOptions *cloogOptions ;
    CloogState *state;
    CloogDomain *domain;
    int i, j;

    struct clast_stmt *root,*split;
    split = NULL;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    state = cloog_state_malloc();
    cloogOptions = cloog_options_malloc(state);

    cloogOptions->fs = malloc (nstmts*sizeof(int));
    cloogOptions->ls = malloc(nstmts*sizeof(int));
    cloogOptions->fs_ls_size = nstmts;

    for (i=0; i<nstmts; i++) {
        cloogOptions->fs[i] = -1;
        cloogOptions->ls[i] = -1;
    }

    cloogOptions->name = "CLooG file produced by PLUTO";
    cloogOptions->compilable = 0;
    cloogOptions->esp = 1;
    cloogOptions->strides = 1;
    cloogOptions->quiet = options->silent;

    /* Generates better code in general */
    cloogOptions->backtrack = options->cloogbacktrack;

    if (options->cloogf >= 1 && options->cloogl >= 1) {
        cloogOptions->f = options->cloogf;
        cloogOptions->l = options->cloogl;
    }else{
        if (cloogf >= 1 && cloogl >= 1) {
            cloogOptions->f = cloogf;
            cloogOptions->l = cloogl;
        }else if (options->distmem) {

            /* We will set f/l statement-wise */

            int nploops = 0;
            Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);

            /* First set for compute statements */
            for (i=0; i<nploops; i++) {
                for (j=0; j<ploops[i]->nstmts; j++) {
                    Stmt *stmt = ploops[i]->stmts[j];
                    assert(stmt->type == ORIG);
                    cloogOptions->fs[ploops[i]->stmts[j]->id] = ploops[i]->depth+2;
                    cloogOptions->ls[i] = prog->num_hyperplanes;
                }
            }

            for (i=0; i<nstmts; i++) {
                if ((stmts[i]->type != FOIFI_COPY_OUT) || (stmts[i]->type != FOIFI_COPY_IN)) {
                    cloogOptions->fs[i] = prog->num_hyperplanes-1; // !!!roshan should find more accurate information
                    cloogOptions->ls[i] = prog->num_hyperplanes;
                }
                else if (stmts[i]->type != ORIG) {
                    assert(stmts[i]->parent_compute_stmt != NULL);
                    const Stmt *pstmt = stmts[i]->parent_compute_stmt;
                    cloogOptions->fs[i] = cloogOptions->fs[pstmt->id];
                    cloogOptions->ls[i] = prog->num_hyperplanes;
                }
            }

            /* Now refine the depths of the original statements that have been
             * set as well as set for those not set */
            for (i=0; i<nstmts; i++) {
                if (stmts[i]->type == ORIG) {
                    cloogOptions->fs[i] =  PLMAX(cloogOptions->fs[i],
                            get_first_point_loop(stmts[i], prog)+1);
                    cloogOptions->ls[i] = prog->num_hyperplanes;
                }
            }

        }else if (options->tile)   {
            for (i=0; i<nstmts; i++) {
                cloogOptions->fs[i] = get_first_point_loop(stmts[i], prog)+1;
                cloogOptions->ls[i] = prog->num_hyperplanes;
            }
        }else{
            /* Default */
            cloogOptions->f = 1;
            /* last level to optimize: number of hyperplanes;
             * since Pluto provides full-ranked transformations */
            cloogOptions->l = prog->num_hyperplanes;
        }
    }

    if (!options->silent)   {
        if (nstmts >= 1 && cloogOptions->fs[0] >= 1) {
            printf("[Pluto] using statement-wise -fs/-ls options: ");
            for (i=0; i<nstmts; i++) {
                printf("S%d(%d,%d), ", i+1, cloogOptions->fs[i], 
                        cloogOptions->ls[i]);
            }
            printf("\n");
        }else{
            printf("[Pluto] using Cloog -f/-l options: %d %d\n", 
                    cloogOptions->f, cloogOptions->l);
        }
    }

    if (options->cloogsh)
        cloogOptions->sh = 1;

    cloogOptions->name = "PLUTO-produced CLooG file";

    fprintf(outfp, "/* Start of CLooG code */\n");
    /* Get the code from CLooG */
    IF_DEBUG(printf("[Pluto] cloog_input_read\n"));
    input = cloog_input_read(cloogfp, cloogOptions) ;
    domain = cloog_domain_copy(input->context);
    IF_DEBUG(printf("[Pluto] cloog_clast_create\n"));
    root = cloog_clast_create_from_input(input, cloogOptions);
    if (options->prevector) {
        pluto_mark_vector(root, prog, cloogOptions);
    }
    if (options->dynschedule) {
        // statement ID doesn't matter since 'root' isn't printed
        splitLoops(root,0,num_common_tiles,&split,prog->nstmts,state,domain);
        assert(split != NULL);
        clast_pprint(outfp, split, 2, cloogOptions);
        cloog_clast_free(split);
    }else{
        if (options->parallel) {
            pluto_mark_parallel(root, prog, cloogOptions);
        }
        clast_pprint(outfp, root, 0, cloogOptions);
    }
    cloog_clast_free(root);

    fprintf(outfp, "/* End of CLooG code */\n");

    if(options->dynschedule) {
        fprintf(outfp,"\n\n}\n\n");
    }

    cloog_options_free(cloogOptions);
    cloog_state_free(state);

    return 0;
}


/* Generate code for a single multicore; the ploog script will insert openmp
 * pragmas later */
int pluto_multicore_codegen(FILE *cloogfp, FILE *outfp, const PlutoProg *prog)
{ 
    if (options->parallel)  {
        fprintf(outfp, "#include <omp.h>\n\n");
    }
    generate_declarations(prog, outfp);

    if (options->multipipe) {
        fprintf(outfp, "\tomp_set_nested(1);\n");
        fprintf(outfp, "\tomp_set_num_threads(2);\n");
    }

    pluto_gen_cloog_code(prog, -1, -1, cloogfp, outfp);

    return 0;
}


/* Decides which loops to mark parallel and generates the corresponding OpenMP
 * pragmas and writes them out to a file. They are later read by a script
 * (ploog) and appropriately inserted into the output Cloog code
 *
 * Returns: the number of parallel loops for which OpenMP pragmas were generated 
 *
 * Generate the #pragma comment -- will be used by a syntactic scanner
 * to put in place -- should implement this with CLast in future */
int pluto_omp_parallelize(PlutoProg *prog)
{
    int i;

    FILE *outfp = fopen(".pragmas", "w");

    if (!outfp) return 1;

    HyperplaneProperties *hProps = prog->hProps;

    int loop;

    /* IMPORTANT: Note that by the time this function is called, pipelined
     * parallelism has already been converted to inner parallelism in
     * tile space (due to a tile schedule) - so we don't need check any
     * PIPE_PARALLEL properties
     */
    /* Detect the outermost sync-free parallel loop - find upto two of them if
     * the multipipe option is set */
    int num_parallel_loops = 0;
    for (loop=0; loop<prog->num_hyperplanes; loop++) {
        if (hProps[loop].dep_prop == PARALLEL && hProps[loop].type != H_SCALAR)   {
            // Remember our loops are 1-indexed (t1, t2, ...)
            fprintf(outfp, "t%d #pragma omp parallel for shared(", loop+1);

            for (i=0; i<loop; i++)  {
                fprintf(outfp, "t%d,", i+1);
            }

            for (i=0; i<num_parallel_loops+1; i++) {
                if (i!=0) fprintf(outfp, ",");
                fprintf(outfp,  "lb%d,ub%d", i+1, i+1);
            }

            fprintf(outfp,  ") private(");

            if (options->prevector) {
                fprintf(outfp,  "ubv,lbv,");
            }

            /* Lower and upper scalars for parallel loops yet to be marked */
            /* NOTE: we extract up to 2 degrees of parallelism
            */
            if (options->multipipe) {
                for (i=num_parallel_loops+1; i<2; i++) {
                    fprintf(outfp,  "lb%d,ub%d,", i+1, i+1);
                }
            }

            for (i=loop; i<prog->num_hyperplanes; i++)  {
                if (i!=loop) fprintf(outfp, ",");
                fprintf(outfp, "t%d", i+1);
            }
            fprintf(outfp, ")\n");

            num_parallel_loops++;

            if (!options->multipipe || num_parallel_loops == 2)   {
                break;
            }
        }
    }

    IF_DEBUG(fprintf(stdout, "[Pluto] marked %d loop(s) parallel\n",
                num_parallel_loops));

    fclose(outfp);

    return num_parallel_loops;
}


// s should be root when called
// !!! only the first occurence of a set of statements at loop level n is split 
// !!! what if there are multiple occurrences of sets of statements at loop level n? e.g. when the loop is split
static int splitLoops(struct clast_stmt *s,int loop_level,int n,
        struct clast_stmt **split, int stmt_num, CloogState *state, 
        CloogDomain *domain)
{
    struct clast_user_stmt *u = NULL;
    CloogStatement *cs = NULL;

    for ( ; s ; s = s->next){
        if (CLAST_STMT_IS_A(s, stmt_for))  {
            struct clast_for *for_s = (struct clast_for *)s;
            if ( loop_level+1 == n )    {
                cs = cloog_statement_alloc(state, stmt_num);
                u = new_clast_user_stmt(domain,cs,NULL); // add a user defined statement
                cloog_statement_free(cs);
                *split = for_s->body;
                for_s->body = &u->stmt;
                break;
            }
            splitLoops(for_s->body, loop_level+1, n, split, stmt_num, state, domain);
        }
        else if (CLAST_STMT_IS_A(s, stmt_guard)) {
            splitLoops(((struct clast_guard *)s)->then, loop_level, n, split, stmt_num, state, domain);
            if (*split != NULL)
                break;
        }
        else if (CLAST_STMT_IS_A(s, stmt_block)) {
            splitLoops(((struct clast_block *)s)->body, loop_level, n, split, stmt_num, state, domain);
            if (*split != NULL)
                break;
        }
    } // end of while
    return 0;
}

void print_dynsched_file(char *srcFileName, FILE *cloogfp, FILE *outfp, PlutoProg* prog)
{
    int i, j, k;
    int srcid, destid, startIteratorIndex, noOfIterators, common_band_num, stmt_common_tile_cnt;
    char sysloogFileName[256];
    char tempString[256], variablesString[256], argumentsString[256];
    FILE *sysloogfp;
    PlutoConstraints *transformedDpolytope;
    Stmt *srcStmt, *destStmt;

    CloogInput *input;
    CloogOptions *cloogOptions;
    CloogState *state;
    CloogDomain *domain;
    struct clast_stmt *root,*split;
    split = NULL;

    Stmt **stmts = prog->stmts;
    Dep **deps = prog->deps;
    int nstmts = prog->nstmts;
    int ndeps = prog->ndeps;
    int npar = prog->npar;
    char **params = prog->params;
    HyperplaneProperties *hProps = prog->hProps;
    int num_hyperplanes = prog->num_hyperplanes;

    /* Number of common tiles is taken as number of common iterators (after transformation) 
     *  in the first non-sequential band.
     *  However the clast functions have assumed that the common tiles start from the first band */
    num_common_tiles = prog->nvar;
    for (j=0;j<num_hyperplanes;j++) {
        if (hProps[j].dep_prop != SEQ) break;
    }
    common_band_num = hProps[j].band_num;
    for (i=0;i<nstmts;i++) {
        stmt_common_tile_cnt = 0;
        j=0;
        while ( (hProps[j].band_num != common_band_num) && (j<num_hyperplanes) ) j++;	
        while ( (hProps[j].band_num == common_band_num) && (j<num_hyperplanes) ) { j++; stmt_common_tile_cnt++; }	
        if (stmt_common_tile_cnt < num_common_tiles)
            num_common_tiles = stmt_common_tile_cnt;
    }
    IF_DEBUG(printf("No of variables = %d , No of common tiles = %d\n", prog->nvar, num_common_tiles));

    state = cloog_state_malloc();
    cloogOptions = cloog_options_malloc(state);

    cloogOptions->name = "CLooG file produced by PLUTO for dynamic scheduling";
    cloogOptions->compilable = 0;
    cloogOptions->esp = 0; // !!!roshan not sure if this needs to be enforced - can 1 be used?
    cloogOptions->quiet = options->silent;
    cloogOptions->backtrack = 1; /* Generates better code in general */

    /* Generate a file that has a function to create DAG
     * 
     * Generate tiled code to split the common tiled loops
     * to enumerate all vertices
     *
     * The loop is written twice, once to estimate the number of vertices
     * and once to actually add the vertices
     */

    /* Assume that the common tile iterators would be named 't1, t2...'*/
    /* Temporarily dag_add_vertex and dag_add_edge assumes 4-tuple vertex in scheduler.h */
    // build the string for intra-tile iterators, assuming the character 't' for it
    strcpy(tempString,"");
    for (j=1;j<num_common_tiles;j++) 
        sprintf(tempString,"%s t%d,", tempString, j);
    sprintf(tempString,"%s t%d",tempString, j);
    strcpy(variablesString, tempString);
    strcat(variablesString, ", ");
    // append 0 for the arguments to satisfy the static 4-tuple vertex in scheduler.h
    /* Careful - temporary stuff*/
    for (j=num_common_tiles+1;j<=4;j++)	
        strcat(tempString,", 0");
    strcpy(argumentsString, tempString);
    strcat(argumentsString, ", ");

    fprintf(outfp,"#define S%d() {dag_add_vertex(%s, 1, &vId);}\n",3,tempString);

    // build the string for inter-tile iterators, assuming the character 's' for it
    strcpy(tempString,"");
    for (j=1;j<num_common_tiles;j++) 
        sprintf(tempString,"%s s%d,", tempString, j);
    sprintf(tempString,"%s s%d",tempString, j);
    strcat(variablesString, tempString);
    // append 0 for the arguments to satisfy the static 4-tuple vertex in scheduler.h
    /* Careful - temporary stuff*/
    for (j=num_common_tiles+1;j<=4;j++)	
        strcat(tempString,", 0");
    strcat(argumentsString, tempString);

    // add_edge statement has ID 1
    // if there is an issue due to this ID overlapping with that of the original code statement,
    // then it can be changed to nstmts+3
    fprintf(outfp,"#define S%d(%s) {dag_add_edge(%s, 0.0);}\n",1,variablesString,argumentsString);
    fprintf(outfp,"#define S%d() {++vertexCount;}\n",2);
    fprintf(outfp,"\nvoid generate_dag()\n{\n");
    fprintf(outfp,"  int %s, vId;\n", variablesString);
    fprintf(outfp,"  int vertexCount=0;\n\n");

    if (options->ft == -1)  {
        if (stmts[0]->num_tiled_loops < 4)   {
            cloogOptions->f = stmts[0]->num_tiled_loops+1;
            cloogOptions->l = stmts[0]->trans->nrows;
        }else{
            cloogOptions->f = stmts[0]->num_tiled_loops+1;
            cloogOptions->l = stmts[0]->trans->nrows;
        }
    }else{
        cloogOptions->f = stmts[0]->num_tiled_loops+options->ft+1;
        cloogOptions->l = stmts[0]->trans->nrows;
    }

    input = cloog_input_read(cloogfp,cloogOptions) ;
    domain = cloog_domain_copy(input->context);
    root = cloog_clast_create_from_input(input, cloogOptions);

    // generate code to estimate the number of vertices
    splitLoops(root,0,num_common_tiles,&split,2, state, domain);
    clast_pprint(outfp, root, 2, cloogOptions);

    assert(split != NULL);
    cloog_clast_free(split);
    split = NULL;

    // initialize the DAG with the estimated number of vertices
    fprintf(outfp,"\n  init_dag(vertexCount);\n\n");

    // generate code to add the vertices
    splitLoops(root,0,num_common_tiles,&split,3, state, domain);
    clast_pprint(outfp, root, 2, cloogOptions);

    assert(split != NULL);
    cloog_clast_free(split);
    cloog_clast_free(root);

    fprintf(outfp, "\n");

    /* This is the core part
     *  Generate inter-tile dependencies (cloog file)
     *  Create a system to capture the dependencies in the transformed tiled space:
     *          - source equalities/inequaities
     *              1. original domain + tile definition (domain in tiled.cloog)
     *              2. transformations (scattering functions in tiled.cloog)
     *          - destination equalities/inequaities
     *              1. original domain + tile definition (domain in tiled.cloog)
     *              2. transformations (scattering functions in tiled.cloog)
     *          - equalities of original h-transformation of dependence
     */
    for (i=0;i<ndeps;i++) {
        /* RAW, WAR, WAW deps matter */
        /*if(!IS_RAR(deps[i]->type)) */
        // !!!roshan not sure if WAR and WAW dependences need to be considered.
        // If they are considered, it leads to slower execution time sometimes (due to increase in size of the DAG?),
        // longer compilation time for generated code of some benchmarks like heat-3d
        // and incorrect results for some benchmarks like heat-2d
        if(deps[i]->type == CANDL_RAW) {
            strcpy(sysloogFileName, srcFileName);
            sysloogFileName[strlen(srcFileName)-2] = '\0';
            sprintf(sysloogFileName, "%s.dep%d.dynsched.sysloog",sysloogFileName,i);

            // Generate the cloog file to generate dag code corresponding to the dep
            sysloogfp = fopen(sysloogFileName, "w+");
            if (!sysloogfp) {
                fprintf(stderr, "Can't open file %s for writing\n", sysloogFileName);
                return;
            }

            fprintf(sysloogfp, "# language: C\n");
            fprintf(sysloogfp, "c\n\n");

            fprintf(sysloogfp, "# Context\n");
            fprintf(sysloogfp, "%d %d\n", npar, npar+2);

            for (j=0; j<npar; j++)  {
                fprintf(sysloogfp, "1 ");
                for (k=0; k<npar; k++)  {
                    if (j == k) {
                        fprintf(sysloogfp, "1 ");
                    }
                    else fprintf(sysloogfp, "0 ");
                }
                fprintf(sysloogfp, "%d\n", -options->context);
            }

            fprintf(sysloogfp, "\n1\n");
            fprintf(sysloogfp, "# Parameter name(s)\n");
            for (j=0; j<npar; j++)  {
                fprintf(sysloogfp, "%s ", params[j]);
            }
            fprintf(sysloogfp, "\n\n");

            fprintf(sysloogfp, "# Number of statements\n");
            fprintf(sysloogfp, "1\n\n");

            fprintf(sysloogfp, "# Iteration Domain\n");
            fprintf(sysloogfp, "1\n\n");

            srcid = deps[i]->src;
            destid = deps[i]->dest;
            srcStmt = stmts[srcid];
            destStmt = stmts[destid];

            // get the transformed dependence polytope, which contains the scattering functions
            transformedDpolytope = pluto_get_transformed_dpoly(deps[i], srcStmt, destStmt);

            assert(transformedDpolytope->ncols == (2 * num_hyperplanes + npar + 1));

            // project out the destination statement intra-tile iterators
            startIteratorIndex = num_hyperplanes + num_common_tiles; 
            noOfIterators = 2 * num_hyperplanes - startIteratorIndex;
            pluto_constraints_project_out(transformedDpolytope, startIteratorIndex, noOfIterators);

            // project out the source statement intra-tile iterators
            startIteratorIndex = num_common_tiles; 
            noOfIterators = num_hyperplanes - startIteratorIndex;
            pluto_constraints_project_out(transformedDpolytope, startIteratorIndex, noOfIterators);

            pluto_constraints_print_polylib(sysloogfp, transformedDpolytope);

            fprintf(sysloogfp, "\n0 0 0\n\n");
            fprintf(sysloogfp, "1\n");
            fprintf(sysloogfp, "# Iterator name(s)\n\n");
            // inter-tile iterators use the character 's', according to assumption above while generating the variable declarations
            for (j=1; j<=num_common_tiles; j++)
                fprintf(sysloogfp, "s%d ",j);
            // intra-tile iterators use the character 't', according to assumption above while generating the variable declarations
            for (j=1; j<=num_common_tiles; j++)
                fprintf(sysloogfp, "t%d ",j);
            fprintf(sysloogfp, "\n\n");
            fprintf(sysloogfp, "# Scattering functions\n");
            fprintf(sysloogfp, "0\n\n");

            rewind(sysloogfp);

            cloogOptions->f = 1;
            cloogOptions->l = num_common_tiles;

            input = cloog_input_read(sysloogfp,cloogOptions) ;
            root = cloog_clast_create_from_input(input, cloogOptions);
            // assuming only one (innermost) statement, whose ID will be 1
            // if it is necessary to change the ID to nstmts+3, 
            // the clast statements should be traversed to find the only clast user statement
            // assuming statement ID of 1 will work, even though the ID overlaps with one of the original code statements
            clast_pprint(outfp, root, 2, cloogOptions);
            cloog_clast_free(root);
            fclose(sysloogfp);
        }
    }

    fprintf(outfp,"\n  update_number_vertices();");
    fprintf(outfp,"\n\n}\n\n");
    fprintf(outfp,"#undef S1\n");
    fprintf(outfp,"#undef S2\n");
    fprintf(outfp,"#undef S3\n\n\n");

    cloog_options_free(cloogOptions);
}
