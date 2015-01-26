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
    PlutoConstraints *ctx = pluto_constraints_dup(prog->context);
    pluto_constraints_intersect(ctx, prog->codegen_context);
    pluto_constraints_print_polylib(fp, ctx);
    pluto_constraints_free(ctx);

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
            pluto_affine_function_print(outfp, stmt->trans->val[j], 
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

    if (options->parallel)   {
        fprintf(outfp, "\tint lb, ub, lbp, ubp, lb2, ub2;\n");
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
    int i;

    struct clast_stmt *root;

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
    IF_DEBUG(printf("[Pluto] cloog_clast_create\n"));
    root = cloog_clast_create_from_input(input, cloogOptions);
    if (options->prevector) {
        pluto_mark_vector(root, prog, cloogOptions);
    }
    if (options->parallel) {
        pluto_mark_parallel(root, prog, cloogOptions);
    }
    clast_pprint(outfp, root, 0, cloogOptions);
    cloog_clast_free(root);

    fprintf(outfp, "/* End of CLooG code */\n");

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
