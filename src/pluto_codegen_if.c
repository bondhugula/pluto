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

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"
#include "program.h"
#include "version.h"


/* Generate and print .cloog file from the transformations computed */
void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog)
{
    int i;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int npar = prog->npar;

    IF_DEBUG(printf("[Pluto] Generating Cloog file\n"));
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
        fprintf(fp, "%d # of domains\n", 1);
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
    fprintf(outfp, "\t#define S%d", stmt->id+1);
    fprintf(outfp, "(");
    for (j=0; j<stmt->dim; j++)  {
        if (j!=0)   fprintf(outfp, ",");
        fprintf(outfp, "%s", stmt->iterators[j]);
    }
    fprintf(outfp, ")\t");

    /* Generate pragmas for Bee/Cl@k */
    if (options->bee)   {
        fprintf(outfp, " schedule");
        for (j=0; j<stmt->trans->nrows; j++)    {
            fprintf(outfp, "[");
            pretty_print_affine_function(outfp, stmt, j);
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
        fprintf(outfp, "\tregister int lb, ub, lb1, ub1, lb2, ub2;\n");
    }
    if (options->prevector)   {
        /* For vectorizable loop bound replacement */
        fprintf(outfp, "\tregister int lbv, ubv;\n\n");
    }

    return 0;
}

/* Call cloog and generate code for the transformed program */
int pluto_gen_cloog_code(const PlutoProg *prog, FILE *cloogfp, FILE *outfp)
{
    CloogProgram *program;
    CloogOptions *cloogOptions ;
    CloogState *state;

    Stmt **stmts = prog->stmts;

    state = cloog_state_malloc();
    cloogOptions = cloog_options_malloc(state);

    cloogOptions->name = "CLooG file produced by PLUTO";
    cloogOptions->compilable = 0;
    cloogOptions->esp = 1;
    cloogOptions->strides = 1;
    cloogOptions->quiet = options->silent;

    /* Generates better code in general */
    cloogOptions->backtrack = 1;

    if (options->cloogf >= 1 && options->cloogl >= 1) {
        cloogOptions->f = options->cloogf;
        cloogOptions->l = options->cloogl;
    }else{
        if (options->tile && stmts[0]->trans != NULL)   {
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
        }else{
            /* Default */
            cloogOptions->f = 1;
            /* last level to optimize: infinity */
            cloogOptions->l = -1;
        }
    }

    if (!options->silent)   {
        printf("[Pluto] using Cloog -f/-l options: %d %d\n", cloogOptions->f, cloogOptions->l);
    }

    cloogOptions->name = "PLUTO-produced CLooG file";

    /* Get the code from CLooG */
    IF_DEBUG(printf("[Pluto] cloog_program_read \n"));
    program = cloog_program_read(cloogfp, cloogOptions) ;
    IF_DEBUG(printf("[Pluto] cloog_program_generate \n"));
    program = cloog_program_generate(program,cloogOptions) ;
    cloog_program_pprint(outfp, program, cloogOptions) ;

    fprintf(outfp, "/* End of CLooG code */\n");

    cloog_options_free(cloogOptions);
    cloog_program_free(program);
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

    pluto_gen_cloog_code(prog, cloogfp, outfp);

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

    IF_DEBUG(fprintf(stdout, "[Pluto] marked %d loop(s) parallel\n", num_parallel_loops));

    fclose(outfp);

    return num_parallel_loops;
}

