/*
 * PLUTO: A automatic parallelizer + locality optimizer (experimental)
 * 
 * Copyright (C) 2007 - 2008 Uday Kumar Bondhugula
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
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pluto.h"

#include "clan/clan.h"
#include "candl/candl.h"

#include "math_support.h"
#include "post_transform.h"
#include "ddg.h"
#include "program.h"

/* Global variables :-( */
int npar;
int nvar;
PlutoOptions *options;

void usage_message();

int main(int argc, char *argv[])
{
    int i;

    FILE *src_fp;

    int option;
    int option_index = 0;

    char srcFileName[256];
    char outFileName[256] = "";

    char cloogFileName[256];
    FILE *cloogfp, *outfp;

    if (argc <= 1)  {
        usage_message();
        return 1;
    }

    options = pluto_options_alloc();

    const struct option pluto_options[] =
    {
        {"tile", no_argument, &options->tile, 1},
        {"notile", no_argument, &options->tile, 0},
        {"debug", no_argument, &options->debug, true},
        {"moredebug", no_argument, &options->moredebug, true},
        {"rar", no_argument, &options->rar, 1},
        {"nofuse", no_argument, &options->fuse, NO_FUSE},
        {"maxfuse", no_argument, &options->fuse, MAXIMAL_FUSE},
        {"smartfuse", no_argument, &options->fuse, SMART_FUSE},
        {"parallel", no_argument, &options->parallel, 1},
        {"parallelize", no_argument, &options->parallel, 1},
        {"unroll", no_argument, &options->unroll, 1},
        {"nounroll", no_argument, &options->unroll, 0},
        {"polyunroll", no_argument, &options->polyunroll, 1},
        {"bee", no_argument, &options->bee, 1},
        {"ufactor", required_argument, 0, 'u'},
        {"prevector", no_argument, &options->prevector, 1},
        {"noprevector", no_argument, &options->prevector, 0},
        {"context", required_argument, 0, 'c'},
        {"cloogf", required_argument, 0, 'F'},
        {"cloogl", required_argument, 0, 'L'},
        {"ft", required_argument, 0, 'f'},
        {"lt", required_argument, 0, 'l'},
        {"multipipe", no_argument, &options->multipipe, 1},
        {"l2tile", no_argument, &options->l2tile, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"indent", no_argument, 0, 'i'},
        {"silent", no_argument, &options->silent, 1},
        {"lastwriter", no_argument, &options->lastwriter, 1},
        {"nobound", no_argument, &options->nobound, 1},
        {"scalpriv", no_argument, &options->scalpriv, 1},
        {0, 0, 0, 0}
    };


    /* Read command-line options */
    while (1) {
        option = getopt_long(argc, argv, "bhiqvf:l:F:L:c:", pluto_options,
                &option_index);

        if (option == -1)   {
            break;
        }

        switch (option) {
            case 0:
                break;
            case 'F':
                options->cloogf = atoi(optarg);
                break;
            case 'L':
                options->cloogl = atoi(optarg);
                break;
            case 'b':
                options->bee = 1;
                break;
            case 'c':
                options->context = atoi(optarg);
                break;
            case 'd':
                break;
            case 'f':
                options->ft = atoi(optarg);
                break;
            case 'g':
                break;
            case 'h':
                usage_message();
                return 1;
            case 'i':
                /* Handled in polycc */
                break;
            case 'l':
                options->lt = atoi(optarg);
                break;
            case 'm':
                break;
            case 'n':
                break;
            case 's':
                break;
            case 'p':
                break;
            case 'q':
                options->silent = 1;
                break;
            case 'u':
                options->ufactor = atoi(optarg);
                break;
            case 'v':
                printf("PLUTO 0.5.0 - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2008  Uday Kumar Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
                return 1;
            default:
                usage_message();
                return 2;
        }
    }


    if (optind <= argc-1)   {
        strncpy(srcFileName, argv[optind], 250);
    }else{
        /* No non-option argument was specified */
        usage_message();
        return 3;
    }

    src_fp  = fopen(srcFileName, "r");

    if (!src_fp)   {
        fprintf(stderr, "pluto: error opening source file: '%s'\n", srcFileName);
        return 5;
    }

    /* Extract polyhedral representation from input program */
    scoplib_scop_p scop;

    clan_options_p clanOptions = clan_options_malloc();

    scop = clan_scop_extract(src_fp, clanOptions);

    if (!scop->statement)   {
        fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                srcFileName);
        return 1;
    }

    /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */

    /* Convert clan scop to Pluto program */
    PlutoProg *prog = scop_to_pluto_prog(scop, options);

    clan_options_free(clanOptions);

    /* Backup irregular program portion in .scop. */
    char* irroption = scoplib_scop_tag_content(scop, "<irregular>",
                                            "</irregular>");

    scoplib_scop_free(scop);

    IF_DEBUG2(deps_print(stdout, prog->deps, prog->ndeps));
    IF_DEBUG2(stmts_print(stdout, prog->stmts, prog->nstmts));

    /* Create the data dependence graph */
    prog->ddg = ddg_create(prog);
    ddg_compute_scc(prog);

    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i].dim;
    }

    /* Make options consistent */
    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

    /* Disable pre-vectorization if tile is not on */
    if (options->tile == 0 && options->prevector == 1) {
        /* If code will not be tiled, pre-vectorization does not make
         * sense */
        if (!options->silent)   {
            fprintf(stdout, "[Pluto] Warning: pre-vectorization does not fit (--tile is off)\n");
        }
        options->prevector = 0;
    }

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Number of statements: %d\n", prog->nstmts);
        fprintf(stdout, "[Pluto] Total number of loops: %d\n", dim_sum);
        fprintf(stdout, "[Pluto] Number of deps: %d\n", prog->ndeps);
        fprintf(stdout, "[Pluto] Maximum domain dimensionality: %d\n", nvar);
        fprintf(stdout, "[Pluto] Number of parameters: %d\n", npar);
    }

    /* Auto transformation */
    pluto_auto_transform(prog);

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Affine transformations [<iter coeff's> <const>]\n\n");
    }

    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;

    /* Print out the transformations */
    if (!options->silent)   {
        for (i=0; i<nstmts; i++) {
            fprintf(stdout, "T(S%d): ", i+1);
            int level;
            printf("(");
            for (level=0; level<prog->num_hyperplanes; level++) {
                if (level > 0) printf(", ");
                pretty_print_affine_function(stdout, &stmts[i], level);
            }
            printf(")\n");

            pluto_matrix_print(stdout, stmts[i].trans);
        }

        print_hyperplane_properties(prog->hProps, prog->num_hyperplanes);
    }

    if (options->tile)   {
        pluto_tile(prog);
    }

    if (options->parallel)   {
        int outermostBandStart, outermostBandEnd;
        getOutermostTilableBand(prog, &outermostBandStart, &outermostBandEnd);

        /* Obtain pipelined parallelization by skewing the tile space */
        bool retval = create_tile_schedule(prog, outermostBandStart, outermostBandEnd);

        /* Even if the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user, but anyway do fine-grained 
         * parallelization
         */
        if (retval && options->tile == 0)   {
            printf("WARNING: --tile is not used and there is pipelined parallelism\n");
            printf("\t This leads to finer grained parallelism; add --tile to the list\n");
            printf("\t of cmd-line options for a better coarse-grained parallelized code.\n");
        }
    }

    if (options->prevector) {
        pre_vectorize(prog);
    }else{
        /* Create an empty .vectorize file */
        fopen(".vectorize", "w");
    }

    if (options->tile && !options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        print_hyperplane_properties(prog->hProps, prog->num_hyperplanes);
    }

    if (options->parallel)  {
        /* Generate meta info for insertion of OpenMP pragmas */
        generate_openmp_pragmas(prog);
    }


    if (options->unroll || options->polyunroll)    {
        /* Will generate a .unroll file */
        /* plann needs a .params */
        FILE *paramsFP = fopen(".params", "w");
        if (paramsFP)   {
            int i;
            for (i=0; i<npar; i++)  {
                fprintf(paramsFP, "%s\n", prog->params[i]);
            }
            fclose(paramsFP);
        }
        detect_unrollable_loops(prog);
    }else{
        /* Create an empty .unroll file */
        fopen(".unroll", "w");
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }

    /* The .cloog file name */
    strcpy(cloogFileName, srcFileName);
    cloogFileName[strlen(srcFileName)-2] = '\0';

    if (options->parallel && options->multipipe)   {
        strcat(cloogFileName, ".par2d.cloog");
    }else if (options->parallel)   {
        strcat(cloogFileName, ".par.cloog");
    }else if (options->tile)  {
        strcat(cloogFileName, ".tiled.cloog");
    }else{
        strcat(cloogFileName, ".opt.cloog");
    }

    cloogfp = fopen(cloogFileName, "w+");

    /* Remove the .c extension and append a new one */
    strcpy(outFileName, srcFileName);
    outFileName[strlen(srcFileName)-2] = '\0';
    strcat(outFileName, ".pluto.c");

    outfp = fopen(outFileName, "w");

    if (!cloogfp)   {
        fprintf(stderr, "Can't open .cloog file: %s\n", cloogFileName);
        return 2;
    }


    /* Generate the .cloog file */
    IF_DEBUG(printf("[Pluto] Generating Cloog file\n"));
    print_cloog_file(cloogfp, prog);
    /* Add the <irregular> tag from clan, if any. */
    if (irroption != NULL)
    {
        fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
        free(irroption);
    }
    rewind(cloogfp);

    if (!outfp) {
        fprintf(stderr, "Can't open file %s for writing\n", outFileName);
        return 1;
    }

    /* Generate code using Cloog and add necessary stuff before/after code */
    pluto_codegen(cloogfp, outfp, prog);

    fclose(cloogfp);

    pluto_options_free(options);

    pluto_prog_free(prog);

    return 0;
}


void usage_message(void)
{
    fprintf(stdout, "Usage: polycc <input.c> [options]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --tile                 Tile for locality\n");
    fprintf(stdout, "       --parallel             Automatically parallelize using OpenMP pragmas\n");
    fprintf(stdout, "       | --parallelize\n");
    fprintf(stdout, "       --l2tile               Tile a second time (typically for L2 cache) - disabled by default \n");
    fprintf(stdout, "       --multipipe            Extract two degrees of pipelined parallelism if possible;\n");
    fprintf(stdout, "                                 by default one degree is extracted (if it exists)\n");
    fprintf(stdout, "       --rar                  Consider RAR dependences too (disabled by default)\n");
    fprintf(stdout, "       --[no]unroll           Unroll-jam (disabled by default)\n");
    fprintf(stdout, "       --ufactor=<factor>     Unroll-jam factor (default is 8)\n");
    fprintf(stdout, "       --[no]prevector        Make code amenable to compiler auto-vectorization (with ICC) - enabled by default\n");
    fprintf(stdout, "       --context=<context>    Parameters are at least as much as <context>\n");
    fprintf(stdout, "       --bee                  Generate pragmas for Bee+Cl@k\n\n");
    fprintf(stdout, "       --indent  | -i         Indent generated code (disabled by default)\n");
    fprintf(stdout, "       --silent  | -q         Silent mode; no output as long as everything goes fine (disabled by default)\n");
    fprintf(stdout, "       --help    | -h         Print this help menu\n");
    fprintf(stdout, "       --version | -v         Display version number\n");
    fprintf(stdout, "\n   Fusion                Options to control fusion heuristic\n");
    fprintf(stdout, "       --nofuse               Do not fuse across SCCs of data dependence graph\n");
    fprintf(stdout, "       --maxfuse              Maximal fusion\n");
    fprintf(stdout, "       --smartfuse [default]  Heuristic (in between nofuse and maxfuse)\n");
    fprintf(stdout, "\n   Debugging\n");
    fprintf(stdout, "       --debug        Verbose output\n");
    fprintf(stdout, "\nTo report bugs, please send an email to the author <bondhugula.1@osu.edu>\n\n");
}
