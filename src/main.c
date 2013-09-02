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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pluto.h"
#include "transforms.h"
#include "math_support.h"
#include "post_transform.h"
#include "program.h"
#include "version.h"

#include "clan/clan.h"
#include "candl/candl.h"

PlutoOptions *options;

void usage_message(void)
{
    fprintf(stdout, "Usage: polycc <input.c> [options] [-o output]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --tile               Tile for locality\n");
    fprintf(stdout, "       --intratileopt       Optimize intra-tile execution order for locality\n");
    fprintf(stdout, "       --parallel             Automatically parallelize using OpenMP pragmas\n");
    fprintf(stdout, "       | --parallelize\n");
    fprintf(stdout, "       --l2tile               Tile a second time (typically for L2 cache) - disabled by default \n");
    fprintf(stdout, "       --multipipe            Extract two degrees of pipelined parallelism if possible;\n");
    fprintf(stdout, "       --distmem            Parallelize for distributed memory (generate MPI)\n");
    fprintf(stdout, "                                 by default one degree is extracted (if it exists)\n");
    fprintf(stdout, "       --rar                  Consider RAR dependences too (disabled by default)\n");
    fprintf(stdout, "       --[no]unroll           Unroll-jam (disabled by default)\n");
    fprintf(stdout, "       --ufactor=<factor>     Unroll-jam factor (default is 8)\n");
    fprintf(stdout, "       --[no]prevector        Make code amenable to compiler auto-vectorization (with ICC) - enabled by default\n");
    fprintf(stdout, "       --context=<context>    Parameters are at least as much as <context>\n");
    fprintf(stdout, "       --forceparallel=<bitvec>  6 bit-vector of depths (1-indexed) to force parallel (0th bit represents depth 1)\n");
    fprintf(stdout, "       --isldep               Use ISL-based dependence tester\n");
    fprintf(stdout, "       --islsolve             Use ISL as ilp solver\n");
    fprintf(stdout, "       --readscoplib          Read input from a scoplib file\n");
    fprintf(stdout, "       --lastwriter           Work with refined dependences (last conflicting access is computed for RAW/WAW)\n");
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
    fprintf(stdout, "       --moredebug    More verbose output\n");
    fprintf(stdout, "\nTo report bugs, please send an email to <pluto-development@googlegroups.com>\n\n");
}

int main(int argc, char *argv[])
{
    int i;

    FILE *src_fp;

    int option;
    int option_index = 0;

    char *srcFileName;

    FILE *cloogfp, *outfp, *dynschedfp;

    dynschedfp = NULL;

    if (argc <= 1)  {
        usage_message();
        return 1;
    }

    options = pluto_options_alloc();

    const struct option pluto_options[] =
    {
        {"tile", no_argument, &options->tile, 1},
        {"notile", no_argument, &options->tile, 0},
        {"intratileopt", no_argument, &options->intratileopt, 1},
        {"debug", no_argument, &options->debug, true},
        {"moredebug", no_argument, &options->moredebug, true},
        {"rar", no_argument, &options->rar, 1},
        {"identity", no_argument, &options->identity, 1},
        {"nofuse", no_argument, &options->fuse, NO_FUSE},
        {"maxfuse", no_argument, &options->fuse, MAXIMAL_FUSE},
        {"smartfuse", no_argument, &options->fuse, SMART_FUSE},
        {"parallel", no_argument, &options->parallel, 1},
        {"parallelize", no_argument, &options->parallel, 1},
        {"innerpar", no_argument, &options->innerpar, 1},
        {"dynschedule", no_argument, &options->dynschedule, 1},
        {"distmem", no_argument, &options->distmem, 1},
#ifdef PLUTO_OPENCL
        {"opencl", no_argument, &options->opencl, 1},
#endif
        {"commopt", no_argument, &options->commopt, 1},
        {"commopt_dep_split", no_argument, &options->commopt_dep_split, 1},
        {"dsfo_pack_foifi", no_argument, &options->dsfo_pack_foifi, 1},
        {"commopt_foifi", no_argument, &options->commopt_foifi, 1},
        {"nocommopt", no_argument, &options->commopt, 0},
        {"commreport", no_argument, &options->commreport, 1},
        {"variables_not_global", no_argument, &options->variables_not_global, 1},
        {"mpiomp", no_argument, &options->mpiomp, 1},
        {"blockcyclic", no_argument, &options->blockcyclic, 1},
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
        {"cloogsh", no_argument, &options->cloogsh, 1},
        {"nocloogbacktrack", no_argument, &options->cloogbacktrack, 0},
        {"cyclesize", required_argument, 0, 'S'},
        {"forceparallel", required_argument, 0, 'p'},
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
        {"isldep", no_argument, &options->isldep, 1},
        {"isldepcompact", no_argument, &options->isldepcompact, 1},
        {"readscoplib", no_argument, &options->readscoplib, 1},
        {"islsolve", no_argument, &options->islsolve, 1},
        {"fusesends", no_argument, &options->fusesends, 1},
        {0, 0, 0, 0}
    };


    /* Read command-line options */
    while (1) {
        option = getopt_long(argc, argv, "bhiqvf:l:F:L:c:o:", pluto_options,
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
            case 'S':
                options->cyclesize = atoi(optarg);
                break;
            case 'b':
                options->bee = 1;
                break;
            case 'c':
                options->context = atoi(optarg);
                break;
            case 'f':
                options->ft = atoi(optarg);
                break;
            case 'g':
                break;
            case 'h':
                usage_message();
                return 2;
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
            case 'o':
                options->out_file = strdup(optarg);
                break;
            case 'p':
                options->forceparallel = atoi(optarg);
                break;
            case 'q':
                options->silent = 1;
                break;
            case 's':
                break;
            case 'u':
                options->ufactor = atoi(optarg);
                break;
            case 'v':
                printf("PLUTO %s - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2008  Uday Kumar Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n", PLUTO_VERSION);
                pluto_options_free(options);
                return 3;
            default:
                usage_message();
                pluto_options_free(options);
                return 4;
        }
    }


    if (optind <= argc-1)   {
        srcFileName = alloca(strlen(argv[optind])+1);
        strcpy(srcFileName, argv[optind]);
    }else{
        /* No non-option argument was specified */
        usage_message();
        pluto_options_free(options);
        return 5;
    }

    src_fp  = fopen(srcFileName, "r");

    if (!src_fp)   {
        fprintf(stderr, "pluto: error opening source file: '%s'\n", srcFileName);
        pluto_options_free(options);
        return 6;
    }

    if (options->fusesends && options->mpiomp) {
        fprintf(stderr, "Error: fusesends should not be used with mpiomp\n");
        return 7;
    }

    if (options->isldepcompact && options->distmem) {
        fprintf(stderr, "Error: shouldn't compact deps for distmem parallelization\n");
        return 7;
    }

    /* Extract polyhedral representation from input program */
    scoplib_scop_p scop;

    clan_options_p clanOptions = clan_options_malloc();

    if (options->readscoplib) scop = scoplib_scop_read(src_fp);
    else scop = clan_scop_extract(src_fp, clanOptions);

    if (!scop || !scop->statement)   {
        fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                srcFileName);
        pluto_options_free(options);
        return 8;
    }
    FILE *srcfp = fopen(".srcfilename", "w");
    if (srcfp)    {
        fprintf(srcfp, "%s\n", srcFileName);
        fclose(srcfp);
    }

    /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */

    /* Convert clan scop to Pluto program */
    PlutoProg *prog = scop_to_pluto_prog(scop, options);

    clan_options_free(clanOptions);

    /* Backup irregular program portion in .scop. */
    char* irroption = scoplib_scop_tag_content(scop, "<irregular>",
            "</irregular>");

    IF_DEBUG2(pluto_deps_print(stdout, prog));
    IF_DEBUG2(pluto_stmts_print(stdout, prog->stmts, prog->nstmts));


    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i]->dim;
    }

    /* Make options consistent */
    if (options->distmem == 1 && options->parallel == 0)    {
        options->parallel = 1;
    }

    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

    if (options->dynschedule && !options->tile)  {
        fprintf(stderr, "[Pluto] WARNING: dynschedule needs tile to be on; turning on tile\n");
        options->tile = 1;
    }

    if (options->dynschedule && (options->parallel || options->multipipe))  {
        fprintf(stderr, "[Pluto] WARNING: --parallel or --multipipe options not needed with --dynschedule; turning off parallel and multipipe\n");
        options->parallel = 0;
        options->multipipe = 0;
    }

    if (options->dynschedule && options->distmem)  {
        fprintf(stderr, "[Pluto] WARNING: dynschedule is not compatible with distributed memory yet; turning off distmem\n");
        options->distmem = 0;
    }

    //reset commopt and commopt_foifi when commopt_dep_split is selected, default of commopt is 1
    if(options->commopt && options->commopt_dep_split) {
        options->commopt = 0;
    }
    if(options->commopt_foifi && options->commopt_dep_split) {
        options->commopt_foifi = 0;
    }

    //reset commopt when commopt_foifi is selected, default of commopt is 1
    if(options->commopt && options->commopt_foifi) {
        options->commopt = 0;
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
        fprintf(stdout, "[Pluto] Maximum domain dimensionality: %d\n", prog->nvar);
        fprintf(stdout, "[Pluto] Number of parameters: %d\n", prog->npar);
    }

    /* Auto transformation */
    if (!options->identity) {
        pluto_auto_transform(prog);
    }
    pluto_detect_transformation_properties(prog);

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Affine transformations [<iter coeff's> <const>]\n\n");
        /* Print out transformations */
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);
    }

    if (options->tile)   {
        pluto_tile(prog);
    }else{
        if (options->intratileopt) {
            int retval = pluto_intra_tile_optimize(prog, 0); 
            if (retval) {
                /* Detect properties again */
                pluto_detect_transformation_properties(prog);
                if (!options->silent) {
                    printf("[Pluto] after intra tile opt\n");
                    pluto_transformations_pretty_print(prog);
                }
            }
        }
    }

    if (options->parallel && !options->tile && !options->identity)   {
        /* Obtain wavefront/pipelined parallelization by skewing if
         * necessary */
        int nbands;
        Band **bands;
        bands = pluto_get_outermost_permutable_bands(prog, &nbands);
        bool retval = create_tile_schedule(prog, bands, nbands);
        pluto_bands_free(bands, nbands);

        /* If the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user */
        if (retval)   {
            printf("[Pluto] WARNING: pipelined parallelism exists and --tile is not used.\n");
            printf("use --tile for better parallelization \n");
            IF_DEBUG(fprintf(stdout, "[Pluto] After skewing:\n"););
            IF_DEBUG(pluto_transformations_pretty_print(prog););
            IF_DEBUG(pluto_print_hyperplane_properties(prog););
        }
    }

    if (options->tile && !options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);
    }

    if (options->unroll || options->polyunroll)    {
        /* Will generate a .unroll file */
        /* plann/plorc needs a .params */
        FILE *paramsFP = fopen(".params", "w");
        if (paramsFP)   {
            int i;
            for (i=0; i<prog->npar; i++)  {
                fprintf(paramsFP, "%s\n", prog->params[i]);
            }
            fclose(paramsFP);
        }
        detect_mark_unrollable_loops(prog);
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }

    int distretval = 1;

#ifdef PLUTO_OPENCL
    if (options->distmem || options->opencl)  {
#else
        if (options->distmem)  {
#endif
            distretval = pluto_distmem_parallelize(prog);
            pluto_compute_dep_satisfaction_complex(prog);
            IF_DEBUG(pluto_transformations_pretty_print(prog));
            IF_DEBUG(pluto_print_hyperplane_properties(prog));
        }


        /* NO MORE TRANSFORMATIONS BEYOND THIS POINT */
        /* Since meta info about loops
         * is printed to be processed by scripts - if transformations are
         * performed, changed loop order/iterator names will be missed  */
        gen_unroll_file(prog);

        char *outFileName;
        char *cloogFileName;
        char *dynschedFileName;
        char *bname, *basec;
        if (options->out_file == NULL)  {
            /* Get basename, remove .c extension and append a new one */
            basec = strdup(srcFileName);
            bname = basename(basec);

            if (strlen(bname) >= 2 && !strcmp(bname+strlen(bname)-2, ".c")) {
                outFileName = malloc(strlen(bname)-2+strlen(".pluto.c")+1);
                strncpy(outFileName, bname, strlen(bname)-2);
                outFileName[strlen(bname)-2] = '\0';
            }else{
                outFileName = malloc(strlen(bname)+strlen(".pluto.c")+1);
                strcpy(outFileName, bname);
            }
            strcat(outFileName, ".pluto.c");
        }else{
            basec = strdup(options->out_file);
            bname = basename(basec);

            outFileName = malloc(strlen(options->out_file)+1);
            strcpy(outFileName, options->out_file);
        }

        if (strlen(bname) >= 2 && !strcmp(bname+strlen(bname)-2, ".c")) {
            cloogFileName = malloc(strlen(bname)-2+strlen(".pluto.cloog")+1);
            dynschedFileName = malloc(strlen(bname)-2+strlen(".pluto.append.c")+1);
            strncpy(cloogFileName, bname, strlen(bname)-2);
            strncpy(dynschedFileName, bname, strlen(bname)-2);
            cloogFileName[strlen(bname)-2] = '\0';
            dynschedFileName[strlen(bname)-2] = '\0';
        }else{
            cloogFileName = malloc(strlen(bname)+strlen(".pluto.cloog")+1);
            dynschedFileName = malloc(strlen(bname)+strlen(".pluto.append.c")+1);
            strcpy(cloogFileName, bname);
            strcpy(dynschedFileName, bname);
        }
        strcat(cloogFileName, ".pluto.cloog");
        strcat(dynschedFileName, ".pluto.append.c");
        free(basec);

        cloogfp = fopen(cloogFileName, "w+");
        if (!cloogfp)   {
            fprintf(stderr, "[Pluto] Can't open .cloog file: '%s'\n", cloogFileName);
            free(cloogFileName);
            pluto_options_free(options);
            pluto_prog_free(prog);
            return 9;
        }
        free(cloogFileName);

        outfp = fopen(outFileName, "w");
        if (!outfp) {
            fprintf(stderr, "[Pluto] Can't open file '%s' for writing\n", outFileName);
            free(outFileName);
            pluto_options_free(options);
            pluto_prog_free(prog);
            fclose(cloogfp);
            return 10;
        }

        if (options->dynschedule) {
            dynschedfp = fopen(dynschedFileName, "w");
            if (!dynschedfp) {
                fprintf(stderr, "[Pluto] Can't open file %s for writing\n", dynschedFileName);
                free(dynschedFileName);
                pluto_options_free(options);
                pluto_prog_free(prog);
                fclose(cloogfp);
                fclose(outfp);
                return 11;
            }
        }

        pluto_detect_scalar_dimensions(prog);
        if (options->moredebug) {
            printf("After scalar dimension detection (final transformations)\n");
            pluto_transformations_pretty_print(prog);
        }

        /* Generate .cloog file */
        pluto_gen_cloog_file(cloogfp, prog);
        /* Add the <irregular> tag from clan, if any */
        if (irroption != NULL) {
            fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
            free(irroption);
        }
        rewind(cloogfp);

        /* Very important: Dont change the order of calls to print_dynsched_file
         * between pluto_gen_cloog_file() and pluto_*_codegen()
         */
        if (options->dynschedule) {
            print_dynsched_file(srcFileName, cloogfp, dynschedfp, prog);
            rewind(cloogfp);
        }

        /* Generate code using Cloog and add necessary stuff before/after code */
        if (options->distmem && !distretval)   {
            pluto_distmem_codegen(prog, cloogfp, outfp);
        }else{
            if (options->distmem) { // no parallel loops to distribute
                // ensure only one processor will output the data
                fprintf(outfp, "#include <mpi.h>\n\n");
                fprintf(outfp, "#define MPI \n\n");
                fprintf(outfp, "\n##ifndef GLOBAL_MY_RANK\n\tint my_rank;\n##endif\n");
                fprintf(outfp, "\tMPI_Init(NULL, NULL);\n");
                fprintf(outfp, "\tMPI_Comm_rank(MPI_COMM_WORLD, &my_rank);\n");
                fprintf(outfp, "\tMPI_Finalize();\n");
            }
            pluto_multicore_codegen(cloogfp, options->dynschedule ? dynschedfp : outfp, prog);
        }

        FILE *tmpfp = fopen(".outfilename", "w");
        if (tmpfp)    {
            fprintf(tmpfp, "%s\n", outFileName);
            fclose(tmpfp);
            printf( "[Pluto] Output written to %s\n", outFileName);
        }
        free(outFileName);

        if (options->dynschedule) {
            tmpfp = fopen(".appendfilename", "w");
            if (tmpfp) {
                fprintf(tmpfp, "%s\n", dynschedFileName);
                fclose(tmpfp);
            }
        }
        free(dynschedFileName);

        /* create main file for dynschedule */
        if (options->dynschedule) {
            fprintf(outfp,"\n\n#include \"scheduler.h\"\n");
            fprintf(outfp,"\n\n\tgenerate_dag();\n");
            fprintf(outfp,"\t##ifdef __STATIC_SCHEDULE__\n");
            fprintf(outfp,"\t\tinit_schedule(atoi(argv[1]));\n");
            fprintf(outfp,"\t\tschedule_dag_using_mcp();\n");
            fprintf(outfp,"\t##endif\n");
            fprintf(outfp,"\tdag_execute();\n");
            fprintf(outfp,"\t##ifdef __STATIC_SCHEDULE__\n");
            fprintf(outfp,"\t\tfree_schedule();\n");
            fprintf(outfp,"\t##endif\n");
            fprintf(outfp,"\tfree_dag();\n");
        }

        if (options->dynschedule) {
            fclose(dynschedfp);
        }
        fclose(cloogfp);
        fclose(outfp);

        pluto_options_free(options);

        pluto_prog_free(prog);

        scoplib_scop_free(scop);

        return 0;
    }
