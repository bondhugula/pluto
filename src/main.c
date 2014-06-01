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


#include "osl/scop.h"
#include "osl/generic.h"
#include "osl/extensions/irregular.h"

#include "pluto.h"
#include "transforms.h"
#include "math_support.h"
#include "post_transform.h"
#include "program.h"
#include "version.h"

#include "clan/clan.h"
#include "candl/candl.h"
#include "candl/scop.h"

#include "pet.h"

PlutoOptions *options;

void usage_message(void)
{
    fprintf(stdout, "Usage: polycc <input.c> [options] [-o output]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --pet                     Use libpet for polyhedral extraction instead of clan [default - clan]\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "       --tile               Tile for locality\n");
    fprintf(stdout, "       --intratileopt       Optimize intra-tile execution order for locality\n");
    fprintf(stdout, "       --l2tile                  Tile a second time (typically for L2 cache) - disabled by default \n");
    fprintf(stdout, "       --parallel                Automatically parallelize (generate OpenMP)\n");
    fprintf(stdout, "       | --parallelize\n");
    fprintf(stdout, "       --lbtile                  Enables full dimensional concurrent start\n");
    fprintf(stdout, "       --partlbtile              Enables one-dimensional concurrent start\n");
    fprintf(stdout, "       --[no]prevector           Transform for and mark loops for (icc) vectorization (enabled by default)\n");
    fprintf(stdout, "       --innerpar                Choose pure inner parallelism over pipelined/wavefront parallelism\n");
    fprintf(stdout, "       --multipipe               Extract two or more degrees of pipelined parallelism if possible;\n");
    fprintf(stdout, "                                 by default one degree is extracted (if it exists)\n");
    fprintf(stdout, "       --variables_not_global    Variables not declared globally (if so, macros provide variable declarations)\n");
    fprintf(stdout, "\n   Runtime               Options related to compilation for a runtime\n");
    fprintf(stdout, "       --dynschedule             Dynamically schedule tasks on processors using Synthesized Runtime Interface\n");
    fprintf(stdout, "                                     (for shared and distributed memory)\n");
    fprintf(stdout, "       --dynschedule_graph       Dynamically schedule tasks on processors using Intel TBB Flow Graph\n");
    fprintf(stdout, "                                     (only for shared-memory)\n");
    fprintf(stdout, "       --dataflow                Alias to --dynschedule\n");
    fprintf(stdout, "\n   Architecture          Options related to compilation for a specific architecture\n");
#ifdef PLUTO_OPENCL
    fprintf(stdout, "       --opencl                  Generate OpenCL code \n");
#endif
    fprintf(stdout, "\n   Analysis              Options related to analyzing generated code (for Runtime or Architecture)\n");
    fprintf(stdout, "       --timereport              Generate code to report communication volume (for distributed-memory only)\n");
    fprintf(stdout, "                                     and analysis of time (for distributed-memory or runtime)\n");
    fprintf(stdout, "\n   Communication code    Options related to communication code generation for distributed-memory\n");
    fprintf(stdout, "       --commopt                 Generate communication code using Flow-Out (FO) scheme (enabled by default)\n");
    fprintf(stdout, "       --commopt_foifi           Generate communication code using Flow-Out Intersection Flow-In (FOIFI) scheme\n");
    fprintf(stdout, "       --commopt_fop             Generate communication code using Flow-Out Partitioning (FOP) scheme (multicast pack by default)\n");
    fprintf(stdout, "       --fop_unicast_runtime     FOP: generate code to choose between unicast and multicast pack at runtime\n\n");
    fprintf(stdout, "       --rar                  Consider RAR dependences too (disabled by default)\n");
    fprintf(stdout, "       --[no]unroll           Unroll-jam (disabled by default)\n");
    fprintf(stdout, "       --ufactor=<factor>     Unroll-jam factor (default is 8)\n");
    fprintf(stdout, "       --[no]prevector           Transform for and mark loops for (icc) vectorization (enabled by default)\n");
    fprintf(stdout, "       --context=<context>    Parameters are at least as much as <context>\n");
    fprintf(stdout, "       --forceparallel=<bitvec>  6 bit-vector of depths (1-indexed) to force parallel (0th bit represents depth 1)\n");
    fprintf(stdout, "       --[no]isldep              Use ISL-based dependence tester (disabled by default)\n");
    fprintf(stdout, "       --islsolve             Use ISL as ilp solver\n");
    fprintf(stdout, "       --readscoplib             Read input from a scoplib file\n");
    fprintf(stdout, "       --[no]lastwriter          Work with refined dependences (last conflicting access is computed for RAW/WAW)\n");
    fprintf(stdout, "                                     (disabled by default)\n");
    fprintf(stdout, "       --bee                  Generate pragmas for Bee+Cl@k\n\n");
    fprintf(stdout, "       --indent  | -i         Indent generated code (disabled by default)\n");
    fprintf(stdout, "       --silent  | -q         Silent mode; no output as long as everything goes fine (disabled by default)\n");
    fprintf(stdout, "       --help    | -h         Print this help menu\n");
    fprintf(stdout, "       --version | -v         Display version number\n");
    fprintf(stdout, "\n   Fusion                Options to control fusion heuristic\n");
    fprintf(stdout, "       --nofuse               Do not fuse across SCCs of data dependence graph\n");
    fprintf(stdout, "       --maxfuse              Maximal fusion\n");
    fprintf(stdout, "       --smartfuse [default]  Heuristic (in between nofuse and maxfuse)\n");
    fprintf(stdout, "\n   Code generation       Options to control Cloog code generation\n");
    fprintf(stdout, "       --nocloogbacktrack        Do not call Cloog with backtrack (default - backtrack)\n");
    fprintf(stdout, "       --cloogsh                 Ask Cloog to use simple convex hull (default - off)\n");
    fprintf(stdout, "\n   Debugging\n");
    fprintf(stdout, "       --debug        Verbose output\n");
    fprintf(stdout, "       --moredebug    More verbose output\n");
    fprintf(stdout, "\nTo report bugs, please email <pluto-development@googlegroups.com>\n\n");
}

int main(int argc, char *argv[])
{
    int i;

    FILE *src_fp;

    struct pet_scop *pscop;
    isl_ctx *pctx = isl_ctx_alloc();

    int option;
    int option_index = 0;

    char *srcFileName;

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
        {"intratileopt", no_argument, &options->intratileopt, 1},
        {"nointratileopt", no_argument, &options->intratileopt, 0},
        {"lbtile", no_argument, &options->lbtile, 1},
        {"pet", no_argument, &options->pet, 1},
        {"partlbtile", no_argument, &options->partlbtile, 1},
        {"dynschedule", no_argument, &options->dynschedule, 1},
        {"dataflow", no_argument, &options->dynschedule, 1},
        {"dynschedule_graph", no_argument, &options->dynschedule_graph, 1},
        {"dynschedule_graph_old", no_argument, &options->dynschedule_graph_old, 1},
        {"dyn_trans_deps_tasks", no_argument, &options->dyn_trans_deps_tasks, 1},
        {"debug", no_argument, &options->debug, 1},
        {"moredebug", no_argument, &options->moredebug, 1},
        {"rar", no_argument, &options->rar, 1},
        {"identity", no_argument, &options->identity, 1},
        {"nofuse", no_argument, &options->fuse, NO_FUSE},
        {"maxfuse", no_argument, &options->fuse, MAXIMAL_FUSE},
        {"smartfuse", no_argument, &options->fuse, SMART_FUSE},
        {"parallel", no_argument, &options->parallel, 1},
        {"parallelize", no_argument, &options->parallel, 1},
        {"innerpar", no_argument, &options->innerpar, 1},
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
        {"nolastwriter", no_argument, &options->nolastwriter, 1},
        {"nobound", no_argument, &options->nobound, 1},
        {"scalpriv", no_argument, &options->scalpriv, 1},
        {"isldep", no_argument, &options->isldep, 1},
        {"noisldep", no_argument, &options->noisldep, 1},
        {"isldepcompact", no_argument, &options->isldepcompact, 1},
        {"readscop", no_argument, &options->readscop, 1},
        {"islsolve", no_argument, &options->islsolve, 1},
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

    PlutoProg *prog; 
    char *irroption = NULL;
    osl_scop_p scop = NULL;

    /* Extract polyhedral representation from input program */
    if (options->pet) {
        pscop = pet_scop_extract_from_C_source(pctx, srcFileName, NULL);

        if (!pscop) {
            fprintf(stdout, "[pluto] No SCoPs extracted or error extracting SCoPs  using pet\n");
            pluto_options_free(options);
            return 12;
        }
        prog = pet_to_pluto_prog(pscop, options);

        FILE *srcfp = fopen(".srcfilename", "w");
        if (srcfp)    {
            fprintf(srcfp, "%s\n", srcFileName);
            fclose(srcfp);
        }
    }else{
    if(!strcmp(srcFileName, "stdin")){  //read from stdin
        src_fp = stdin;
        osl_interface_p registry = osl_interface_get_default_registry();
        scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
        }else{  // read from regular file

      src_fp  = fopen(srcFileName, "r");
  
      if (!src_fp)   {
                fprintf(stderr, "pluto: error opening source file: '%s'\n", 
                        srcFileName);
          pluto_options_free(options);
          return 6;
      }
  
      clan_options_p clanOptions = clan_options_malloc();
  
      if (options->readscop){
        osl_interface_p registry = osl_interface_get_default_registry();
        scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
      }else{
        scop = clan_scop_extract(src_fp, clanOptions);
      }
  
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
  
      clan_options_free(clanOptions);

      /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */
    }

    /* Convert clan scop to Pluto program */
        prog = scop_to_pluto_prog(scop, options);

    /* Backup irregular program portion in .scop. */
    osl_irregular_p irreg_ext = NULL;
    irreg_ext = osl_generic_lookup(scop->extension, OSL_URI_IRREGULAR);
    if(irreg_ext!=NULL)
      irroption = osl_irregular_sprint(irreg_ext);  //TODO: test it
    osl_irregular_free(irreg_ext);
    }
    IF_MORE_DEBUG(pluto_prog_print(stdout, prog));

    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i]->dim;
    }

    /* Make options consistent */
    if (options->noisldep == 1) {
        options->isldep = 0;
    }

    if (options->nolastwriter == 1) {
        options->lastwriter = 0;
    }

    if (options->identity == 1) {
        options->partlbtile = 0;
        options->lbtile = 0;
    }

    if (options->dynschedule && options->dynschedule_graph) {
        options->dynschedule_graph = 0;
    }

    if (options->dynschedule && options->dynschedule_graph_old) {
        options->dynschedule_graph_old = 0;
    }

    if (options->dynschedule_graph && options->dynschedule_graph_old) {
        options->dynschedule_graph_old = 0;
    }

    if (options->partlbtile == 1 && options->lbtile == 0)    {
        options->lbtile = 1;
    }

    if (options->lbtile == 1 && options->tile == 0)    {
        options->tile = 1;
    }

    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

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

        if (options->lbtile) {
            pluto_reschedule_tile(prog);
        }
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
        bool retval = pluto_create_tile_schedule(prog, bands, nbands);
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
        pluto_detect_mark_unrollable_loops(prog);
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }

    if(!options->pet && !strcmp(srcFileName, "stdin")){  
        //input stdin == output stdout
      pluto_populate_scop(scop, prog, options);
      osl_scop_print(stdout, scop);
    }else{  // do the usual Pluto stuff
  
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
            strncpy(cloogFileName, bname, strlen(bname)-2);
            cloogFileName[strlen(bname)-2] = '\0';
        }else{
            cloogFileName = malloc(strlen(bname)+strlen(".pluto.cloog")+1);
            strcpy(cloogFileName, bname);
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
  
        if (!options->pet) pluto_detect_scalar_dimensions(prog);
        if (options->moredebug) {
            printf("After scalar dimension detection (final transformations)\n");
            pluto_transformations_pretty_print(prog);
        }

      /* Generate .cloog file */
      pluto_gen_cloog_file(cloogfp, prog);
      /* Add the <irregular> tag from clan, if any */
        if(!options->pet) {
            if (irroption) {
          fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
          free(irroption);
      }
        }
        rewind(cloogfp);
  
        /* Very important: Dont change the order of calls to print_dynsched_file
         * between pluto_gen_cloog_file() and pluto_*_codegen()
         */
  
      /* Generate code using Cloog and add necessary stuff before/after code */
        pluto_multicore_codegen(cloogfp, outfp, prog);
  
      FILE *tmpfp = fopen(".outfilename", "w");
      if (tmpfp)    {
          fprintf(tmpfp, "%s\n", outFileName);
          fclose(tmpfp);
          printf( "[Pluto] Output written to %s\n", outFileName);
      }
        free(outFileName);
 
      fclose(cloogfp);
      fclose(outfp);
    }

    if (!options->pet) osl_scop_free(scop);

    pluto_options_free(options);

    pluto_prog_free(prog);

    return 0;
}
