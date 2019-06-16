/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2015 Uday Bondhugula
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
 * Top-level file for 'pluto' executable.
 */
#include <assert.h>
#include <getopt.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "osl/extensions/irregular.h"
#include "osl/generic.h"
#include "osl/scop.h"

#include "math_support.h"
#include "pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"
#include "version.h"

#include "candl/candl.h"
#include "candl/scop.h"
#include "clan/clan.h"

#include "pet.h"

PlutoOptions *options;

void usage_message(void) {
  fprintf(stdout, "Usage: polycc <input.c> [options] [-o output]\n");
  fprintf(stdout, "\nOptions:\n");
  fprintf(stdout, "       --pet                     Use libpet for polyhedral "
                  "extraction instead of clan [default - clan]\n");
  fprintf(stdout, "       --isldep                  Use ISL-based dependence "
                  "tester (enabled by default)\n");
  fprintf(
      stdout,
      "       --candldep                Use Candl as the dependence tester\n");
  fprintf(stdout, "       --[no]lastwriter          Remove transitive "
                  "dependences (last conflicting access is computed for "
                  "RAW/WAW)\n");
  fprintf(stdout, "                                 (disabled by default)\n");
  fprintf(stdout,
          "       --islsolve [default]      Use ISL as ILP solver (default)\n");
  fprintf(stdout, "       --pipsolve                Use PIP as ILP solver\n");
#ifdef GLPK
  fprintf(stdout, "       --glpk                    Use GLPK as ILP solver "
                  "(default in case of pluto-lp and pluto-dfp)\n");
#endif
#if defined GLPK || defined GUROBI
  fprintf(stdout,
          "       --lp                      Solve MIP instead of ILP\n");
  fprintf(stdout, "       --dfp                     Use Pluto-lp-dfp instead "
                  "of pluto-ilp [disabled by default]\n");
  fprintf(stdout, "       --ilp                     Use ILP in pluto-lp-dfp "
                  "instead of LP\n");
  fprintf(stdout, "       --lpcolor                 Color FCG based on the "
                  "solutions of the lp-problem [disabled by default]\n");
  fprintf(stdout, "       --clusterscc              Cluster the statemtns of "
                  "an SCC. This is supported only availabe with decoupled "
                  "approach [disabled by default]\n");
#endif
  fprintf(stdout, "\n");
#ifdef GUROBI
  fprintf(stdout,
          "       --gurobi                  Use Gurobi as ILP solver\n");
#endif
  fprintf(stdout, "\n");
  fprintf(stdout,
          "\n  Optimizations          Options related to optimization\n");
  fprintf(stdout, "       --tile                    Tile for locality "
                  "[disabled by default]\n");
  fprintf(stdout, "       --[no]intratileopt        Optimize intra-tile "
                  "execution order for locality [enabled by default]\n");
  fprintf(stdout, "       --l2tile                  Tile a second time "
                  "(typically for L2 cache) [disabled by default] \n");
  fprintf(stdout, "       --parallel                Automatically parallelize "
                  "(generate OpenMP pragmas) [disabled by default]\n");
  fprintf(stdout, "    or --parallelize\n");
  fprintf(stdout, "       --[no]diamond-tile        Performs diamond tiling "
                  "(enabled by default)\n");
  fprintf(stdout, "       --full-diamond-tile       Enables full-dimensional "
                  "concurrent start\n");
  fprintf(stdout, "       --[no]prevector           Mark loops for (icc/gcc) "
                  "vectorization (enabled by default)\n");
  fprintf(stdout, "       --multipar                Extract all degrees of "
                  "parallelism [disabled by default];\n");
  fprintf(stdout, "                                    by default one degree "
                  "is extracted within any schedule sub-tree (if it exists)\n");
  fprintf(stdout, "       --innerpar                Choose pure inner "
                  "parallelism over pipelined/wavefront parallelism [disabled "
                  "by default]\n");
  fprintf(stdout,
          "\n   Fusion                Options to control fusion heuristic\n");
  fprintf(stdout, "       --nofuse                  Do not fuse across SCCs of "
                  "data dependence graph\n");
  fprintf(stdout, "       --maxfuse                 Maximal fusion\n");
  fprintf(stdout, "       --smartfuse [default]     Heuristic (in between "
                  "nofuse and maxfuse)\n");
  fprintf(stdout, "       --typedfuse               Typed fusion. Fuses SCCs "
                  "only when there is no loss of parallelism\n");
  fprintf(stdout, "       --hybridfuse              Typed fusion at outer "
                  "levels and max fuse at inner level\n");
  fprintf(stdout, "       --delayedcut              Delays the cut between "
                  "SCCs of different dimensionalities in dfp approach\n");
  fprintf(stdout, "\n   Index Set Splitting        \n");
  fprintf(stdout, "       --iss                  \n");
  fprintf(
      stdout,
      "\n   Code generation       Options to control Cloog code generation\n");
  fprintf(stdout, "       --nocloogbacktrack        Do not call Cloog with "
                  "backtrack (default - backtrack)\n");
  fprintf(stdout, "       --cloogsh                 Ask Cloog to use simple "
                  "convex hull (default - off)\n");
  fprintf(stdout, "       --codegen-context=<value> Parameters are at least as "
                  "much as <value>\n");
  fprintf(stdout, "\n   Miscellaneous\n");
  fprintf(stdout, "       --rar                     Consider RAR dependences "
                  "too (disabled by default)\n");
  fprintf(
      stdout,
      "       --[no]unroll              Unroll-jam (disabled by default)\n");
  fprintf(
      stdout,
      "       --ufactor=<factor>        Unroll-jam factor (default is 8)\n");
  fprintf(stdout, "       --forceparallel=<bitvec>  6 bit-vector of depths "
                  "(1-indexed) to force parallel (0th bit represents depth "
                  "1)\n");
  fprintf(stdout,
          "       --readscop                Read input from a scoplib file\n");
  fprintf(stdout,
          "       --bee                     Generate pragmas for Bee+Cl@k\n\n");
  fprintf(stdout, "       --indent  | -i            Indent generated code "
                  "(disabled by default)\n");
  fprintf(stdout, "       --silent  | -q            Silent mode; no output as "
                  "long as everything goes fine (disabled by default)\n");
  fprintf(stdout, "       --help    | -h            Print this help menu\n");
  fprintf(stdout, "       --version | -v            Display version number\n");
  fprintf(stdout, "\n   Debugging\n");
  fprintf(stdout, "       --debug                   Verbose/debug output\n");
  fprintf(stdout,
          "       --moredebug               More verbose/debug output\n");
  fprintf(stdout, "\nTo report bugs, please email "
                  "<pluto-development@googlegroups.com>\n\n");
}

static double rtclock() {
  struct timeval Tp;
  int stat = gettimeofday(&Tp, NULL);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

int main(int argc, char *argv[]) {
  if (argc <= 1) {
    usage_message();
    return 1;
  }

  double t_start_all = rtclock();

  options = pluto_options_alloc();

  int option_index = 0;
  int nolastwriter = 0;

  const struct option pluto_options[] = {
    {"fast-lin-ind-check", no_argument, &options->flic, 1},
    {"flic", no_argument, &options->flic, 1},
    {"tile", no_argument, &options->tile, 1},
    {"notile", no_argument, &options->tile, 0},
    {"noparallel", no_argument, &options->parallel, 0},
    {"intratileopt", no_argument, &options->intratileopt, 1},
    {"nointratileopt", no_argument, &options->intratileopt, 0},
    {"pet", no_argument, &options->pet, 1},
    {"diamond-tile", no_argument, &options->diamondtile, 1},
    {"nodiamond-tile", no_argument, &options->diamondtile, 0},
    {"full-diamond-tile", no_argument, &options->fulldiamondtile, 1},
    {"debug", no_argument, &options->debug, true},
    {"moredebug", no_argument, &options->moredebug, true},
    {"rar", no_argument, &options->rar, 1},
    {"identity", no_argument, &options->identity, 1},
    {"nofuse", no_argument, (int *)&options->fuse, kNoFuse},
    {"maxfuse", no_argument, (int *)&options->fuse, kMaximalFuse},
    {"smartfuse", no_argument, (int *)&options->fuse, kSmartFuse},
    {"typedfuse", no_argument, (int *)&options->fuse, kTypedFuse},
    {"hybridfuse", no_argument, &options->hybridcut, 1},
    {"delayedcut", no_argument, &options->delayed_cut, 1},
    {"parallel", no_argument, &options->parallel, 1},
    {"parallelize", no_argument, &options->parallel, 1},
    {"innerpar", no_argument, &options->innerpar, 1},
    {"iss", no_argument, &options->iss, 1},
    {"unroll", no_argument, &options->unroll, 1},
    {"nounroll", no_argument, &options->unroll, 0},
    {"bee", no_argument, &options->bee, 1},
    {"ufactor", required_argument, 0, 'u'},
    {"prevector", no_argument, &options->prevector, 1},
    {"noprevector", no_argument, &options->prevector, 0},
    {"codegen-context", required_argument, 0, 'c'},
    {"coeff-bound", required_argument, 0, 'C'},
    {"cloogf", required_argument, 0, 'F'},
    {"cloogl", required_argument, 0, 'L'},
    {"cloogsh", no_argument, &options->cloogsh, 1},
    {"nocloogbacktrack", no_argument, &options->cloogbacktrack, 0},
    {"cyclesize", required_argument, 0, 'S'},
    {"forceparallel", required_argument, 0, 'p'},
    {"ft", required_argument, 0, 'f'},
    {"lt", required_argument, 0, 'l'},
    {"multipar", no_argument, &options->multipar, 1},
    {"l2tile", no_argument, &options->l2tile, 1},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {"indent", no_argument, 0, 'i'},
    {"silent", no_argument, &options->silent, 1},
    {"lastwriter", no_argument, &options->lastwriter, 1},
    {"nolastwriter", no_argument, &nolastwriter, 1},
    {"nodepbound", no_argument, &options->nodepbound, 1},
    {"scalpriv", no_argument, &options->scalpriv, 1},
    {"isldep", no_argument, &options->isldep, 1},
    {"candldep", no_argument, &options->candldep, 1},
    {"isldepaccesswise", no_argument, &options->isldepaccesswise, 1},
    {"isldepstmtwise", no_argument, &options->isldepaccesswise, 0},
    {"isldepcoalesce", no_argument, &options->isldepcoalesce, 1},
    {"readscop", no_argument, &options->readscop, 1},
    {"pipsolve", no_argument, &options->pipsolve, 1},
#ifdef GLPK
    {"glpk", no_argument, &options->glpk, 1},
#endif
#ifdef GUROBI
    {"gurobi", no_argument, &options->gurobi, 1},
#endif
#if defined GLPK || defined GUROBI
    {"lp", no_argument, &options->lp, 1},
    {"dfp", no_argument, &options->dfp, 1},
    {"ilp", no_argument, &options->ilp, 1},
    {"lpcolor", no_argument, &options->lpcolour, 1},
    {"clusterscc", no_argument, &options->scc_cluster, 1},
#endif
    {"islsolve", no_argument, &options->islsolve, 1},
    {"time", no_argument, &options->time, 1},
    {0, 0, 0, 0}
  };

  /* Read command-line options */
  while (1) {
    int option = getopt_long(argc, argv, "bhiqvf:l:F:L:c:o:", pluto_options,
                             &option_index);

    if (option == -1) {
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
      options->codegen_context = atoi(optarg);
      break;
    case 'C':
      options->coeff_bound = atoi(optarg);
      if (options->coeff_bound <= 0) {
        printf("ERROR: coeff-bound should be at least 1\n");
        return 2;
      }
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
      printf(
          "PLUTO version %s - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2015  Uday Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n",
          PLUTO_VERSION);
      pluto_options_free(options);
      return 3;
    default:
      usage_message();
      pluto_options_free(options);
      return 4;
    }
  }

  char *srcFileName;
  if (optind <= argc - 1) {
    srcFileName = (char *)alloca(strlen(argv[optind]) + 1);
    strcpy(srcFileName, argv[optind]);
  } else {
    /* No non-option argument was specified */
    usage_message();
    pluto_options_free(options);
    return 5;
  }

  /* Make options consistent. */
  if (options->isldep && options->candldep) {
    printf("[pluto] ERROR: only one of isldep and candldep should be "
           "specified)\n");
    pluto_options_free(options);
    usage_message();
    return 1;
  }

  /* isldep is the default */
  if (!options->isldep && !options->candldep) {
    options->isldep = 1;
  }

  if (options->lastwriter && options->candldep) {
    printf("[pluto] ERROR: --lastwriter is only supported with --isldep\n");
    pluto_options_free(options);
    usage_message();
    return 1;
  }

  if (options->lastwriter && nolastwriter) {
    printf("[pluto] WARNING: both --lastwriter, --nolastwriter are on\n");
    printf("[pluto] disabling --lastwriter\n");
    options->lastwriter = 0;
  }

  if (options->diamondtile == 1 && options->tile == 0) {
    options->diamondtile = 0;
  }
  if (options->fulldiamondtile == 1 && options->tile == 0) {
    options->diamondtile = 0;
    options->fulldiamondtile = 0;
  }

  if (options->multipar == 1 && options->parallel == 0) {
    fprintf(stdout,
            "Warning: multipar needs parallel to be on; turning on parallel\n");
    options->parallel = 1;
  }

  if (options->multipar == 1 && options->parallel == 0) {
    fprintf(stdout,
            "Warning: multipar needs parallel to be on; turning on parallel\n");
    options->parallel = 1;
  }

  if (options->gurobi) {
    options->islsolve = 0;
  }
#ifdef GLPK
  if (options->lp && !(options->glpk || options->gurobi)) {
    if (!options->silent) {
      printf("[pluto] LP option available with a LP solver only. Using GLPK "
             "for lp solving\n");
    }
    options->glpk = 1;
  }

  /* By default Pluto-dfp uses lp. */
  if (options->dfp && !options->ilp) {
    options->lp = 1;
  }

  if (options->dfp && !(options->glpk || options->gurobi)) {
    if (!options->silent) {
      printf("[pluto] Dfp framework is currently supported with GLPK and "
             "Gurobi solvers.\n");
      printf("[pluto] Using GLPK for constraint solving [default]. Use "
             "--gurobi to use Gurobi instead of GLPK.\n");
    }
    options->glpk = 1;
  }

  if (options->glpk) {
    /* Turn off islsolve */
    options->islsolve = 0;
    options->pipsolve = 0;
  }
#endif

  // If --pipsolve is provided, disable islsolve.
  if (options->pipsolve)
    options->islsolve = 0;

  if (options->dfp && !(options->glpk || options->gurobi)) {
    printf("[pluto] ERROR: DFP framework is currently supported with GLPK or "
           "GUROBI solvers only. Run ./configure --help to for more "
           "information on using different solvers with Pluto.\n");
    pluto_options_free(options);
    usage_message();
    return 1;
  }
  if (options->scc_cluster && !options->dfp) {
    printf("[pluto] Warning: SCC clustering heuristics available with dfp "
           "option (FCG based approach) only. Disabling clustering \n");
  }

  if (options->fuse == kTypedFuse && !options->dfp) {
    printf("[Pluto] WARNING: Typed fusion available with DFP framework only. "
           "Turning off typed fusion.\n");
    options->fuse = kSmartFuse;
  }

  /* Make lastwriter default with dfp. This removes transitive dependences and
   * hence reduces FCG construction time */
  if (options->dfp && !options->lastwriter) {
    if (!options->silent) {
      printf("[pluto] Enabling lastwriter dependence analysis with DFP\n");
    }
    options->lastwriter = 1;
  }
  /* Typed fuse is available with clustered FCG approach only */
  if (options->fuse == kTypedFuse && options->dfp && !options->scc_cluster) {
    if (!options->silent) {
      printf("[pluto] Typed fuse supported only with clustered FCG approach. "
             "Turning on SCC clustering\n");
    }
    options->scc_cluster = 1;
  }

  /* Extract polyhedral representation from osl scop */
  PlutoProg *prog = NULL;

  osl_scop_p scop = NULL;
  char *irroption = NULL;

  double t_d;

  /* Extract polyhedral representation from input program */
  if (options->pet) {
    // Extract using PET.
    isl_ctx *pctx = isl_ctx_alloc_with_pet_options();
    struct pet_scop *pscop =
        pet_scop_extract_from_C_source(pctx, srcFileName, NULL);

    if (!pscop) {
      fprintf(
          stdout,
          "[pluto] No SCoPs extracted or error extracting SCoPs  using pet\n");
      pluto_options_free(options);
      isl_ctx_free(pctx);
      return 12;
    }
    double t_start = rtclock();
    prog = pet_to_pluto_prog(pscop, pctx, options);
    t_d = rtclock() - t_start;

    pet_scop_free(pscop);
    isl_ctx_free(pctx);

    FILE *srcfp = fopen(".srcfilename", "w");
    if (srcfp) {
      fprintf(srcfp, "%s\n", srcFileName);
      fclose(srcfp);
    }
  } else {
    // Extract polyhedral representation using Clan.
    FILE *src_fp;
    if (!strcmp(srcFileName, "stdin")) {
      // Read from stdin.
      src_fp = stdin;
      osl_interface_p registry = osl_interface_get_default_registry();
      double t_start = rtclock();
      scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
      t_d = rtclock() - t_start;
    } else {
      // Read from regular file.
      src_fp = fopen(srcFileName, "r");
      if (!src_fp) {
        fprintf(stderr, "pluto: error opening source file: '%s'\n",
                srcFileName);
        pluto_options_free(options);
        return 6;
      }

      clan_options_p clanOptions = clan_options_malloc();

      if (options->readscop) {
        osl_interface_p registry = osl_interface_get_default_registry();
        double t_start = rtclock();
        scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
        t_d = rtclock() - t_start;
      } else {
        double t_start = rtclock();
        scop = clan_scop_extract(src_fp, clanOptions);
        t_d = rtclock() - t_start;
      }
      fclose(src_fp);

      if (!scop || !scop->statement) {
        fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                srcFileName);
        osl_scop_free(scop);
        pluto_options_free(options);
        return 8;
      }
      FILE *srcfp = fopen(".srcfilename", "w");
      if (srcfp) {
        fprintf(srcfp, "%s\n", srcFileName);
        fclose(srcfp);
      }

      clan_options_free(clanOptions);
    }

    /* Convert clan scop to Pluto program */
    prog = scop_to_pluto_prog(scop, options);

    /* Backup irregular program portion in .scop. */
    osl_irregular_p irreg_ext = NULL;
    irreg_ext =
        (osl_irregular_p)osl_generic_lookup(scop->extension, OSL_URI_IRREGULAR);
    if (irreg_ext != NULL)
      // TODO: test it
      irroption = osl_irregular_sprint(irreg_ext);
    osl_irregular_free(irreg_ext);
  }
  IF_MORE_DEBUG(pluto_prog_print(stdout, prog));

  int dim_sum = 0;
  for (unsigned i = 0; i < prog->nstmts; i++) {
    dim_sum += prog->stmts[i]->dim;
  }

  if (!options->silent) {
    fprintf(stdout, "[pluto] Number of statements: %d\n", prog->nstmts);
    fprintf(stdout, "[pluto] Total number of loops: %d\n", dim_sum);
    fprintf(stdout, "[pluto] Number of deps: %d\n", prog->ndeps);
    fprintf(stdout, "[pluto] Maximum domain dimensionality: %d\n", prog->nvar);
    fprintf(stdout, "[pluto] Number of parameters: %d\n", prog->npar);
  }

  if (options->iss) {
    pluto_iss_dep(prog);
  }

  double t_start = rtclock();
  /* Auto transformation */
  if (!options->identity) {
    pluto_auto_transform(prog);
  }
  double t_t = rtclock() - t_start;

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(
        stdout,
        "[pluto] Affine transformations [<iter coeff's> <param> <const>]\n\n");
    /* Print out transformations */
    pluto_transformations_pretty_print(prog);
  }

  if (options->tile) {
    pluto_tile(prog);
  } else {
    if (options->intratileopt) {
      pluto_intra_tile_optimize(prog, 0);
    }
  }

  if (options->parallel && !options->tile && !options->identity) {
    /* Obtain wavefront/pipelined parallelization by skewing if
     * necessary */
    unsigned nbands;
    Band **bands;
    pluto_compute_dep_satisfaction(prog);
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    bool retval = pluto_create_tile_schedule(prog, bands, nbands);
    pluto_bands_free(bands, nbands);

    /* If the user hasn't supplied --tile and there is only pipelined
     * parallelism, we will warn the user */
    if (retval) {
      printf("[pluto] WARNING: pipelined parallelism exists and --tile is not "
             "used.\n");
      printf("\tUse --tile for better parallelization \n");
      fprintf(stdout, "[pluto] After skewing:\n");
      pluto_transformations_pretty_print(prog);
    }
  }

  if (options->unroll) {
    /* Will generate a .unroll file */
    /* plann/plorc needs a .params */
    FILE *paramsFP = fopen(".params", "w");
    if (paramsFP) {
      int i;
      for (i = 0; i < prog->npar; i++) {
        fprintf(paramsFP, "%s\n", prog->params[i]);
      }
      fclose(paramsFP);
    }
    pluto_detect_mark_unrollable_loops(prog);
  }

  double t_c = 0.0;

  if (!options->pet && !strcmp(srcFileName, "stdin")) {
    // input stdin == output stdout
    pluto_populate_scop(scop, prog, options);
    osl_scop_print(stdout, scop);
  } else {
    // Do the usual Pluto stuff.

    /* NO MORE TRANSFORMATIONS BEYOND THIS POINT */
    /* Since meta info about loops is printed to be processed by scripts - if
     * transformations are performed, changed loop order/iterator names will
     * be missed. */
    gen_unroll_file(prog);

    char *basec, *bname;
    char *outFileName;
    if (options->out_file == NULL) {
      /* Get basename, remove .c extension and append a new one */
      basec = strdup(srcFileName);
      bname = basename(basec);

      /* max size when tiled.* */
      outFileName = (char *)malloc(strlen(bname) + strlen(".pluto.c") + 1);

      if (strlen(bname) >= 2 && !strcmp(bname + strlen(bname) - 2, ".c")) {
        memcpy(outFileName, bname, strlen(bname) - 2);
        outFileName[strlen(bname) - 2] = '\0';
      } else {
        outFileName = (char *)malloc(strlen(bname) + strlen(".pluto.c") + 1);
        strcpy(outFileName, bname);
      }
      strcat(outFileName, ".pluto.c");
    } else {
      basec = strdup(options->out_file);
      bname = basename(basec);

      outFileName = (char *)malloc(strlen(options->out_file) + 1);
      strcpy(outFileName, options->out_file);
    }

    char *cloogFileName;
    if (strlen(bname) >= 2 && !strcmp(bname + strlen(bname) - 2, ".c")) {
      cloogFileName =
          (char *)malloc(strlen(bname) - 2 + strlen(".pluto.cloog") + 1);
      strcpy(cloogFileName, bname);
      cloogFileName[strlen(bname) - 2] = '\0';
    } else {
      cloogFileName =
          (char *)malloc(strlen(bname) + strlen(".pluto.cloog") + 1);
      strcpy(cloogFileName, bname);
    }
    strcat(cloogFileName, ".pluto.cloog");
    free(basec);

    FILE *cloogfp = fopen(cloogFileName, "w+");
    if (!cloogfp) {
      fprintf(stderr, "[Pluto] Can't open .cloog file: '%s'\n", cloogFileName);
      free(cloogFileName);
      pluto_options_free(options);
      pluto_prog_free(prog);
      return 9;
    }
    free(cloogFileName);

    FILE *outfp = fopen(outFileName, "w");
    if (!outfp) {
      fprintf(stderr, "[Pluto] Can't open file '%s' for writing\n",
              outFileName);
      free(outFileName);
      pluto_options_free(options);
      pluto_prog_free(prog);
      fclose(cloogfp);
      return 10;
    }

    /* Generate .cloog file */
    pluto_gen_cloog_file(cloogfp, prog);
    /* Add the <irregular> tag from clan, if any */
    if (!options->pet) {
      if (irroption) {
        fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
        free(irroption);
      }
    }
    rewind(cloogfp);

    /* Generate code using Cloog and add necessary stuff before/after code */
    t_start = rtclock();
    pluto_multicore_codegen(cloogfp, outfp, prog);
    t_c = rtclock() - t_start;

    FILE *tmpfp = fopen(".outfilename", "w");
    if (tmpfp) {
      fprintf(tmpfp, "%s\n", outFileName);
      fclose(tmpfp);
      PLUTO_MESSAGE(printf("[Pluto] Output written to %s\n", outFileName););
    }
    free(outFileName);

    fclose(cloogfp);
    fclose(outfp);
  }
  osl_scop_free(scop);

  double t_all = rtclock() - t_start_all;

  if (options->time && !options->silent) {
    printf("\n[pluto] Timing statistics\n[pluto] SCoP extraction + dependence "
           "analysis time: %0.6lfs\n",
           t_d);
    printf("[pluto] Auto-transformation time: %0.6lfs\n", t_t);
    if (options->dfp) {
      printf("[pluto] \tFCG construction time: %0.6lfs\n",
             prog->fcg_const_time);
      printf("[pluto] \tFCG colouring time: %0.6lfs\n", prog->fcg_colour_time);
      printf("[pluto] \tscaling + shifting time: %0.6lfs\n",
             prog->fcg_dims_scale_time);
      printf("[pluto] \tskew determination time: %0.6lfs\n", prog->skew_time);
    }
    printf("[pluto] \t\tTotal constraint solving time (LP/MIP/ILP) time: "
           "%0.6lfs\n",
           prog->mipTime);
    printf("[pluto] Code generation time: %0.6lfs\n", t_c);
    printf("[pluto] Other/Misc time: %0.6lfs\n", t_all - t_c - t_t - t_d);
    printf("[pluto] Total time: %0.6lfs\n", t_all);
    printf("[pluto] All times: %0.6lf %0.6lf %.6lf %.6lf\n", t_d, t_t, t_c,
           t_all - t_c - t_t - t_d);
  }

  pluto_prog_free(prog);
  pluto_options_free(options);

  return 0;
}
