/******************************************************************************
 *               libpluto -  A library version of Pluto                       *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2012 Uday Bondhugula                                         *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License version 2.1 as published by the Free Software Foundation.          *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * A copy of the GNU Lesser General Public Licence can be found in the file
 * `LICENSE.LGPL2' in the top-level directory of this distribution.
 *
 */

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "candl/candl.h"
#include "constraints.h"
#include "pluto.h"
#include "pluto/libpluto.h"
#include "program.h"
#include "isl/map.h"
#include "isl/space.h"

#include "candl/scop.h"

PlutoOptions *options;

static double rtclock() {
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday(&Tp, &Tzp);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

/* Temporary data structure used inside extra_stmt_domains
 *
 * stmts points to the array of Stmts being constructed
 * index is the index of the next stmt in the array
 */
struct pluto_extra_stmt_info {
  Stmt **stmts;
  int index;
};

static isl_stat extract_basic_set(__isl_take isl_basic_set *bset, void *user) {
  Stmt **stmts;
  Stmt *stmt;
  PlutoConstraints *bcst;
  struct pluto_extra_stmt_info *info;

  info = (struct pluto_extra_stmt_info *)user;

  stmts = info->stmts;
  stmt = stmts[info->index];

  bcst = isl_basic_set_to_pluto_constraints(bset);
  if (stmt->domain) {
    stmt->domain = pluto_constraints_unionize_simple(stmt->domain, bcst);
    pluto_constraints_free(bcst);
  } else {
    stmt->domain = bcst;
  }

  isl_basic_set_free(bset);
  return isl_stat_ok;
}

/* Used by libpluto interface */
static isl_stat extract_stmt(__isl_take isl_set *set, void *user) {
  isl_stat r;
  Stmt **stmts;
  int id, i;

  stmts = (Stmt **)user;

  int dim = isl_set_dim(set, isl_dim_all);
  int npar = isl_set_dim(set, isl_dim_param);
  PlutoMatrix *trans = pluto_matrix_alloc(dim - npar, dim + 1);
  pluto_matrix_set(trans, 0);
  trans->nrows = 0;

  /* A statement's domain (isl_set) should be named S_%d */
  const char *name = isl_set_get_tuple_name(set);
  assert(name);
  assert(strlen(name) >= 3);
  assert(name[0] == 'S');
  assert(name[1] == '_');
  assert(isdigit(name[2]));
  id = atoi(isl_set_get_tuple_name(set) + 2);

  stmts[id] = pluto_stmt_alloc(dim - npar, NULL, trans);

  Stmt *stmt = stmts[id];
  stmt->type = ORIG;
  stmt->id = id;

  for (i = 0; i < stmt->dim; i++) {
    char *iter = (char *)malloc(13);
    sprintf(iter, "i%d", i);
    stmt->iterators[i] = iter;
  }

  struct pluto_extra_stmt_info info = {stmts, id};
  r = isl_set_foreach_basic_set(set, &extract_basic_set, &info);

  pluto_constraints_set_names_range(stmt->domain, stmt->iterators, 0, 0,
                                    stmt->dim);

  for (i = 0; i < npar; i++) {
    char *param = (char *)malloc(13);
    sprintf(param, "p%d", i);
    stmt->domain->names[stmt->dim + i] = param;
  }

  pluto_matrix_free(trans);

  int j;
  for (j = 0; j < stmt->dim; j++) {
    stmt->is_orig_loop[j] = true;
  }

  isl_set_free(set);

  return r;
}

/* Used by libpluto interface */
static int extract_stmts(__isl_keep isl_union_set *domains, Stmt **stmts) {
  isl_union_set_foreach_set(domains, &extract_stmt, stmts);

  return 0;
}

/*
 * Pluto tiling modifies the domain; move the information from the
 * domain to the schedules (inequalities will be added to the schedules)
 * Polly for LLVM expects domains to remain unchanged
 * (space/dimensionality-wise)
 */
PlutoConstraints *normalize_domain_schedule(Stmt *stmt, PlutoProg *prog) {
  int del, r, c, j, snodes, nrows;
  PlutoConstraints *domain;

  PlutoConstraints *sched = pluto_stmt_get_schedule(stmt);

  domain = stmt->domain;
  snodes = stmt->dim - stmt->dim_orig;

  while (domain != NULL) {
    del = 0;
    nrows = domain->nrows;
    for (r = 0; r < nrows; r++) {
      for (c = 0; c < snodes; c++) {
        if (domain->val[r - del][c] != 0)
          break;
      }
      if (c < snodes) {
        PlutoConstraints *cut = pluto_constraints_select_row(domain, r - del);
        for (j = 0; j < stmt->trans->nrows; j++) {
          pluto_constraints_add_dim(cut, 0, NULL);
        }
        pluto_constraints_add(sched, cut);
        pluto_constraints_remove_row(domain, r - del);
        del++;
      }
    }
    domain = domain->next;
  }

  for (c = 0; c < snodes; c++) {
    pluto_stmt_remove_dim(stmt, 0, prog);
  }
  return sched;
}

/*
 * Output schedules are isl relations that have dims in the order
 * isl_dim_out, isl_dim_in, div, param, const
 */
__isl_give isl_union_map *pluto_schedule(isl_union_set *domains,
                                         isl_union_map *dependences,
                                         PlutoOptions *options_l) {
  int i, j, retval;
  unsigned nbands, n_ibands;
  isl_ctx *ctx;
  isl_space *space;
  double t_t, t_all, t_start;

  ctx = isl_union_set_get_ctx(domains);
  space = isl_union_set_get_space(domains);

  // isl_union_set_dump(domains);
  // isl_union_map_dump(dependences);

  PlutoProg *prog = pluto_prog_alloc();
  prog->options = options_l;

  /* global var */
  options = options_l;

  prog->nvar = -1;
  prog->nstmts = isl_union_set_n_set(domains);

  if (prog->nstmts >= 1) {
    prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
  } else {
    prog->stmts = NULL;
  }

  for (i = 0; i < prog->nstmts; i++) {
    prog->stmts[i] = NULL;
  }

  extract_stmts(domains, prog->stmts);

  for (i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, prog->stmts[i]->dim);
  }

  if (prog->nstmts >= 1) {
    Stmt *stmt = prog->stmts[0];
    prog->npar = stmt->domain->ncols - stmt->dim - 1;
    prog->params = (char **)malloc(sizeof(char *) * prog->npar);
  } else
    prog->npar = 0;

  for (i = 0; i < prog->npar; i++) {
    char *param = (char *)malloc(13);
    sprintf(param, "p%d", i);
    prog->params[i] = param;
  }

  prog->ndeps = 0;
  isl_union_map_foreach_map(dependences, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  extract_deps(prog->deps, 0, prog->stmts, dependences, OSL_DEPENDENCE_RAW);

  IF_DEBUG(pluto_prog_print(stdout, prog););

  t_start = rtclock();
  retval = pluto_auto_transform(prog);
  t_t = rtclock() - t_start;

  if (retval) {
    /* Failure */
    pluto_prog_free(prog);
    isl_space_free(space);

    if (!options->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }

    return NULL;
  }

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(stdout, "[pluto] Affine transformations\n\n");
    /* Print out transformations */
    pluto_transformations_pretty_print(prog);
  }

  Band **bands, **ibands;
  bands = pluto_get_outermost_permutable_bands(prog, &nbands);
  ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
  printf("Outermost tilable bands: %d bands\n", nbands);
  pluto_bands_print(bands, nbands);
  printf("Innermost tilable bands: %d bands\n", n_ibands);
  pluto_bands_print(ibands, n_ibands);

  if (options->tile) {
    pluto_tile(prog);
  } else {
    if (options->intratileopt) {
      pluto_intra_tile_optimize(prog, 0);
    }
  }

  if ((options->parallel) && !options->tile) {
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
      printf("[pluto] WARNING: use --tile for better parallelization \n");
      fprintf(stdout, "[pluto] After skewing:\n");
      pluto_transformations_pretty_print(prog);
      /* IF_DEBUG(pluto_print_hyperplane_properties(prog);); */
    }
  }

  if (options->parallel && !options->silent) {
    unsigned nploops;
    Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);
    printf("[pluto_mark_parallel] %d parallel loops\n", nploops);
    pluto_loops_print(ploops, nploops);
    printf("\n");
    pluto_loops_free(ploops, nploops);
  }

  /* Construct isl_union_map for pluto schedules */
  isl_union_map *schedules = isl_union_map_empty(isl_space_copy(space));
  for (i = 0; i < prog->nstmts; i++) {
    isl_basic_map *bmap;
    isl_map *map;
    Stmt *stmt = prog->stmts[i];
    PlutoConstraints *sched = normalize_domain_schedule(stmt, prog);

    bmap = isl_basic_map_from_pluto_constraints(
        ctx, sched, sched->ncols - stmt->trans->nrows - prog->npar - 1,
        stmt->trans->nrows, prog->npar);
    bmap = isl_basic_map_project_out(bmap, isl_dim_in, 0,
                                     sched->ncols - stmt->trans->nrows -
                                         stmt->domain->ncols);
    char name[20];
    snprintf(name, sizeof(name), "S_%d", i);
    map = isl_map_from_basic_map(bmap);

    /* Copy ids of the original parameter dimensions  */
    for (j = 0; j < isl_space_dim(space, isl_dim_param); j++) {
      isl_id *id = isl_space_get_dim_id(space, isl_dim_param, j);
      map = isl_map_set_dim_id(map, isl_dim_param, j, id);
    }

    map = isl_map_set_tuple_name(map, isl_dim_in, name);
    schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));

    pluto_constraints_free(sched);
  }

  pluto_prog_free(prog);
  isl_space_free(space);

  t_all = rtclock() - t_start;

  if (options->time && !options->silent) {
    printf("[pluto] Auto-transformation time: %0.6lfs\n", t_t);
    printf("[pluto] Other/Misc time: %0.6lfs\n", t_all - t_t);
    printf("[pluto] Total time: %0.6lfs\n", t_all);
  }

  return schedules;
}

Remapping *pluto_get_remapping(isl_union_set *domains,
                               isl_union_map *dependences,
                               PlutoOptions *options_l) {

  unsigned i, nbands, n_ibands;
  int retval;
  isl_space *space;

  space = isl_union_set_get_space(domains);

  // isl_union_set_dump(domains);
  // isl_union_map_dump(dependences);

  PlutoProg *prog = pluto_prog_alloc();
  prog->options = options_l;

  /* global var */
  options = options_l;

  prog->nvar = -1;
  prog->nstmts = isl_union_set_n_set(domains);

  if (prog->nstmts >= 1) {
    prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
  } else {
    prog->stmts = NULL;
  }

  for (i = 0; i < prog->nstmts; i++) {
    prog->stmts[i] = NULL;
  }

  extract_stmts(domains, prog->stmts);

  for (i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, prog->stmts[i]->dim);
  }

  if (prog->nstmts >= 1) {
    Stmt *stmt = prog->stmts[0];
    prog->npar = stmt->domain->ncols - stmt->dim - 1;
    prog->params = (char **)malloc(sizeof(char *) * prog->npar);
  } else
    prog->npar = 0;

  for (i = 0; i < prog->npar; i++) {
    char *param = (char *)malloc(13);
    sprintf(param, "p%d", i);
    prog->params[i] = param;
  }

  prog->ndeps = 0;
  isl_union_map_foreach_map(dependences, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  extract_deps(prog->deps, 0, prog->stmts, dependences, OSL_DEPENDENCE_RAW);

  IF_DEBUG(pluto_prog_print(stdout, prog););

  retval = pluto_auto_transform(prog);

  if (retval) {
    /* Failure */
    pluto_prog_free(prog);
    isl_space_free(space);

    if (!options->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }

    return NULL;
  }

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(stdout, "[pluto] Affine transformations\n\n");
    /* Print out transformations */
    pluto_transformations_pretty_print(prog);
  }

  Band **bands, **ibands;
  bands = pluto_get_outermost_permutable_bands(prog, &nbands);
  ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
  printf("Outermost tilable bands: %d bands\n", nbands);
  pluto_bands_print(bands, nbands);
  printf("Innermost tilable bands: %d bands\n", n_ibands);
  pluto_bands_print(ibands, n_ibands);

  if (options->tile) {
    pluto_tile(prog);
  } else {
    if (options->intratileopt) {
      pluto_intra_tile_optimize(prog, 0);
    }
  }

  Remapping *remapping = (Remapping *)malloc(sizeof(Remapping));
  remapping->nstmts = prog->nstmts;
  remapping->stmt_inv_matrices =
      (PlutoMatrix **)malloc(sizeof(PlutoMatrix *) * prog->nstmts);
  remapping->stmt_divs = (int **)malloc(sizeof(int *) * prog->nstmts);

  for (i = 0; i < prog->nstmts; i++) {
    remapping->stmt_inv_matrices[i] =
        pluto_stmt_get_remapping(prog->stmts[i], &remapping->stmt_divs[i]);
  }
  return remapping;
}

/*
 * Performs pluto-scheduling on an osl_scop.
 * If dependences are not found in Extensions, it'll recalculate them.
 * The osl_scop's statements' domains and scattering are replaced
 * by new ones.
 */
int pluto_schedule_osl(osl_scop_p scop, PlutoOptions *options_l) {
  int i = 0;

  if (!scop || !scop->statement) {
    fprintf(stderr, "Empty Scop passed\n");
    return EXIT_FAILURE;
  }

  options = options_l;

  /* Convert clan scop to Pluto program */
  PlutoProg *prog = scop_to_pluto_prog(scop, options);

  int dim_sum = 0;
  for (i = 0; i < prog->nstmts; i++) {
    dim_sum += prog->stmts[i]->dim;
  }

  /* Make options consistent */
  if (options->multipar == 1 && options->parallel == 0) {
    fprintf(stdout,
            "Warning: multipar needs parallel to be on; turning on parallel\n");
    options->parallel = 1;
  }

  /* Disable pre-vectorization if tile is not on */
  if (options->tile == 0 && options->prevector == 1) {
    /* If code will not be tiled, pre-vectorization does not make
     * sense */
    if (!options->silent) {
      fprintf(
          stdout,
          "[Pluto] Warning: pre-vectorization does not fit (--tile is off)\n");
    }
    options->prevector = 0;
  }

  if (!options->silent) {
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

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(stdout,
            "[Pluto] Affine transformations [<iter coeff's> <const>]\n\n");
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
      printf("[pluto] WARNING: use --tile for better parallelization \n");
      fprintf(stdout, "[pluto] After skewing:\n");
      pluto_transformations_pretty_print(prog);
      /* IF_DEBUG(pluto_print_hyperplane_properties(prog);); */
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

  /* NO MORE TRANSFORMATIONS BEYOND THIS POINT */

  /* Replace the osl_scop's original domains and scatterings
   * by ones newly created by pluto
   */
  pluto_populate_scop(scop, prog, options);

  pluto_prog_free(prog);

  return EXIT_SUCCESS;
}

/* Pluto_schedule method to get schedule, parallel loops and remapping
 *  all in one function
 */
__isl_give isl_union_map *pluto_parallel_schedule_with_remapping(
    isl_union_set *domains, isl_union_map *dependences, Ploop ***ploops,
    unsigned *nploops, Remapping **remap, PlutoOptions *options_l) {
  unsigned i, j, nbands, n_ibands;
  int retval;
  isl_ctx *ctx;
  isl_space *space;
  double t_t, t_all, t_start;

  ctx = isl_union_set_get_ctx(domains);
  space = isl_union_set_get_space(domains);

  PlutoProg *prog = pluto_prog_alloc();
  prog->options = options_l;

  /* global var */
  options = options_l;

  prog->nvar = -1;
  prog->nstmts = isl_union_set_n_set(domains);

  if (prog->nstmts >= 1) {
    prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
  } else {
    prog->stmts = NULL;
  }

  for (i = 0; i < prog->nstmts; i++) {
    prog->stmts[i] = NULL;
  }

  extract_stmts(domains, prog->stmts);

  for (i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, prog->stmts[i]->dim);
  }

  if (prog->nstmts >= 1) {
    Stmt *stmt = prog->stmts[0];
    prog->npar = stmt->domain->ncols - stmt->dim - 1;
    prog->params = (char **)malloc(sizeof(char *) * prog->npar);
  } else
    prog->npar = 0;

  for (i = 0; i < prog->npar; i++) {
    char *param = (char *)malloc(11);
    sprintf(param, "p%d", i);
    prog->params[i] = param;
  }

  prog->ndeps = 0;
  isl_union_map_foreach_map(dependences, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  extract_deps(prog->deps, 0, prog->stmts, dependences, OSL_DEPENDENCE_RAW);

  IF_DEBUG(pluto_prog_print(stdout, prog););

  t_start = rtclock();
  retval = pluto_auto_transform(prog);
  t_t = rtclock() - t_start;

  if (retval) {
    /* Failure */
    pluto_prog_free(prog);
    isl_space_free(space);

    if (!options->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }

    return NULL;
  }

  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(stdout, "[pluto] Affine transformations\n\n");
    /* Print out transformations */
    pluto_transformations_pretty_print(prog);
  }

  Band **bands, **ibands;
  bands = pluto_get_outermost_permutable_bands(prog, &nbands);
  ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
  printf("Outermost tilable bands: %d bands\n", nbands);
  pluto_bands_print(bands, nbands);
  printf("Innermost tilable bands: %d bands\n", n_ibands);
  pluto_bands_print(ibands, n_ibands);

  if (options->tile) {
    pluto_tile(prog);
  } else {
    if (options->intratileopt) {
      pluto_intra_tile_optimize(prog, 0);
    }
  }

  if ((options->parallel) && !options->tile && !options->identity) {
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
      printf("[pluto] WARNING: use --tile for better parallelization \n");
      fprintf(stdout, "[pluto] After skewing:\n");
      pluto_transformations_pretty_print(prog);
      /* IF_DEBUG(pluto_print_hyperplane_properties(prog);); */
    }
  }

  if (options->parallel) {
    *ploops = pluto_get_parallel_loops(prog, nploops);
    if (!options->silent) {
      printf("[pluto_mark_parallel] %d parallel loops\n", *nploops);
      pluto_loops_print(*ploops, *nploops);
      printf("\n");
    }
  } else {
    *nploops = 0;
  }

  // Constructing remapping Matrix
  Remapping *remapping = (Remapping *)malloc(sizeof(Remapping));
  remapping->nstmts = prog->nstmts;
  remapping->stmt_inv_matrices =
      (PlutoMatrix **)malloc(sizeof(PlutoMatrix *) * prog->nstmts);
  remapping->stmt_divs = (int **)malloc(sizeof(int *) * prog->nstmts);

  *remap = remapping;

  for (i = 0; i < prog->nstmts; i++) {
    remapping->stmt_inv_matrices[i] =
        pluto_stmt_get_remapping(prog->stmts[i], &remapping->stmt_divs[i]);
    if (!options->silent) {
      printf("[libpluto] Statement %d Id- %d\n", i, prog->stmts[i]->id);
      pluto_matrix_print(stdout, remapping->stmt_inv_matrices[i]);
    }
  }

  /* Construct isl_union_map for pluto schedules */
  isl_union_map *schedules = isl_union_map_empty(isl_space_copy(space));
  for (i = 0; i < prog->nstmts; i++) {
    isl_basic_map *bmap;
    isl_map *map;
    Stmt *stmt = prog->stmts[i];
    PlutoConstraints *sched = normalize_domain_schedule(stmt, prog);

    bmap = isl_basic_map_from_pluto_constraints(
        ctx, sched, sched->ncols - stmt->trans->nrows - prog->npar - 1,
        stmt->trans->nrows, prog->npar);
    bmap = isl_basic_map_project_out(bmap, isl_dim_in, 0,
                                     sched->ncols - stmt->trans->nrows -
                                         stmt->domain->ncols);
    char name[20];
    snprintf(name, sizeof(name), "S_%d", i);
    map = isl_map_from_basic_map(bmap);

    /* Copy ids of the original parameter dimensions  */
    for (j = 0; j < isl_space_dim(space, isl_dim_param); j++) {
      isl_id *id = isl_space_get_dim_id(space, isl_dim_param, j);
      map = isl_map_set_dim_id(map, isl_dim_param, j, id);
    }

    map = isl_map_set_tuple_name(map, isl_dim_in, name);
    schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));

    pluto_constraints_free(sched);
  }

  pluto_prog_free(prog);
  isl_space_free(space);

  t_all = rtclock() - t_start;

  if (options->time && !options->silent) {
    printf("[pluto] Auto-transformation time: %0.6lfs\n", t_t);
    printf("[pluto] Other/Misc time: %0.6lfs\n", t_all - t_t);
    printf("[pluto] Total time: %0.6lfs\n", t_all);
  }

  return schedules;
}

/* pluto_schedule_str is a wrapper method around
 * pluto_parallel_schedule_with_remapping().
 * This method accepts domain, dependence and PlutoOptions as string
 * and returns a transformed schedule, remapping and parallel loops.
 */
void pluto_schedule_str(const char *domains_str, const char *dependences_str,
                        char **schedules_str_buffer_ptr, char **p_loops,
                        Remapping **remapping_ptr, PlutoOptions *options) {

  isl_ctx *ctx = isl_ctx_alloc();
  Ploop **ploop;
  unsigned nploop = 0, i;
  Remapping *remapping;

  isl_union_set *domains = isl_union_set_read_from_str(ctx, domains_str);
  isl_union_map *dependences =
      isl_union_map_read_from_str(ctx, dependences_str);

  assert(remapping_ptr != NULL);
  isl_union_map *schedule = pluto_parallel_schedule_with_remapping(
      domains, dependences, &ploop, &nploop, &remapping, options);

  if (options->parallel) {
    if (!options->silent) {
      pluto_loops_print(ploop, nploop);
    }

    // NOTE: assuming max 4 digits
    // number of parallel loops
    char *str = (char *)(malloc(sizeof(char) * 5));
    sprintf(str, "%d", nploop);

    // NOTE: assuming max 4 digits per integer, and 1 char for comma
    // 1 place for 'nploop' int itself, and nploop places for the rest
    p_loops[0] = (char *)malloc(sizeof(char) * 6 * (nploop + 1));
    strcpy(p_loops[0], str);

    for (i = 1; i < nploop + 1; i++) {
      // the result is a csv list
      strcat(p_loops[0], ",");
      // add the i'th parallel loop dim
      sprintf(str, "%d", ploop[i - 1]->depth + 1);
      strcat(p_loops[0], str);
    }
  }

  *remapping_ptr = remapping;

  isl_printer *printer = isl_printer_to_str(ctx);
  isl_printer_print_union_map(printer, schedule);

  *schedules_str_buffer_ptr = isl_printer_get_str(printer);
  assert(*schedules_str_buffer_ptr != NULL && "isl printer providing empty"
                                              " string");

  isl_printer_free(printer);
  isl_union_set_free(domains);
  isl_union_map_free(dependences);
  isl_union_map_free(schedule);

  isl_ctx_free(ctx);
}

void pluto_get_remapping_str(const char *domains_str,
                             const char *dependences_str,
                             Remapping **remapping_ptr, PlutoOptions *options) {

  isl_ctx *ctx = isl_ctx_alloc();
  isl_union_set *domains = isl_union_set_read_from_str(ctx, domains_str);
  isl_union_map *dependences =
      isl_union_map_read_from_str(ctx, dependences_str);

  assert(remapping_ptr != NULL);
  Remapping *remapping = pluto_get_remapping(domains, dependences, options);
  *remapping_ptr = remapping;
};

void pluto_remapping_free(Remapping *remapping) {
  assert(remapping != NULL);
  free(remapping);
};

void pluto_schedules_strbuf_free(char *schedules_str_buffer) {
  free(schedules_str_buffer);
}
