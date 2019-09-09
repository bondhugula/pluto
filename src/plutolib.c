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
 * This file is part of libpluto.
 */

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "isl/map.h"
#include "isl/space.h"
#include "isl/union_set.h"

#include "constraints.h"
#include "math_support.h"
#include "pluto.h"
#include "pluto/pluto.h"
#include "program.h"

PlutoOptions *options = NULL;

/// Run Pluto on the supplied PlutoProg.
static int pluto_schedule_prog(PlutoProg *prog) {
  // Set global var for rest of Pluto.
  options = prog->options;

  unsigned dim_sum = 0;
  for (int i = 0; i < prog->nstmts; i++) {
    dim_sum += prog->stmts[i]->dim;
  }

  if (!options->silent) {
    fprintf(stdout, "[Pluto] Number of statements: %d\n", prog->nstmts);
    fprintf(stdout, "[Pluto] Total number of loops: %u\n", dim_sum);
    fprintf(stdout, "[Pluto] Number of deps: %d\n", prog->ndeps);
    fprintf(stdout, "[Pluto] Maximum domain dimensionality: %d\n", prog->nvar);
    fprintf(stdout, "[Pluto] Number of parameters: %d\n", prog->npar);
  }

  /* Auto transformation */
  if (!options->identity) {
    if (pluto_auto_transform(prog))
      return 1;
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
  return 0;
}

/* Temporary data structure used inside extra_stmt_domains.
 *
 * stmts points to the array of Stmts being constructed
 * index is the index of the next stmt in the array
 */
struct pluto_extra_stmt_info {
  Stmt **stmts;
  unsigned index;
};

static isl_stat extract_basic_set(__isl_take isl_basic_set *bset, void *user) {
  struct pluto_extra_stmt_info *info = (struct pluto_extra_stmt_info *)user;

  Stmt **stmts = info->stmts;
  Stmt *stmt = stmts[info->index];

  PlutoConstraints *bcst = isl_basic_set_to_pluto_constraints(bset);
  if (stmt->domain) {
    stmt->domain = pluto_constraints_unionize_simple(stmt->domain, bcst);
    pluto_constraints_free(bcst);
  } else {
    stmt->domain = bcst;
  }

  isl_basic_set_free(bset);
  return isl_stat_ok;
}

/* Used by libpluto interface. Expects that stmts[idx] have been initialized to
 * NULL.  */
static isl_stat extract_stmt(__isl_take isl_set *set, void *user) {
  Stmt **stmts = (Stmt **)user;

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
  unsigned id = atoi(isl_set_get_tuple_name(set) + 2);

  assert(stmts[id] == NULL && "stmt already initialized - repeated domain?");

  stmts[id] = pluto_stmt_alloc(dim - npar, NULL, trans);

  Stmt *stmt = stmts[id];
  stmt->type = ORIG;
  stmt->id = id;

  for (unsigned i = 0; i < stmt->dim; i++) {
    char *iter = (char *)malloc(13);
    sprintf(iter, "i%d", i);
    stmt->iterators[i] = iter;
  }

  struct pluto_extra_stmt_info info = {stmts, id};
  isl_stat r = isl_set_foreach_basic_set(set, &extract_basic_set, &info);

  pluto_constraints_set_names_range(stmt->domain, stmt->iterators, 0, 0,
                                    stmt->dim);

  for (int i = 0; i < npar; i++) {
    char *param = (char *)malloc(13);
    sprintf(param, "p%d", i);
    stmt->domain->names[stmt->dim + i] = param;
  }

  pluto_matrix_free(trans);

  for (unsigned j = 0; j < stmt->dim; j++) {
    stmt->is_orig_loop[j] = true;
  }

  isl_set_free(set);

  return r;
}

/* Used by libpluto interface */
static void extract_stmts(__isl_keep isl_union_set *domains, PlutoProg *prog) {
  prog->nstmts = isl_union_set_n_set(domains);

  if (prog->nstmts > 0) {
    prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
  } else {
    prog->stmts = NULL;
  }

  for (int i = 0; i < prog->nstmts; i++) {
    prog->stmts[i] = NULL;
  }

  isl_union_set_foreach_set(domains, &extract_stmt, prog->stmts);

  for (int i = 0; i < prog->nstmts; i++)
    assert(prog->stmts[i] && "statement extraction failed; invalid domains?");
}

/*
 * Pluto tiling modifies the domain; move the information from the
 * domain to the schedules (inequalities will be added to the schedules)
 * Polly for LLVM expects domains to remain unchanged
 * (space/dimensionality-wise)
 */
static PlutoConstraints *normalize_domain_schedule(Stmt *stmt,
                                                   PlutoProg *prog) {
  PlutoConstraints *sched = pluto_stmt_get_schedule(stmt);

  PlutoConstraints *domain = stmt->domain;
  unsigned snodes = stmt->dim - stmt->dim_orig;

  while (domain != NULL) {
    unsigned del = 0;
    unsigned nrows = domain->nrows;
    for (unsigned r = 0; r < nrows; r++) {
      unsigned c;
      for (c = 0; c < snodes; c++) {
        if (domain->val[r - del][c] != 0)
          break;
      }
      if (c < snodes) {
        PlutoConstraints *cut = pluto_constraints_select_row(domain, r - del);
        for (unsigned j = 0; j < stmt->trans->nrows; j++) {
          pluto_constraints_add_dim(cut, 0, NULL);
        }
        pluto_constraints_add(sched, cut);
        pluto_constraints_free(cut);
        pluto_constraints_remove_row(domain, r - del);
        del++;
      }
    }
    domain = domain->next;
  }

  for (unsigned c = 0; c < snodes; c++) {
    pluto_stmt_remove_dim(stmt, 0, prog);
  }
  return sched;
}

// Helper structure for callback.
struct stmt_acc_info {
  unsigned stmtId;
  isl_union_map **reads_or_writes;
};

// Callback to extract a specific statement's accesses.
static isl_stat stmt_filter(__isl_take isl_map *map, void *user) {
  const char *name = isl_map_get_tuple_name(map, isl_dim_in);
  if (atoi(name + 2) == ((struct stmt_acc_info *)user)->stmtId) {
    isl_union_map **reads_or_writes =
        ((struct stmt_acc_info *)user)->reads_or_writes;
    *reads_or_writes =
        isl_union_map_union(*reads_or_writes, isl_union_map_from_map(map));
    return isl_stat_ok;
  }
  isl_map_free(map);
  return isl_stat_ok;
}

/// Get schedules in ISL format from PlutoProg's transformations.
// TODO: remove the need from isl_space argument.
static __isl_give isl_union_map *
pluto_schedules_to_isl(PlutoProg *prog, __isl_take isl_space *space,
                       isl_ctx *ctx) {
  /* Construct isl_union_map for pluto schedules */
  isl_union_map *schedules = isl_union_map_empty(isl_space_copy(space));
  for (int i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    PlutoConstraints *sched = normalize_domain_schedule(stmt, prog);

    isl_basic_map *bmap = isl_basic_map_from_pluto_constraints(
        ctx, sched, sched->ncols - stmt->trans->nrows - prog->npar - 1,
        stmt->trans->nrows, prog->npar);
    bmap = isl_basic_map_project_out(bmap, isl_dim_in, 0,
                                     sched->ncols - stmt->trans->nrows -
                                         stmt->domain->ncols);
    char name[20];
    snprintf(name, sizeof(name), "S_%d", i);
    isl_map *map = isl_map_from_basic_map(bmap);

    /* Copy ids of the original parameter dimensions  */
    for (int j = 0, e = isl_space_dim(space, isl_dim_param); j < e; j++) {
      isl_id *id = isl_space_get_dim_id(space, isl_dim_param, j);
      map = isl_map_set_dim_id(map, isl_dim_param, j, id);
    }

    map = isl_map_set_tuple_name(map, isl_dim_in, name);
    schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));

    pluto_constraints_free(sched);
  }
  isl_space_free(space);
  return schedules;
}

/// Construct PlutoProg given domains and dependences in ISL structures.
/// `reads' and `writes' are optional. PlutoAccesses are constructed when they
/// are not set to NULL.
static PlutoProg *pluto_prog_from_isl_domains_dependences(
    __isl_take isl_union_set *domains, __isl_take isl_union_map *dependences,
    __isl_take __isl_take isl_union_map *reads, isl_union_map *writes,
    PlutoOptions *options_l) {
  PlutoProg *prog = pluto_prog_alloc();
  /* Global var */
  options = options_l;
  prog->options = options_l;

  prog->nvar = -1;

  extract_stmts(domains, prog);
  isl_union_set_free(domains);

  for (int i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, (int)prog->stmts[i]->dim);
  }

  if (prog->nstmts > 0) {
    Stmt *stmt = prog->stmts[0];
    prog->npar = stmt->domain->ncols - stmt->dim - 1;
    prog->params = (char **)malloc(sizeof(char *) * prog->npar);
  } else {
    prog->npar = 0;
  }

  for (int i = 0; i < prog->npar; i++) {
    char *param = (char *)malloc(13);
    sprintf(param, "p%d", i);
    prog->params[i] = param;
  }

  // Extract access functions.
  if (reads && writes) {
    for (unsigned i = 0; i < prog->nstmts; ++i) {
      // Extract this statement's reads and writes from 'reads' and 'writes'.
      isl_union_map *reads_s =
          isl_union_map_empty(isl_union_map_get_space(reads));
      isl_union_map *writes_s =
          isl_union_map_empty(isl_union_map_get_space(writes));
      struct stmt_acc_info r_info = {i, &reads_s};
      struct stmt_acc_info w_info = {i, &writes_s};
      isl_union_map_foreach_map(reads, stmt_filter, &r_info);
      isl_union_map_foreach_map(writes, stmt_filter, &w_info);

      // Extract and populate into statements.
      extract_accesses_for_pluto_stmt(prog->stmts[i], reads_s, writes_s);

      isl_union_map_free(reads_s);
      isl_union_map_free(writes_s);
    }
    isl_union_map_free(reads);
    isl_union_map_free(writes);
  }

  prog->ndeps = 0;
  isl_union_map_foreach_map(dependences, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (int i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  extract_deps_from_isl_union_map(dependences, prog->deps, 0, prog->stmts,
                                  PLUTO_DEP_RAW);
  isl_union_map_free(dependences);

  return prog;
}

/// Run the Pluto transformation algorithm on the provided domains and
/// dependences. Read and writes accesses can be optionally provided (NULL
/// otherwise); if they are provided, they are exploited for certain late
/// transformations (for intra-tile optimization in particular). Returns the
/// schedules as an isl_union_map, ownership of which is with the caller.
__isl_give isl_union_map *pluto_transform(__isl_take isl_union_set *domains,
                                          __isl_take isl_union_map *dependences,
                                          __isl_take isl_union_map *reads,
                                          __isl_take isl_union_map *writes,
                                          PlutoOptions *options_l) {
  isl_ctx *ctx = isl_union_set_get_ctx(domains);
  isl_space *space = isl_union_set_get_space(domains);

  PlutoProg *prog = pluto_prog_from_isl_domains_dependences(
      domains, dependences, reads, writes, options_l);

  IF_MORE_DEBUG(printf("Extracted PlutoProg\n"));
  IF_MORE_DEBUG(pluto_prog_print(stdout, prog));

  int retval = pluto_schedule_prog(prog);

  if (retval) {
    /* Failure */
    pluto_prog_free(prog);
    isl_space_free(space);

    if (!options_l->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }
    return NULL;
  }

  isl_union_map *schedules = pluto_schedules_to_isl(prog, space, ctx);

  pluto_prog_free(prog);

  return schedules;
}

/// Use the Pluto transformation algorithm on the domains cum schedules provided
/// in `schedules' and with the read and writes access relations provided in
/// `reads' and `writes' respectively. Returns the schedules as an
/// isl_union_map, ownership of which is with the caller. Note that the returned
/// schedules encode both the mapping and the set information.
__isl_give isl_union_map *pluto_schedule(__isl_take isl_union_map *schedules,
                                         __isl_take isl_union_map *reads,
                                         __isl_take isl_union_map *writes,
                                         PlutoOptions *options_l) {
  isl_union_map *dep_raw, *dep_war, *dep_waw, *dep_rar;

  options = options_l;

  isl_union_map *empty = isl_union_map_empty(isl_union_map_get_space(reads));

  compute_deps_isl(reads, writes, schedules, empty, &dep_raw, &dep_war,
                   &dep_waw, &dep_rar, NULL, NULL);
  isl_union_map_free(empty);

  isl_union_map *dependences = isl_union_map_union(dep_raw, dep_war);
  dependences = isl_union_map_union(dependences, dep_waw);
  isl_union_set *domains = isl_union_map_domain(schedules);

  isl_union_map *pluto_schedules = pluto_transform(
      isl_union_set_copy(domains), dependences, reads, writes, options_l);
  pluto_schedules = isl_union_map_intersect_domain(pluto_schedules, domains);

  return pluto_schedules;
}

/// Construct a remapping matrix for the transformations in prog's statements
/// and return that.
static Remapping *pluto_get_remapping(PlutoProg *prog) {
  Remapping *remapping = (Remapping *)malloc(sizeof(Remapping));
  remapping->nstmts = prog->nstmts;
  remapping->stmt_inv_matrices =
      (PlutoMatrix **)malloc(sizeof(PlutoMatrix *) * prog->nstmts);
  remapping->stmt_divs = (int **)malloc(sizeof(int *) * prog->nstmts);

  for (int i = 0; i < prog->nstmts; i++) {
    remapping->stmt_inv_matrices[i] =
        pluto_stmt_get_remapping(prog->stmts[i], &remapping->stmt_divs[i]);
    IF_DEBUG(printf("[libpluto] Statement %d Id- %d\n", i, prog->stmts[i]->id));
    IF_DEBUG(pluto_matrix_print(stdout, remapping->stmt_inv_matrices[i]));
  }
  return remapping;
}

static Remapping *
pluto_schedule_and_get_remapping(__isl_take isl_union_set *domains,
                                 __isl_take isl_union_map *dependences,
                                 PlutoOptions *options_l) {
  PlutoProg *prog = pluto_prog_from_isl_domains_dependences(
      domains, dependences, NULL, NULL, options_l);

  IF_MORE_DEBUG(printf("Extracted PlutoProg\n"));
  IF_DEBUG(pluto_prog_print(stdout, prog););

  if (pluto_schedule_prog(prog)) {
    /* Failure */
    pluto_prog_free(prog);

    if (!options_l->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }
    return NULL;
  }

  Remapping *remapping = pluto_get_remapping(prog);
  pluto_prog_free(prog);
  return remapping;
}

/// A version of pluto_schedule to get the schedule, parallel loops and
/// remapping all in one shot.
__isl_give isl_union_map *pluto_parallel_schedule_with_remapping(
    __isl_take isl_union_set *domains, __isl_take isl_union_map *dependences,
    Ploop ***ploops, unsigned *nploops, Remapping **remap,
    PlutoOptions *options_l) {
  isl_ctx *ctx = isl_union_set_get_ctx(domains);
  isl_space *space = isl_union_set_get_space(domains);

  PlutoProg *prog = pluto_prog_from_isl_domains_dependences(
      domains, dependences, NULL, NULL, options_l);

  IF_MORE_DEBUG(printf("Extracted PlutoProg\n"));
  IF_DEBUG(pluto_prog_print(stdout, prog));

  if (pluto_schedule_prog(prog)) {
    /* Failure */
    pluto_prog_free(prog);

    if (!options_l->silent) {
      printf("[libpluto] failure, returning NULL schedules\n");
    }
    return NULL;
  }

  Remapping *remapping = pluto_get_remapping(prog);
  *remap = remapping;

  isl_union_map *schedules = pluto_schedules_to_isl(prog, space, ctx);

  pluto_prog_free(prog);
  isl_space_free(space);

  return schedules;
}

/// A wrapper around pluto_parallel_schedule_with_remapping(). This method
/// accepts domain, dependence and PlutoOptions as string and returns a
/// transformed schedule, remapping and parallel loops.
// This function is a HACK. The reason this exists is to allow for easy FFI
// between PolyMage and Pluto. Sending isl objects between PyIsl to libpluto is
// hard (because PyIsl does not seem to have a way to access the underlying C
// object pointer). Hence, the solution is to convert everything to strings,
// and return the generated schedule as a string as well, which is then
// converted back to an ISL object.
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
  *remapping_ptr = remapping;

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
                             const char *dependences_str, PlutoOptions *options,
                             Remapping *remapping) {
  isl_ctx *ctx = isl_ctx_alloc();
  isl_union_set *domains = isl_union_set_read_from_str(ctx, domains_str);
  isl_union_map *dependences =
      isl_union_map_read_from_str(ctx, dependences_str);

  Remapping *remapping_l =
      pluto_schedule_and_get_remapping(domains, dependences, options);
  *remapping = *remapping_l;
  free(remapping_l);
  isl_ctx_free(ctx);
};

void pluto_remapping_free(Remapping remapping) {
  for (int i = 0; i < remapping.nstmts; ++i) {
    pluto_matrix_free(remapping.stmt_inv_matrices[i]);
  }
  free(remapping.stmt_inv_matrices);
  for (int i = 0; i < remapping.nstmts; ++i) {
    free(remapping.stmt_divs[i]);
  }
  free(remapping.stmt_divs);
}

void pluto_schedules_strbuf_free(char *schedules_str_buffer) {
  free(schedules_str_buffer);
}
