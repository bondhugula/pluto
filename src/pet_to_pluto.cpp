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
 * pet_to_pluto.cpp
 *
 * This file contains functions that do the job interfacing the PLUTO
 * core to the frontend and related matters
 *
 */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "constraints.h"
#include "math_support.h"
#include "pluto.h"
#include "program.h"

#include <isl/flow.h>
#include <isl/id.h>
#include <isl/map.h>
#include <isl/mat.h>
#include <isl/set.h>
#include <isl/space.h>
#include <isl/union_map.h>
#include <isl/val.h>

#include "pet.h"

struct acc_info {
  char *prefix;
  int acc_num;
  isl_union_map **new_maps;
  isl_union_map **schedule;
  isl_map *base_schedule;
};

static isl_stat set_tuple_name(__isl_take isl_map *map, void *usr) {
  char *name;

  struct acc_info *info = (struct acc_info *)usr;
  name = (char *)malloc(strlen(info->prefix) + 4);
  sprintf(name, "%s%d", info->prefix, info->acc_num);
  map = isl_map_set_tuple_name(map, isl_dim_in, name);
  info->acc_num++;

  *info->new_maps =
      isl_union_map_union(*info->new_maps, isl_union_map_from_map(map));
  isl_map *schedule_i = isl_map_copy(info->base_schedule);
  schedule_i = isl_map_set_tuple_name(schedule_i, isl_dim_in, name);
  *info->schedule =
      isl_union_map_union(*info->schedule, isl_union_map_from_map(schedule_i));
  free(name);

  return isl_stat_ok;
}

/* Compute dependences based on the domain, scheduling, and access
 * information in "pscop", and put the result in "prog".
 */
static void compute_deps_pet(struct pet_scop *pscop,
                             isl_map **stmt_wise_schedules, PlutoProg *prog,
                             PlutoOptions *options) {
  isl_union_map *dep_raw, *dep_war, *dep_waw, *dep_rar;

  // These are only going to be computed under lastwriter.
  isl_union_map *trans_dep_war = NULL;
  isl_union_map *trans_dep_waw = NULL;

  if (!options->silent) {
    printf("[pluto] compute_deps (isl%s)\n",
           options->lastwriter ? " with lastwriter" : "");
  }

  isl_space *space = isl_set_get_space(pscop->context);
  isl_union_map *empty = isl_union_map_empty(isl_space_copy(space));

  isl_union_map *reads = isl_union_map_copy(empty);
  isl_union_map *writes = isl_union_map_copy(empty);
  isl_union_map *schedule = isl_union_map_copy(empty);

  isl_union_map *schedules = isl_schedule_get_map(pscop->schedule);

  for (int i = 0; i < prog->nstmts; i++) {
    struct pet_stmt *pstmt = prog->stmts[i]->pstmt;
    Stmt *stmt = prog->stmts[i];

    isl_union_map *s_umap = isl_union_map_intersect_domain(
        isl_union_map_copy(schedules),
        isl_union_set_from_set(isl_set_copy(pstmt->domain)));

    isl_map *s_map = isl_map_from_union_map(s_umap);

    isl_union_map *lreads = pet_stmt_collect_accesses(
        pstmt, pet_expr_access_may_read, 0, isl_space_copy(space));
    isl_union_map *lwrites = pet_stmt_collect_accesses(
        pstmt, pet_expr_access_may_write, 0, isl_space_copy(space));

    char name[20];
    sprintf(name, "S_%d_r", stmt->id);
    struct acc_info rinfo = {name, 0, &reads, &schedule, s_map};
    isl_union_map_foreach_map(lreads, &set_tuple_name, &rinfo);
    sprintf(name, "S_%d_w", stmt->id);
    struct acc_info winfo = {name, 0, &writes, &schedule, s_map};
    isl_union_map_foreach_map(lwrites, &set_tuple_name, &winfo);

    isl_map_free(s_map);
    isl_union_map_free(lreads);
    isl_union_map_free(lwrites);
  }
  isl_union_map_free(schedules);
  isl_space_free(space);

  compute_deps_isl(reads, writes, schedule, empty, &dep_raw, &dep_war, &dep_waw,
                   &dep_rar, &trans_dep_war, &trans_dep_waw);

  prog->ndeps = 0;
  isl_union_map_foreach_map(dep_raw, &isl_map_count, &prog->ndeps);
  isl_union_map_foreach_map(dep_war, &isl_map_count, &prog->ndeps);
  isl_union_map_foreach_map(dep_waw, &isl_map_count, &prog->ndeps);

  prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
  for (int i = 0; i < prog->ndeps; i++) {
    prog->deps[i] = pluto_dep_alloc();
  }
  prog->ndeps = 0;
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_raw,
                              OSL_DEPENDENCE_RAW);
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_war,
                              OSL_DEPENDENCE_WAR);
  prog->ndeps += extract_deps(prog->deps, prog->ndeps, prog->stmts, dep_waw,
                              OSL_DEPENDENCE_WAW);
  prog->transdeps = NULL;
  prog->ntransdeps = 0;

  isl_union_map_free(dep_raw);
  isl_union_map_free(dep_war);
  isl_union_map_free(dep_waw);
  isl_union_map_free(trans_dep_waw);
  isl_union_map_free(trans_dep_war);

  isl_union_map_free(writes);
  isl_union_map_free(reads);
  isl_union_map_free(schedule);
  isl_union_map_free(empty);
}

/* Removes certain trivial dead code - in particular, all writes to variables
 * that have been marked as killed through special kill statements. */
static void mark_trivial_dead_code(struct pet_scop *pscop,
                                   std::vector<bool> *dead) {
  dead->resize(pscop->n_stmt, false);
  for (int s = 0; s < pscop->n_stmt; s++) {
    struct pet_stmt *pstmt = pscop->stmts[s];
    if (!pet_stmt_is_kill(pstmt)) {
      continue;
    }
    (*dead)[s] = true;
    isl_space *space = isl_set_get_space(pscop->context);
    isl_union_map *writes_s =
        pet_stmt_collect_accesses(pstmt, pet_expr_access_killed, 0, space);

    if (isl_union_map_n_map(writes_s) != 1) {
      isl_union_map_free(writes_s);
      continue;
    }

    isl_map *write_s = isl_map_from_union_map(writes_s);
    isl_space *acc_space_s = isl_map_get_space(write_s);
    const char *killed_name =
        isl_space_get_tuple_name(acc_space_s, isl_dim_out);
    isl_map_free(write_s);

    // Mark any other writes to the same variable name dead.
    // This is a HACK to get rid of old IV init's and increments.
    for (int j = 0; j < pscop->n_stmt; j++) {
      if ((*dead)[j])
        continue;
      struct pet_stmt *opstmt = pscop->stmts[j];
      isl_space *space = isl_set_get_space(pscop->context);
      isl_union_map *writes = pet_stmt_collect_accesses(
          opstmt, pet_expr_access_may_write, 0, space);

      if (isl_union_map_n_map(writes) != 1) {
        isl_union_map_free(writes);
        continue;
      }
      isl_map *write = isl_map_from_union_map(writes);
      isl_space *acc_space = isl_map_get_space(write);

      const char *name = isl_space_get_tuple_name(acc_space, isl_dim_out);
      if (!strcmp(name, killed_name))
        (*dead)[j] = true;
      isl_space_free(acc_space);
      isl_map_free(write);
    }
    isl_space_free(acc_space_s);
  }
}

/* Read statement info from pet structures (nvar: max domain dim) */
static Stmt **pet_to_pluto_stmts(
    struct pet_scop *pscop, isl_map **stmt_wise_schedules,
    const std::unordered_map<struct pet_stmt *, char *> &stmtTextMap,
    int *nstmts) {
  int i, j, s, t;
  Stmt **stmts;
  int nvar, npar, max_sched_rows;
  char **params;

  IF_DEBUG(printf("[pluto] Pet to Pluto stmts\n"););

  npar = isl_set_dim(pscop->context, isl_dim_all);
  *nstmts = pscop->n_stmt;

  /* This takes cares of marking trivial statements such as original iterator
   * assignments and increments as dead code */
  std::vector<bool> dead;
  mark_trivial_dead_code(pscop, &dead);

  if (*nstmts == 0)
    return NULL;

  IF_DEBUG(printf("[pluto] Pet SCoP context\n"););
  IF_DEBUG(isl_set_dump(pscop->context););

  params = NULL;
  if (npar >= 1) {
    params = (char **)malloc(sizeof(char *) * npar);
  }
  isl_space *cspace = isl_set_get_space(pscop->context);
  for (i = 0; i < npar; i++) {
    params[i] = strdup(isl_space_get_dim_name(cspace, isl_dim_param, i));
  }
  isl_space_free(cspace);

  /* Max dom dimensionality */
  nvar = -1;

  *nstmts = 0;

  for (s = 0; s < pscop->n_stmt; s++) {
    struct pet_stmt *pstmt = pscop->stmts[s];
    int stmt_dim = isl_set_dim(pstmt->domain, isl_dim_set);
    nvar = PLMAX(nvar, stmt_dim);
    if (!dead[s])
      (*nstmts)++;
  }

  stmts = (Stmt **)malloc(*nstmts * sizeof(Stmt *));

  isl_union_map *s_umap = isl_schedule_get_map(pscop->schedule);

  max_sched_rows = 0;

  for (s = 0, t = 0; s < pscop->n_stmt; s++) {
    if (dead[s])
      continue;
    struct pet_stmt *pstmt = pscop->stmts[s];
    PlutoConstraints *domain = isl_set_to_pluto_constraints(pstmt->domain);

    isl_map *s_map = isl_map_from_union_map(isl_union_map_intersect_domain(
        isl_union_map_copy(s_umap),
        isl_union_set_from_set(isl_set_copy(pstmt->domain))));

    int nrows = isl_map_dim(s_map, isl_dim_out);
    max_sched_rows = PLMAX(max_sched_rows, nrows);

    PlutoMatrix *trans = isl_map_to_pluto_func(
        s_map, isl_set_dim(pstmt->domain, isl_dim_set), npar);

    isl_map_free(s_map);

    stmts[t] = pluto_stmt_alloc(isl_set_dim(pstmt->domain, isl_dim_set), domain,
                                trans);
    pluto_constraints_free(domain);
    pluto_matrix_free(trans);

    Stmt *stmt = stmts[t];

    stmt->id = t;
    stmt->type = ORIG;

    for (unsigned j = 0; j < stmt->dim; j++) {
      stmt->is_orig_loop[j] = true;
    }

    /* Tile it if it's tilable unless turned off by .fst/.precut file */
    stmt->tile = 1;

    /* Store the iterator names*/
    isl_space *dspace = isl_set_get_space(pstmt->domain);
    for (unsigned j = 0; j < stmt->dim; j++) {
      stmt->iterators[j] =
          strdup(isl_space_get_dim_name(dspace, isl_dim_set, j));
    }
    isl_space_free(dspace);

    pluto_constraints_set_names_range(stmt->domain, stmt->iterators, 0, 0,
                                      stmt->dim);
    pluto_constraints_set_names_range(stmt->domain, params, stmt->dim, 0, npar);

    // Copy the body of the statement found by print_user.
    const auto &entry = stmtTextMap.find(pstmt);
    if (entry != stmtTextMap.end()) {
      stmt->text = strdup(entry->second);
      // The string from the ISL printer has a trailing new line; wipe that out.
      stmt->text[strlen(stmt->text) - 1] = '\0';
    } else {
      stmt->text = strdup("/* kill statement */");
    }

    isl_space *space = isl_set_get_space(pscop->context);
    isl_union_map *reads = pet_stmt_collect_accesses(
        pstmt, pet_expr_access_may_read, 0, isl_space_copy(space));
    isl_union_map *writes =
        pet_stmt_collect_accesses(pstmt, pet_expr_access_may_write, 0, space);
    isl_union_map_foreach_map(reads, &isl_map_count, &stmt->nreads);
    isl_union_map_foreach_map(writes, &isl_map_count, &stmt->nwrites);

    struct pluto_access_meta_info e_reads = {&stmt->reads, 0, stmt->dim, npar};
    struct pluto_access_meta_info e_writes = {&stmt->writes, 0, stmt->dim,
                                              npar};

    if (stmt->nreads >= 1) {
      stmt->reads =
          (PlutoAccess **)malloc(stmt->nreads * sizeof(PlutoAccess *));
    }
    if (stmt->nwrites >= 1) {
      stmt->writes =
          (PlutoAccess **)malloc(stmt->nwrites * sizeof(PlutoAccess *));
    }
    for (j = 0; j < stmt->nreads; j++) {
      stmt->reads[j] = NULL;
    }
    for (j = 0; j < stmt->nwrites; j++) {
      stmt->writes[j] = NULL;
    }

    isl_union_map_foreach_map(reads, &isl_map_extract_access_func, &e_reads);
    isl_union_map_foreach_map(writes, &isl_map_extract_access_func, &e_writes);

    isl_union_map_free(reads);
    isl_union_map_free(writes);

    stmts[t]->pstmt = pstmt;
    t++;
  }

  isl_union_map_free(s_umap);

  for (s = 0; s < *nstmts; s++) {
    /* Pad with all zero rows */
    int curr_sched_rows = stmts[s]->trans->nrows;
    for (j = curr_sched_rows; j < max_sched_rows; j++) {
      pluto_stmt_add_hyperplane(stmts[s], H_SCALAR, j);
    }
  }

  for (j = 0; j < npar; j++) {
    free(params[j]);
  }
  free(params);

  return stmts;
}

/* Find the element in scop->stmts that has the given "name".
 *  */
static struct pet_stmt *find_stmt(struct pet_scop *scop, const char *name) {
  int i;

  for (i = 0; i < scop->n_stmt; ++i) {
    struct pet_stmt *stmt = scop->stmts[i];
    const char *name_i;

    name_i = isl_set_get_tuple_name(stmt->domain);
    if (!strcmp(name, name_i))
      return stmt;
  }
  return NULL;
}

/* Find the element in scop->stmts that the same name
 *  * as the function call by the given user node.
 *   * These names are determined by the names of the domains
 *    * of the schedule constructed in transform().
 *     */
static struct pet_stmt *extract_pet_stmt(__isl_keep isl_ast_node *node,
                                         struct pet_scop *scop) {
  isl_ast_expr *expr, *arg;
  isl_id *id;
  struct pet_stmt *stmt;

  expr = isl_ast_node_user_get_expr(node);
  arg = isl_ast_expr_get_op_arg(expr, 0);
  isl_ast_expr_free(expr);
  id = isl_ast_expr_get_id(arg);
  isl_ast_expr_free(arg);
  stmt = find_stmt(scop, isl_id_get_name(id));
  isl_id_free(id);

  return stmt;
}

/* Index transformation callback for pet_stmt_build_ast_exprs.
 * "index" expresses the array indices in terms of statement iterators
 * "iterator_map" expresses the statement iterators in terms of
 * AST loop iterators.
 *
 * The result expresses the array indices in terms of
 * AST loop iterators.
 */
static __isl_give isl_multi_pw_aff *
pullback_index(__isl_take isl_multi_pw_aff *index, __isl_keep isl_id *id,
               void *user) {
  isl_pw_multi_aff *iterator_map = (isl_pw_multi_aff *)user;

  iterator_map = isl_pw_multi_aff_copy(iterator_map);
  return isl_multi_pw_aff_pullback_pw_multi_aff(index, iterator_map);
}

/* Set the iterator names using schedule map of the statement*/
static __isl_give isl_id_list *generate_names(isl_ctx *ctx,
                                              struct pet_stmt *stmt) {
  int i;
  isl_id_list *names;
  isl_id *id;

  names = isl_id_list_alloc(ctx, isl_set_dim(stmt->domain, isl_dim_set));

  for (i = 0; i < isl_set_dim(stmt->domain, isl_dim_set); i++) {
    id = isl_id_alloc(ctx, isl_set_get_dim_name(stmt->domain, isl_dim_set, i),
                      NULL);
    names = isl_id_list_add(names, id);
  }

  return names;
}

static __isl_give void free_isl_id_to_ast_expr(void *user) {
  isl_id_to_ast_expr_free((isl_id_to_ast_expr *)user);
}

/* Transform the accesses in the statement associated to the domain
 * called by "node" to refer to the AST loop iterators, construct
 * corresponding AST expressions using "build" and attach them
 * to the node.
 */
static __isl_give isl_ast_node *at_each_domain(__isl_take isl_ast_node *node,
                                               __isl_keep isl_ast_build *build,
                                               void *user) {
  struct pet_stmt *stmt;
  isl_ctx *ctx;
  isl_id *id;
  isl_map *map;
  isl_pw_multi_aff *iterator_map;
  isl_id_to_ast_expr *ref2expr;
  struct pet_scop *scop = (struct pet_scop *)user;

  ctx = isl_ast_node_get_ctx(node);

  stmt = extract_pet_stmt(node, scop);
  if (!stmt)
    isl_die(ctx, isl_error_internal, "cannot find statement",
            isl_ast_node_free(node);
            node = NULL);

  map = isl_map_from_union_map(isl_ast_build_get_schedule(build));
  map = isl_map_reverse(map);
  iterator_map = isl_pw_multi_aff_from_map(map);

  isl_id_list *iterators = generate_names(ctx, stmt);
  build = isl_ast_build_set_iterators(build, iterators);

  ref2expr = pet_stmt_build_ast_exprs(stmt, build, &pullback_index,
                                      iterator_map, NULL, NULL);
  isl_pw_multi_aff_free(iterator_map);

  id = isl_id_alloc(ctx, NULL, ref2expr);
  id = isl_id_set_free_user(id, &free_isl_id_to_ast_expr);

  return isl_ast_node_set_annotation(node, id);
}

struct print_stmt_user_info {
  struct pet_scop *scop;
  // A map to hold source text corresponding to the statement.
  std::unordered_map<struct pet_stmt *, char *> *stmtTextMap;
};

/*
 * Print the statement corresponding to "node" to "p", and also set
 * pet_stmt's text to that.
 *
 * We look for the statement in the pet_scop passed through "user".
 * The AST expressions for all references in the statement
 * have been attached to the node by at_each_domain().
 *
 * Note that p may either be a file printer or string printer
 */
static __isl_give isl_printer *
print_stmt(__isl_take isl_printer *p,
           __isl_take isl_ast_print_options *print_options,
           __isl_keep isl_ast_node *node, void *user) {
  isl_id_to_ast_expr *ref2expr;
  isl_id *id;
  struct pet_stmt *pstmt;
  struct print_stmt_user_info *info = (struct print_stmt_user_info *)user;
  struct pet_scop *scop = info->scop;
  auto *stmtTextMap = info->stmtTextMap;

  /* Printer just for the stmt text */
  isl_printer *p_l;

  p_l = isl_printer_to_str(isl_printer_get_ctx(p));
  p_l = isl_printer_set_output_format(p_l, ISL_FORMAT_C);

  pstmt = extract_pet_stmt(node, scop);

  id = isl_ast_node_get_annotation(node);
  ref2expr = (isl_id_to_ast_expr *)isl_id_get_user(id);
  isl_id_free(id);

  /* Print to both p and p_l. */
  /* Printing to 'p' is just for debugging purposes - so that we could see the
   * AST */
  p = pet_stmt_print_body(pstmt, p, ref2expr);
  p_l = pet_stmt_print_body(pstmt, p_l, ref2expr);
  (*stmtTextMap)[pstmt] = isl_printer_get_str(p_l);
  isl_printer_free(p_l);

  isl_ast_print_options_free(print_options);

  return p;
}

/*
 * Collect the iteration domains of the statements in "scop",
 * skipping kill statements.
 */
static __isl_give isl_union_set *
collect_non_kill_domains(struct pet_scop *scop,
                         int (*pred)(struct pet_stmt *stmt)) {
  int i;
  isl_set *domain_i;
  isl_union_set *domain;

  if (!scop)
    return NULL;

  domain = isl_union_set_empty(isl_set_get_space(scop->context));

  for (i = 0; i < scop->n_stmt; ++i) {
    struct pet_stmt *stmt = scop->stmts[i];

    if (pred(stmt))
      continue;

    if (stmt->n_arg > 0)
      isl_die(isl_union_set_get_ctx(domain), isl_error_unsupported,
              "data dependent conditions not supported",
              return isl_union_set_free(domain));

    domain_i = isl_set_copy(scop->stmts[i]->domain);
    domain = isl_union_set_add_set(domain, domain_i);
  }

  return domain;
}

// Code generate the scop 'scop' and print the corresponding C code to 'p'.
static __isl_give isl_printer *construct_stmt_body(
    struct pet_scop *scop, __isl_take isl_printer *p,
    std::unordered_map<struct pet_stmt *, char *> *stmtTextMap) {
  isl_ctx *ctx = isl_printer_get_ctx(p);
  isl_union_set *domain_set;
  isl_ast_build *build;
  isl_ast_print_options *print_options;
  isl_ast_node *tree;

  domain_set = collect_non_kill_domains(scop, &pet_stmt_is_kill);
  isl_schedule *sched = isl_schedule_intersect_domain(
      isl_schedule_copy(scop->schedule), domain_set);

  build = isl_ast_build_from_context(isl_set_copy(scop->context));
  build = isl_ast_build_set_at_each_domain(build, &at_each_domain, scop);

  tree = isl_ast_build_node_from_schedule(build, sched);
  isl_ast_build_free(build);

  print_options = isl_ast_print_options_alloc(ctx);

  struct print_stmt_user_info info = {scop, stmtTextMap};
  print_options =
      isl_ast_print_options_set_print_user(print_options, &print_stmt, &info);
  p = isl_ast_node_print(tree, p, print_options);

  isl_ast_node_free(tree);

  return p;
}

/*
 * Extract necessary information from pet_scop to create PlutoProg - a
 * representation of the program sufficient to be used throughout Pluto.
 * PlutoProg also includes dependences; uses isl.
 */
PlutoProg *pet_to_pluto_prog(struct pet_scop *pscop, isl_ctx *ctx,
                             PlutoOptions *options) {
  int i, max_sched_rows, npar;

  if (pscop == NULL)
    return NULL;

  pet_scop_align_params(pscop);

  PlutoProg *prog = pluto_prog_alloc();

  /* Program parameters */
  npar = isl_set_dim(pscop->context, isl_dim_all);

  isl_space *cspace = isl_set_get_space(pscop->context);
  for (i = 0; i < npar; i++) {
    pluto_prog_add_param(prog, isl_space_get_dim_name(cspace, isl_dim_param, i),
                         prog->npar);
  }
  isl_space_free(cspace);

  pluto_constraints_free(prog->context);
  prog->context = isl_set_to_pluto_constraints(pscop->context);
  IF_DEBUG(printf("[pluto] Pet SCoP context\n"));
  IF_DEBUG(isl_set_dump(pscop->context););
  IF_DEBUG(pluto_constraints_compact_print(stdout, prog->context));

  if (options->codegen_context != -1) {
    for (i = 0; i < prog->npar; i++) {
      pluto_constraints_add_inequality(prog->codegen_context);
      prog->codegen_context->val[i][i] = 1;
      prog->codegen_context->val[i][prog->codegen_context->ncols - 1] =
          -options->codegen_context;
    }
  }

  read_codegen_context_from_file(prog->codegen_context);

  prog->options = options;

  prog->nvar = -1;
  max_sched_rows = 0;

  for (i = 0; i < pscop->n_stmt; i++) {
    struct pet_stmt *pstmt = pscop->stmts[i];

    int stmt_dim = isl_set_dim(pstmt->domain, isl_dim_set);
    prog->nvar = PLMAX(prog->nvar, stmt_dim);

    isl_union_map *sched_map = isl_schedule_get_map(pscop->schedule);
    isl_union_map *stmt_sched_umap = isl_union_map_intersect_domain(
        sched_map, isl_union_set_from_set(isl_set_copy(pstmt->domain)));
    isl_map *stmt_sched_map = isl_map_from_union_map(stmt_sched_umap);

    int nrows = isl_map_dim(stmt_sched_map, isl_dim_out);
    max_sched_rows = PLMAX(max_sched_rows, nrows);
    isl_map_free(stmt_sched_map);
  }

  isl_printer *p = isl_printer_to_str(ctx);
  // A map to hold source text corresponding to statements.
  std::unordered_map<struct pet_stmt *, char *> stmtTextMap;
  p = construct_stmt_body(pscop, p, &stmtTextMap);
  isl_printer_free(p);

  prog->stmts = pet_to_pluto_stmts(pscop, NULL, stmtTextMap, &prog->nstmts);

  // Free the source strings.
  for (const auto &entry : stmtTextMap) {
    free(entry.second);
  }

  /* Compute dependences */
  compute_deps_pet(pscop, NULL, prog, options);

  /* Add hyperplanes */
  if (prog->nstmts >= 1) {
    for (i = 0; i < max_sched_rows; i++) {
      pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_UNKNOWN);
      prog->hProps[prog->num_hyperplanes - 1].type =
          (i % 2) ? H_LOOP : H_SCALAR;
    }
  }

  /* Hack for linearized accesses */
  FILE *lfp = fopen(".linearized", "r");
  FILE *nlfp = fopen(".nonlinearized", "r");
  char tmpstr[256];
  char linearized[256];
  if (lfp && nlfp) {
    for (i = 0; i < prog->nstmts; i++) {
      rewind(lfp);
      rewind(nlfp);
      while (!feof(lfp) && !feof(nlfp)) {
        fgets(tmpstr, 256, nlfp);
        fgets(linearized, 256, lfp);
        if (strstr(tmpstr, prog->stmts[i]->text)) {
          prog->stmts[i]->text = (char *)realloc(
              prog->stmts[i]->text, sizeof(char) * (strlen(linearized) + 1));
          strcpy(prog->stmts[i]->text, linearized);
        }
      }
    }
    fclose(lfp);
    fclose(nlfp);
  }

  return prog;
}
