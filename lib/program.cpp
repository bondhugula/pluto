/*
 * Pluto: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE in the
 * top-level directory for details.
 *
 * This file is part of libpluto.
 *
 * This file contains functions that do the job interfacing the PLUTO
 * core to the frontend and related matters
 *
 */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "constraints.h"
#include "ddg.h"
#include "isl_support.h"
#include "math_support.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "program.h"

#include <isl/aff.h>
#include <isl/flow.h>
#include <isl/map.h>
#include <isl/space.h>
#include <isl/union_map.h>

void pluto_add_dep(PlutoProg *prog, Dep *dep) {
  dep->id = prog->ndeps;
  prog->ndeps++;
  prog->deps = (Dep **)realloc(prog->deps, sizeof(Dep *) * prog->ndeps);
  prog->deps[prog->ndeps - 1] = dep;
}

void pluto_dep_print(FILE *fp, const Dep *dep) {
  fprintf(fp,
          "--- Dep %d from S%d to S%d; satisfied: %d, sat level: %d; Type: ",
          dep->id + 1, dep->src + 1, dep->dest + 1, dep->satisfied,
          dep->satisfaction_level);

  switch (dep->type) {
  case PLUTO_DEP_UNDEFINED:
    fprintf(fp, "UNSET");
    break;
  case PLUTO_DEP_RAW:
    fprintf(fp, "RAW");
    break;
  case PLUTO_DEP_WAR:
    fprintf(fp, "WAR");
    break;
  case PLUTO_DEP_WAW:
    fprintf(fp, "WAW");
    break;
  case PLUTO_DEP_RAR:
    fprintf(fp, "RAR");
    break;
  default:
    fprintf(fp, "unknown");
    break;
  }

  fprintf(fp, "\n");
  if (dep->src_acc) {
    fprintf(fp, "on variable: %s\n", dep->src_acc->name);
  }

  fprintf(fp, "Dependence polyhedron\n");
  pluto_constraints_compact_print(fp, dep->dpolytope);
  fprintf(fp, "\n");
}

void pluto_deps_print(FILE *fp, PlutoProg *prog) {
  int i;
  if (prog->ndeps == 0)
    printf("** No dependences **\n\n");
  for (i = 0; i < prog->ndeps; i++) {
    pluto_dep_print(fp, prog->deps[i]);
  }
}

void pluto_access_print(FILE *fp, const PlutoAccess *acc, const Stmt *stmt) {
  int j, npar;

  if (!acc) {
    fprintf(fp, "access is NULL\n");
    return;
  }

  npar = acc->mat->ncols - (int)stmt->dim - 1;

  fprintf(fp, "%s", acc->name);
  for (unsigned i = 0; i < acc->mat->nrows; i++) {
    fprintf(fp, "[");
    const char **vars =
        (const char **)malloc((stmt->dim + npar) * sizeof(char *));
    for (unsigned j = 0; j < stmt->dim; j++) {
      vars[j] = stmt->iterators[j];
    }
    for (j = 0; j < npar; j++) {
      if (stmt->domain->names && stmt->domain->names[stmt->dim + j]) {
        vars[stmt->dim + j] = stmt->domain->names[stmt->dim + j];
      } else {
        vars[stmt->dim + j] = "p?";
      }
    }
    pluto_affine_function_print(stdout, acc->mat->val[i], stmt->dim + npar,
                                vars);
    fprintf(fp, "]");
    free(vars);
  }
  fprintf(fp, "\n");
}

void pluto_stmt_print(FILE *fp, const Stmt *stmt) {
  fprintf(fp, "S%d \"%s\"\n", stmt->id + 1, stmt->text);

  fprintf(fp, "ndims: %d; orig_depth: %d\n", stmt->dim, stmt->dim_orig);
  fprintf(fp, "Index set\n");
  pluto_constraints_compact_print(fp, stmt->domain);
  pluto_stmt_transformation_print(stmt);

  if (stmt->nreads == 0) {
    fprintf(fp, "No Read accesses\n");
  } else {
    fprintf(fp, "Read accesses\n");
    for (int i = 0; i < stmt->nreads; i++) {
      pluto_access_print(fp, stmt->reads[i], stmt);
    }
  }

  if (stmt->nwrites == 0) {
    fprintf(fp, "No write access\n");
  } else {
    fprintf(fp, "Write accesses\n");
    for (int i = 0; i < stmt->nwrites; i++) {
      pluto_access_print(fp, stmt->writes[i], stmt);
    }
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    printf("Original loop: %d -> %s\n", i,
           stmt->is_orig_loop[i] ? "yes" : "no");
  }

  fprintf(fp, "\n");
}

/*
 * Checks whether the transformation at 'level' corresponds to a tiled
 * dimension (an affine function of domain dimension divided by a constant),
 * and if yes, returns the function that was tiled (as a function of the
 * domain iterators/param) aka tiling hyperplane, tile_size will point to the
 * tile size / divisor; if no, returns NULL
 *
 * t_i = z is a tile dimension if there exist the following two constraints
 * in the domain:
 *
 * T*z <= f(i_S, p) <= T*z + T - 1, where T is the tile size / divisor
 *
 * Thus, t_i = f(i_S, p)/T
 *
 * pos: position of the supernode in the domain
 *
 */
static int64_t *pluto_check_supernode(const Stmt *stmt, unsigned pos,
                                      int *tile_size) {
  int lb_pos, ub_pos;
  int64_t *tile_hyp;

  PlutoConstraints *dom = stmt->domain;

  lb_pos = -1;
  ub_pos = -1;
  *tile_size = -1;

  for (unsigned r = 0; r < dom->nrows; r++) {
    if (dom->val[r][pos] >= 1 &&
        dom->val[r][dom->ncols - 1] == dom->val[r][pos] - 1) {
      ub_pos = r;
      *tile_size = dom->val[r][pos];
    }

    if (dom->val[r][pos] <= -1 && dom->val[r][dom->ncols - 1] == 0) {
      lb_pos = r;
    }
  }

  if (ub_pos == -1 || lb_pos == -1)
    return NULL;

  int c;
  for (c = 0; c < (int)dom->ncols - 1; c++) {
    if (dom->val[ub_pos][c] != -dom->val[lb_pos][c])
      break;
  }
  if (c < (int)dom->ncols - 1)
    return NULL;

  tile_hyp = (int64_t *)malloc(dom->ncols * sizeof(int64_t));

  for (unsigned c = 0; c < dom->ncols; c++) {
    if (c == pos)
      tile_hyp[c] = 0;
    else
      tile_hyp[c] = dom->val[lb_pos][c];
  }

  return tile_hyp;
}

/// Returns true if the input loop is a tile space loop.
bool is_tile_space_loop(Ploop *loop, const PlutoProg *prog) {
  int tile_size;
  int64_t *tiling_hyp =
      pluto_check_supernode(loop->stmts[0], loop->depth, &tile_size);
  if (tiling_hyp == NULL) {
    return false;
  }
  free(tiling_hyp);
  return true;
}

static int is_skewed(int64_t *func, int len) {
  int count, i;

  count = 0;

  for (i = 0; i < len; i++) {
    if (func[i] != 0)
      count++;
  }

  return count <= 1 ? 0 : 1;
}

/*
 * Prints the 1-d affine function/hyperplane for 'stmt' at depth 'level'
 */
void pluto_stmt_print_hyperplane(FILE *fp, const Stmt *stmt, int level) {
  int npar, j;

  npar = stmt->domain->ncols - (int)stmt->dim - 1;

  char **vars = (char **)malloc((stmt->dim + npar) * sizeof(char *));

  for (unsigned j = 0; j < stmt->dim; j++) {
    vars[j] = strdup(stmt->iterators[j]);
  }
  for (j = 0; j < npar; j++) {
    if (stmt->domain->names && stmt->domain->names[stmt->dim + j]) {
      vars[stmt->dim + j] = stmt->domain->names[stmt->dim + j];
    } else {
      vars[stmt->dim + j] = (char *)"p?";
    }
  }

  for (unsigned j = 0; j < stmt->dim; j++) {
    /* Detect if this dimension is an affine function of other dimensions
     * divided by a constant -- useful to print tiled hyperplanes, the
     * dividing constant being the tile size */
    int div;
    int64_t *super_func;
    super_func = pluto_check_supernode(stmt, j, &div);
    if (super_func) {
      char *tmp;
      tmp = pluto_affine_function_sprint(super_func, stmt->dim + npar,
                                         (const char **)vars);
      free(vars[j]);
      vars[j] = tmp;
      if (is_skewed(super_func, stmt->domain->ncols)) {
        vars[j] = (char *)realloc(vars[j],
                                  1 + strlen(vars[j]) + 2 + log10(div) + 1 + 1);
        sprintf(vars[j], "(%s", tmp = strdup(vars[j]));
        free(tmp);
        sprintf(vars[j] + strlen(vars[j]), ")/%d", div);
      } else {
        vars[j] =
            (char *)realloc(vars[j], strlen(vars[j]) + 1 + log10(div) + 1 + 1);
        sprintf(vars[j] + strlen(vars[j]), "/%d", div);
      }
      free(super_func);
    }
  }

  pluto_affine_function_print(fp, stmt->trans->val[level], stmt->dim + npar,
                              (const char **)vars);

  for (unsigned j = 0; j < stmt->dim; j++) {
    free(vars[j]);
  }

  free(vars);
}

void pluto_stmts_print(FILE *fp, Stmt **stmts, int nstmts) {
  int i;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_print(fp, stmts[i]);
  }
}

void pluto_prog_print(FILE *fp, PlutoProg *prog) {
  int i;

  fprintf(fp, "nvar = %d, npar = %d\n", prog->nvar, prog->npar);
  fprintf(fp, "Parameters: ");

  for (i = 0; i < prog->npar; i++) {
    fprintf(fp, "%s ", prog->params[i]);
  }
  fprintf(fp, "\n");

  pluto_stmts_print(fp, prog->stmts, prog->nstmts);
  pluto_deps_print(fp, prog);
  pluto_transformations_pretty_print(prog);
}

void pluto_dep_free(Dep *dep) {
  pluto_constraints_free(dep->dpolytope);
  pluto_constraints_free(dep->bounding_poly);
  pluto_constraints_free(dep->depsat_poly);
  if (dep->dirvec) {
    free(dep->dirvec);
  }
  if (dep->dirvec) {
    free(dep->satvec);
  }
  pluto_constraints_free(dep->cst);
  pluto_constraints_free(dep->bounding_cst);
  free(dep);
}

/* Temporary data structure used inside extract_deps.
 *
 * deps points to the array of Deps being constructed
 * type is the type of the next Dep
 * index is the index of the next Dep in the array.
 */
struct pluto_extra_dep_info {
  Dep **deps;
  Stmt **stmts;
  PlutoDepType type;
  int index;
  PlutoContext *context;
};

/* Convert an isl_basic_map describing part of a dependence to a Dep.
 * The names of the input and output spaces are of the form S_d or S_d_e
 * with d an integer identifying the statement, e identifying the access
 * (relative to the statement). If it's of the form S_d_e and read/write
 * accesses for the statement are available, source and target accesses
 * are set for the dependence, otherwise not.
 *
 * isl divs are removed; so this is an over-approximation in some cases
 */
static isl_stat basic_map_extract_dep(__isl_take isl_basic_map *bmap,
                                      void *user) {
  Stmt **stmts;
  Dep *dep;
  struct pluto_extra_dep_info *info;

  info = (struct pluto_extra_dep_info *)user;

  stmts = info->stmts;

  bmap = isl_basic_map_remove_divs(bmap);

  dep = info->deps[info->index];

  dep->id = info->index;

  dep->dirvec = NULL;
  dep->type = info->type;
  dep->src = atoi(isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2);

  /* The range space can be wrapped, if we didn't use lastwriter.  For
   * example, [T, N] -> { S_1_w0[t, i] -> [S_0_r2[t', i'] -> a[o2]] : i'
   * = -1 + i and o2 = i and t >= 0 and i >= 3 and i <= -2 + N and t' > t
   * and t' < T } */
  isl_space *space = isl_basic_map_get_space(bmap);
  if (isl_space_range_is_wrapping(space)) {
    bmap = isl_basic_map_range_factor_domain(bmap);
    isl_space_free(space);
    space = isl_basic_map_get_space(bmap);
  }
  dep->dest = atoi(isl_space_get_tuple_name(space, isl_dim_out) + 2);
  dep->dpolytope = isl_basic_map_to_pluto_constraints(bmap, info->context);
  dep->bounding_poly = pluto_constraints_dup(dep->dpolytope);
  isl_space_free(space);

  /* Inconsistent dependence if this assertion fails */
  assert(dep->dpolytope->ncols == stmts[dep->src]->dim + stmts[dep->dest]->dim +
                                      stmts[dep->src]->domain->ncols -
                                      stmts[dep->src]->dim);

  pluto_constraints_set_names_range(dep->dpolytope, stmts[dep->src]->iterators,
                                    0, 0, stmts[dep->src]->dim);

  /* Suffix the destination iterators with a '*/
  char **dnames = (char **)malloc(stmts[dep->dest]->dim * sizeof(char *));
  for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
    dnames[j] = (char *)malloc(strlen(stmts[dep->dest]->iterators[j]) + 2);
    strcpy(dnames[j], stmts[dep->dest]->iterators[j]);
    strcat(dnames[j], "'");
  }
  pluto_constraints_set_names_range(
      dep->dpolytope, dnames, stmts[dep->src]->dim, 0, stmts[dep->dest]->dim);
  for (unsigned j = 0; j < stmts[dep->dest]->dim; j++) {
    free(dnames[j]);
  }
  free(dnames);

  /* parameters */
  pluto_constraints_set_names_range(
      dep->dpolytope, stmts[dep->dest]->domain->names,
      stmts[dep->src]->dim + stmts[dep->dest]->dim, stmts[dep->dest]->dim,
      stmts[dep->dest]->domain->ncols - stmts[dep->dest]->dim - 1);

  if (info->context->options->isldepaccesswise) {
    /* Extract access function information */
    int src_acc_num, dest_acc_num;
    char src_type, dest_type;
    const char *src_name, *dest_name;
    src_name = isl_basic_map_get_tuple_name(bmap, isl_dim_in) + 2;
    while (*src_name != '\0' && *(src_name++) != '_')
      ;
    if (*src_name != '\0') {
      src_type = *src_name;
      src_acc_num = atoi(src_name + 1);
    } else
      assert(0); // access function num not encoded in dependence

    dest_name = isl_basic_map_get_tuple_name(bmap, isl_dim_out) + 2;
    while (*dest_name != '\0' && *(dest_name++) != '_')
      ;
    if (*dest_name != '\0') {
      dest_type = *dest_name;
      dest_acc_num = atoi(dest_name + 1);
    } else {
      assert(0 && "access function num not encoded in dependence");
    }

    switch (info->type) {
    case PLUTO_DEP_RAW:
      dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
      break;
    case PLUTO_DEP_WAW:
      dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
      break;
    /*
     * Sometimes dep_war from isl are not WAR deps; there are WAW deps
     * included in the may deps and we can't assume that they are all WAR may
     * deps. Mark them correctly.
     */
    case PLUTO_DEP_WAR:
      dep->dest_acc = stmts[dep->dest]->writes[dest_acc_num];
      if (src_type == 'w' && dest_type == 'w') {
        /* Fix the type */
        dep->type = PLUTO_DEP_WAW;
        /* This is really a WAW dep */
        dep->src_acc = stmts[dep->src]->writes[src_acc_num];
      } else {
        dep->src_acc = stmts[dep->src]->reads[src_acc_num];
      }
      break;
    case PLUTO_DEP_RAR:
      dep->src_acc = stmts[dep->src]->reads[src_acc_num];
      dep->dest_acc = stmts[dep->dest]->reads[dest_acc_num];
      break;
    default:
      assert(0 && "Unexpected dependence type");
    }
  } else {
    dep->src_acc = NULL;
    dep->dest_acc = NULL;
  }

  info->index++;
  isl_basic_map_free(bmap);
  return isl_stat_ok;
}

/* Extract Pluto dependences from an isl_map */
static isl_stat map_extract_dep(__isl_take isl_map *map, void *user) {
  isl_stat r = isl_map_foreach_basic_map(map, &basic_map_extract_dep, user);
  isl_map_free(map);
  return r;
}

/// Extract a Pluto access function from an access encoded in an isl_basic_map.
static isl_stat
isl_basic_map_extract_access_func(__isl_take isl_basic_map *bmap, void *user) {
  isl_map *map = isl_map_from_basic_map(bmap);

  int dim = isl_map_dim(map, isl_dim_out);
  int ncols =
      isl_map_dim(map, isl_dim_in) + isl_map_dim(map, isl_dim_param) + 1;

  struct pluto_access_meta_info *info = (struct pluto_access_meta_info *)user;

  PlutoMatrix *acc_mat = pluto_matrix_alloc(0, ncols, info->context);

  for (int i = 0; i < dim; i++) {
    PlutoMatrix *func_onedim = NULL;
    if (isl_map_dim_is_single_valued(map, i)) {
      isl_pw_aff *pw_aff = isl_pw_aff_from_map_dim(map, i);
      /* Best effort: Gets it from the last piece */
      struct pluto_mat_context_info m_info = {&func_onedim, info->context};
      isl_pw_aff_foreach_piece(pw_aff, isl_aff_to_pluto_func, &m_info);
      pluto_matrix_add(acc_mat, func_onedim);
      pluto_matrix_free(func_onedim);
      isl_pw_aff_free(pw_aff);
    } else {
      // Multi-valued map for access function? Setting this to zero.
      pluto_matrix_add_row(acc_mat, 0);
      pluto_matrix_zero_row(acc_mat, 0);
    }
  }
  (*info->accs)[info->index] = (PlutoAccess *)malloc(sizeof(PlutoAccess));
  PlutoAccess *acc = (*info->accs)[info->index];
  acc->name = strdup(isl_basic_map_get_tuple_name(bmap, isl_dim_out));
  acc->mat = acc_mat;
  info->index++;
  isl_map_free(map);
  return isl_stat_ok;
}

/// Extract a Pluto access function from an isl_map.
isl_stat isl_map_extract_access_func(__isl_take isl_map *map, void *user) {
  /* Extract a PlutoAccess from every isl_basic_map */
  isl_stat r =
      isl_map_foreach_basic_map(map, &isl_basic_map_extract_access_func, user);
  isl_map_free(map);
  return r;
}

/// Extract read and writes accesses encoded in isl_union_map structures and
/// populate Pluto stmt's accesses. The access relation maps statement domain
/// instances to elements in a data space. The name on the output tuple provides
/// array name.
void extract_accesses_for_pluto_stmt(Stmt *stmt,
                                     __isl_keep isl_union_map *reads,
                                     __isl_keep isl_union_map *writes,
                                     PlutoContext *context) {
  isl_union_map_foreach_map(reads, &isl_map_count, &stmt->nreads);
  isl_union_map_foreach_map(writes, &isl_map_count, &stmt->nwrites);

  int npar = stmt->domain->ncols - stmt->dim - 1;

  struct pluto_access_meta_info e_reads = {&stmt->reads, 0, stmt->dim, npar,
                                           context};
  struct pluto_access_meta_info e_writes = {&stmt->writes, 0, stmt->dim, npar,
                                            context};

  if (stmt->nreads > 0) {
    stmt->reads = (PlutoAccess **)malloc(stmt->nreads * sizeof(PlutoAccess *));
  }
  if (stmt->nwrites > 0) {
    stmt->writes =
        (PlutoAccess **)malloc(stmt->nwrites * sizeof(PlutoAccess *));
  }
  for (int j = 0; j < stmt->nreads; j++) {
    stmt->reads[j] = NULL;
  }
  for (int j = 0; j < stmt->nwrites; j++) {
    stmt->writes[j] = NULL;
  }

  isl_union_map_foreach_map(reads, &isl_map_extract_access_func, &e_reads);
  isl_union_map_foreach_map(writes, &isl_map_extract_access_func, &e_writes);

  // Verify correct extraction.
  for (int j = 0; j < stmt->nreads; j++) {
    assert(stmt->reads[j]->mat->ncols == stmt->dim + npar + 1 &&
           "Inconsistent access matrix");
  }
  for (int j = 0; j < stmt->nwrites; j++) {
    // Verify correct extraction.
    assert(stmt->writes[j]->mat->ncols == stmt->dim + npar + 1 &&
           "Inconsistent access matrix");
  }
}

/// Extract deps from isl union map into Pluto Deps.
int extract_deps_from_isl_union_map(__isl_keep isl_union_map *umap, Dep **deps,
                                    int first, Stmt **stmts, PlutoDepType type,
                                    PlutoContext *context) {
  struct pluto_extra_dep_info info = {deps, stmts, type, first, context};

  isl_union_map_foreach_map(umap, &map_extract_dep, &info);

  return info.index - first;
}

// Compute dependences using ISL.
// If options->lastwriter is false, then
//       RAW deps are those from any earlier write to a read
//       WAW deps are those from any earlier write to a write
//       WAR deps are those from any earlier read to a write
//       RAR deps are those from any earlier read to a read
//  If options->lastwriter is true, then
//       RAW deps are those from the last write to a read
//       WAW deps are those from the last write to a write
//       WAR deps are those from any earlier read not masked by an intermediate
//       write to a write
//       RAR deps are those from the last read to a read
//
//  The RAR deps are only computed if options->rar is set.
void compute_deps_isl(__isl_keep isl_union_map *reads,
                      __isl_keep isl_union_map *writes,
                      __isl_keep isl_union_map *schedule,
                      __isl_keep isl_union_map *empty, isl_union_map **dep_raw,
                      isl_union_map **dep_war, isl_union_map **dep_waw,
                      isl_union_map **dep_rar, isl_union_map **trans_dep_war,
                      isl_union_map **trans_dep_waw, PlutoOptions *options) {
  assert(options && "options not set");

  if (options->lastwriter) {
    // Compute RAW dependences with last writer (no transitive dependences).
    isl_union_map_compute_flow(
        isl_union_map_copy(reads), isl_union_map_copy(writes),
        isl_union_map_copy(empty), isl_union_map_copy(schedule), dep_raw, NULL,
        NULL, NULL);
    // Compute WAW and WAR dependences without transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(writes),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), dep_waw,
        dep_war, NULL, NULL);
    // Compute WAR dependences with transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL,
        trans_dep_war, NULL, NULL);
    // Compute WAW dependences with transitive dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL,
        trans_dep_waw, NULL, NULL);
    if (options->rar) {
      // Compute RAR dependences without transitive dependences.
      isl_union_map_compute_flow(
          isl_union_map_copy(reads), isl_union_map_copy(reads),
          isl_union_map_copy(empty), isl_union_map_copy(schedule), dep_rar,
          NULL, NULL, NULL);
    }
  } else {
    // Without lastwriter, compute transitive dependences.
    // RAW dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(reads), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL, dep_raw,
        NULL, NULL);
    // WAR dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL, dep_war,
        NULL, NULL);
    // WAW dependences.
    isl_union_map_compute_flow(
        isl_union_map_copy(writes), isl_union_map_copy(empty),
        isl_union_map_copy(writes), isl_union_map_copy(schedule), NULL, dep_waw,
        NULL, NULL);
    if (options->rar) {
      // RAR dependences.
      isl_union_map_compute_flow(
          isl_union_map_copy(reads), isl_union_map_copy(empty),
          isl_union_map_copy(reads), isl_union_map_copy(schedule), NULL,
          dep_rar, NULL, NULL);
    }
  }

  if (options->isldepcoalesce) {
    *dep_raw = isl_union_map_coalesce(*dep_raw);
    *dep_war = isl_union_map_coalesce(*dep_war);
    *dep_waw = isl_union_map_coalesce(*dep_waw);

    if (options->lastwriter) {
      *trans_dep_war = isl_union_map_coalesce(*trans_dep_war);
      *trans_dep_waw = isl_union_map_coalesce(*trans_dep_waw);
    }
  }
}

int read_codegen_context_from_file(PlutoConstraints *codegen_context) {
  FILE *fp = fopen("codegen.context", "r");
  PlutoContext *context = codegen_context->context;

  if (fp) {
    IF_DEBUG(printf("[Pluto] Reading from codegen.context\n"););
    PlutoConstraints *cc = pluto_constraints_read(fp, context);
    if (cc && cc->ncols == codegen_context->ncols) {
      pluto_constraints_add(codegen_context, cc);
      return 0;
    }
    IF_DEBUG(printf("[WARNING] Failed to read from codegen.context\n"););
  }

  return 1;
}

/* Get an upper bound for transformation coefficients to prevent spurious
 * transformations that represent shifts or skews proportional to trip counts:
 * this happens when loop bounds are constants
 */
int pluto_prog_get_largest_const_in_domains(const PlutoProg *prog) {
  int max = 0;
  for (unsigned i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    for (unsigned r = 0; r < stmt->domain->nrows; r++) {
      max = PLMAX(max, stmt->domain->val[r][stmt->domain->ncols - 1]);
    }
  }

  return max - 1;
}

PlutoProg *pluto_prog_alloc(PlutoContext *context) {
  assert(context->options && "options null");

  PlutoProg *prog = (PlutoProg *)malloc(sizeof(PlutoProg));

  prog->context = context;
  prog->nstmts = 0;
  prog->stmts = NULL;
  prog->npar = 0;
  prog->nvar = 0;
  prog->params = NULL;
  prog->param_context = pluto_constraints_alloc(1, prog->npar + 1, context);
  prog->codegen_context = pluto_constraints_alloc(1, prog->npar + 1, context);
  prog->deps = NULL;
  prog->ndeps = 0;
  prog->transdeps = NULL;
  prog->ntransdeps = 0;
  prog->ddg = NULL;
  prog->fcg = NULL;
  prog->hProps = NULL;
  prog->num_hyperplanes = 0;
  prog->decls = (char *)malloc(16384 * 9);
  prog->data_names = NULL;
  prog->num_data = 0;

  strcpy(prog->decls, "");

  prog->globcst = NULL;

  prog->num_parameterized_loops = -1;

  return prog;
}

void pluto_prog_free(PlutoProg *prog) {
  /* Free dependences */
  for (int i = 0; i < prog->ndeps; i++) {
    pluto_dep_free(prog->deps[i]);
  }
  free(prog->deps);

  for (int i = 0; i < prog->ntransdeps; i++) {
    pluto_dep_free(prog->transdeps[i]);
  }
  free(prog->transdeps);

  /* Free DDG */
  if (prog->ddg != NULL) {
    graph_free(prog->ddg);
  }

  if (prog->hProps != NULL) {
    free(prog->hProps);
  }

  for (int i = 0; i < prog->npar; i++) {
    free(prog->params[i]);
  }
  if (prog->npar >= 1) {
    free(prog->params);
  }

  /* Statements */
  for (unsigned i = 0; i < prog->nstmts; i++) {
    pluto_stmt_free(prog->stmts[i]);
  }
  if (prog->nstmts >= 1) {
    free(prog->stmts);
  }

  pluto_constraints_free(prog->param_context);
  pluto_constraints_free(prog->codegen_context);

  pluto_constraints_free(prog->globcst);

  free(prog->decls);

  for (int i = 0; i < prog->num_data; i++) {
    free(prog->data_names[i]);
  }
  free(prog->data_names);

  free(prog);
}

PlutoOptions *pluto_options_alloc() {
  PlutoOptions *options;

  options = (PlutoOptions *)malloc(sizeof(PlutoOptions));

  /* Initialize to default */
  options->flic = 0;
  options->tile = 1;
  options->intratileopt = 1;
  options->dynschedule = 0;
  options->dynschedule_graph = 0;
  options->dynschedule_graph_old = 0;
  options->dyn_trans_deps_tasks = 0;
  options->debug = 0;
  options->moredebug = 0;
  options->scancount = 0;
  options->parallel = 1;
  options->innerpar = 0;
  options->identity = 0;

  options->pet = 0;

  /* Enable one dimension of concurrent startup by default */
  options->diamondtile = 1;
  options->fulldiamondtile = 0;

  options->per_cc_obj = 0;

  options->iss = 0;
  options->unrolljam = 1;

  /* Unroll/jam factor */
  options->ufactor = 8;

  /* Ignore input deps */
  options->rar = 0;

  /* Override for first and last levels to tile */
  options->ft = -1;
  options->lt = -1;

  /* Override for first and last cloog options */
  options->cloogf = -1;
  options->cloogl = -1;

  options->cloogsh = 0;

  options->cloogbacktrack = 1;

  options->multipar = 0;
  options->second_level_tile = 0;
  options->prevector = 1;
  options->fuse = kSmartFuse;

  /* Experimental */
  options->delayed_cut = 0;
  options->hybridcut = 0;

  /* Default context is no context */
  options->codegen_context = -1;

  options->coeff_bound = -1;

  options->forceparallel = 0;

  options->bee = 0;

  options->isldep = 0;
  options->isldepaccesswise = 1;
  /* Disabled due to a potential bug in coalescing. Reproduce with
   * examples/heat-2d/heat-2d.c - coalescing dep_raw leads to no hyperplanes
   * being found. */
  options->isldepcoalesce = 0;

  options->candldep = 0;

  options->pipsolve = 0;
  options->islsolve = 1;
  options->glpk = 0;
  options->gurobi = 0;

  options->lp = 0;
  options->dfp = 0;
  options->ilp = 0;

  options->lpcolour = 0;
  options->scc_cluster = 0;

  options->readscop = 0;

  options->lastwriter = 0;

  options->nodepbound = 0;

  options->scalpriv = 0;

  options->silent = 0;

  options->out_file = NULL;

  options->time = 1;

  return options;
}

PlutoContext *pluto_context_alloc() {
  PlutoContext *context = (PlutoContext *)malloc(sizeof(PlutoContext));
  context->options = pluto_options_alloc();
  return context;
}

void pluto_context_free(PlutoContext *context) {
  pluto_options_free(context->options);
  free(context);
}

/* Add global/program parameter at position 'pos' */
void pluto_prog_add_param(PlutoProg *prog, const char *param, int pos) {
  for (unsigned i = 0; i < prog->nstmts; i++) {
    Stmt *stmt = prog->stmts[i];
    pluto_constraints_add_dim(
        stmt->domain, stmt->domain->ncols - 1 - prog->npar + pos, param);
    pluto_matrix_add_col(stmt->trans,
                         stmt->trans->ncols - 1 - prog->npar + pos);

    for (int j = 0; j < stmt->nwrites; j++) {
      pluto_matrix_add_col(stmt->writes[j]->mat, stmt->dim + pos);
    }
    for (int j = 0; j < stmt->nreads; j++) {
      pluto_matrix_add_col(stmt->reads[j]->mat, stmt->dim + pos);
    }
  }
  for (int i = 0; i < prog->ndeps; i++) {
    pluto_constraints_add_dim(
        prog->deps[i]->dpolytope,
        prog->deps[i]->dpolytope->ncols - 1 - prog->npar + pos, NULL);
  }
  pluto_constraints_add_dim(prog->param_context,
                            prog->param_context->ncols - 1 - prog->npar + pos,
                            param);
  pluto_constraints_add_dim(prog->codegen_context,
                            prog->codegen_context->ncols - 1 - prog->npar + pos,
                            param);

  prog->params =
      (char **)realloc(prog->params, sizeof(char *) * (prog->npar + 1));

  for (int i = prog->npar - 1; i >= pos; i--) {
    prog->params[i + 1] = prog->params[i];
  }

  prog->params[pos] = strdup(param);
  prog->npar++;
}

void pluto_options_free(PlutoOptions *options) {
  if (options->out_file != NULL) {
    free(options->out_file);
  }
  free(options);
}

/* pos: position of domain iterator
 * time_pos: position of time iterator; iter: domain iterator; supply -1
 * if you don't want a scattering function row added for it */
void pluto_stmt_add_dim(Stmt *stmt, unsigned pos, int time_pos,
                        const char *iter, PlutoHypType hyp_type,
                        PlutoProg *prog) {
  int i, npar;

  npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(pos <= stmt->dim);
  assert(time_pos <= (int)stmt->trans->nrows);
  assert(stmt->dim + npar + 1 == stmt->domain->ncols);

  pluto_constraints_add_dim(stmt->domain, pos, NULL);
  stmt->dim++;
  stmt->iterators =
      (char **)realloc(stmt->iterators, stmt->dim * sizeof(char *));
  for (i = stmt->dim - 2; i >= (int)pos; i--) {
    stmt->iterators[i + 1] = stmt->iterators[i];
  }
  stmt->iterators[pos] = strdup(iter);

  /* Stmt should always have a transformation */
  assert(stmt->trans != NULL);
  pluto_matrix_add_col(stmt->trans, pos);

  if (time_pos != -1) {
    pluto_matrix_add_row(stmt->trans, time_pos);
    stmt->trans->val[time_pos][pos] = 1;

    stmt->hyp_types = (PlutoHypType *)realloc(
        stmt->hyp_types, sizeof(PlutoHypType) * stmt->trans->nrows);
    for (i = stmt->trans->nrows - 2; i >= time_pos; i--) {
      stmt->hyp_types[i + 1] = stmt->hyp_types[i];
    }
    stmt->hyp_types[time_pos] = hyp_type;
  }

  /* Update is_orig_loop */
  stmt->is_orig_loop =
      (bool *)realloc(stmt->is_orig_loop, sizeof(bool) * stmt->dim);
  for (i = stmt->dim - 2; i >= (int)pos; i--) {
    stmt->is_orig_loop[i + 1] = stmt->is_orig_loop[i];
  }
  stmt->is_orig_loop[pos] = true;

  for (i = 0; i < stmt->nwrites; i++) {
    pluto_matrix_add_col(stmt->writes[i]->mat, pos);
  }
  for (i = 0; i < stmt->nreads; i++) {
    pluto_matrix_add_col(stmt->reads[i]->mat, pos);
  }

  /* Update dependences */
  for (i = 0; i < prog->ndeps; i++) {
    if (prog->deps[i]->src == stmt->id) {
      pluto_constraints_add_dim(prog->deps[i]->dpolytope, pos, NULL);
      pluto_constraints_add_dim(prog->deps[i]->bounding_poly, pos, NULL);
    }
    if (prog->deps[i]->dest == stmt->id) {
      pluto_constraints_add_dim(prog->deps[i]->dpolytope,
                                prog->stmts[prog->deps[i]->src]->dim + pos,
                                NULL);
      pluto_constraints_add_dim(prog->deps[i]->bounding_poly,
                                prog->stmts[prog->deps[i]->src]->dim + pos,
                                NULL);
    }
  }

  for (i = 0; i < prog->ntransdeps; i++) {
    assert(prog->transdeps[i] != NULL);
    if (prog->transdeps[i]->src == stmt->id) {
      pluto_constraints_add_dim(prog->transdeps[i]->dpolytope, pos, NULL);
      pluto_constraints_add_dim(prog->transdeps[i]->bounding_poly, pos, NULL);
    }
    if (prog->transdeps[i]->dest == stmt->id) {
      pluto_constraints_add_dim(prog->transdeps[i]->dpolytope,
                                prog->stmts[prog->transdeps[i]->src]->dim + pos,
                                NULL);
      pluto_constraints_add_dim(prog->transdeps[i]->bounding_poly,
                                prog->stmts[prog->transdeps[i]->src]->dim + pos,
                                NULL);
    }
  }
}

/* Warning: use it only to knock off a dummy dimension (unrelated to
 * anything else */
void pluto_stmt_remove_dim(Stmt *stmt, unsigned pos, PlutoProg *prog) {
  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(pos <= stmt->dim);
  assert(stmt->dim + npar + 1 == stmt->domain->ncols);

  pluto_constraints_remove_dim(stmt->domain, pos);
  stmt->dim--;

  if (stmt->iterators != NULL) {
    free(stmt->iterators[pos]);
    for (int i = pos; i <= (int)stmt->dim - 1; i++) {
      stmt->iterators[i] = stmt->iterators[i + 1];
    }
    stmt->iterators =
        (char **)realloc(stmt->iterators, stmt->dim * sizeof(char *));
  }

  pluto_matrix_remove_col(stmt->trans, pos);

  /* Update is_orig_loop */
  for (int i = pos; i <= (int)stmt->dim - 1; i++) {
    stmt->is_orig_loop[i] = stmt->is_orig_loop[i + 1];
  }
  stmt->is_orig_loop =
      (bool *)realloc(stmt->is_orig_loop, sizeof(bool) * stmt->dim);

  for (int i = 0; i < stmt->nwrites; i++) {
    pluto_matrix_remove_col(stmt->writes[i]->mat, pos);
  }

  for (int i = 0; i < stmt->nreads; i++) {
    pluto_matrix_remove_col(stmt->reads[i]->mat, pos);
  }

  /* Update deps */
  for (int i = 0; i < prog->ndeps; i++) {
    if (prog->deps[i]->src == stmt->id) {
      pluto_constraints_remove_dim(prog->deps[i]->dpolytope, pos);
    }
    if (prog->deps[i]->dest == stmt->id) {
      pluto_constraints_remove_dim(prog->deps[i]->dpolytope,
                                   prog->stmts[prog->deps[i]->src]->dim + pos);
    }
  }

  for (int i = 0; i < prog->ntransdeps; i++) {
    assert(prog->transdeps[i] != NULL);
    if (prog->transdeps[i]->src == stmt->id) {
      pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope, pos);
    }
    if (prog->transdeps[i]->dest == stmt->id) {
      pluto_constraints_remove_dim(prog->transdeps[i]->dpolytope,
                                   prog->stmts[prog->transdeps[i]->src]->dim +
                                       pos);
    }
  }
}

void pluto_stmt_add_hyperplane(Stmt *stmt, PlutoHypType type, unsigned pos) {
  assert(pos <= stmt->trans->nrows);

  pluto_matrix_add_row(stmt->trans, pos);

  stmt->hyp_types = (PlutoHypType *)realloc(
      stmt->hyp_types, sizeof(PlutoHypType) * stmt->trans->nrows);
  for (int i = stmt->trans->nrows - 2; i >= (int)pos; i--) {
    stmt->hyp_types[i + 1] = stmt->hyp_types[i];
  }
  stmt->hyp_types[pos] = type;

  if (stmt->first_tile_dim >= (int)pos)
    stmt->first_tile_dim++;
  if (stmt->last_tile_dim >= (int)pos)
    stmt->last_tile_dim++;
}

void pluto_prog_add_hyperplane(PlutoProg *prog, int pos,
                               PlutoHypType hyp_type) {
  int i;

  prog->num_hyperplanes++;
  prog->hProps = (HyperplaneProperties *)realloc(
      prog->hProps, prog->num_hyperplanes * sizeof(HyperplaneProperties));

  for (i = prog->num_hyperplanes - 2; i >= pos; i--) {
    prog->hProps[i + 1] = prog->hProps[i];
  }
  /* Initialize some */
  prog->hProps[pos].unroll = NO_UNROLL;
  prog->hProps[pos].prevec = 0;
  prog->hProps[pos].band_num = -1;
  prog->hProps[pos].dep_prop = UNKNOWN;
  prog->hProps[pos].type = hyp_type;
}

/*
 * Create a new statement (see also pluto_stmt_dup)
 */
Stmt *pluto_create_stmt(int dim, const PlutoConstraints *domain,
                        const PlutoMatrix *trans, char **iterators,
                        const char *text, PlutoStmtType type) {
  Stmt *stmt = pluto_stmt_alloc(dim, domain, trans);

  stmt->type = type;

  stmt->text = strdup(text);

  for (unsigned i = 0; i < stmt->dim; i++) {
    stmt->iterators[i] = strdup(iterators[i]);
  }

  pluto_constraints_set_names_range(stmt->domain, stmt->iterators, 0, 0,
                                    stmt->dim);

  /* TODO: Set names for parameters */

  return stmt;
}

/* Pad statement transformations so that they all equal number
 * of rows */
void pluto_pad_stmt_transformations(PlutoProg *prog) {
  int i, nstmts;

  nstmts = prog->nstmts;
  Stmt **stmts = prog->stmts;

  /* Pad all trans if necessary with zeros */
  unsigned max_nrows = 0;
  for (i = 0; i < nstmts; i++) {
    if (stmts[i]->trans != NULL) {
      max_nrows = PLMAX(max_nrows, stmts[i]->trans->nrows);
    }
  }

  if (max_nrows >= 1) {
    for (i = 0; i < nstmts; i++) {
      if (stmts[i]->trans == NULL) {
        stmts[i]->trans = pluto_matrix_alloc(
            max_nrows, stmts[i]->dim + prog->npar + 1, prog->context);
        stmts[i]->trans->nrows = 0;
      }

      unsigned curr_rows = stmts[i]->trans->nrows;

      /* Add all zero rows */
      for (unsigned j = curr_rows; j < max_nrows; j++) {
        pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
      }
    }

    unsigned old_hyp_num = prog->num_hyperplanes;
    for (unsigned i = old_hyp_num; i < max_nrows; i++) {
      /* This is not really H_SCALAR, but this is the best we can do */
      pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);
    }
  }
}

/* Add statement to program; can't reuse arg stmt pointer any more */
void pluto_add_given_stmt(PlutoProg *prog, Stmt *stmt) {
  prog->stmts =
      (Stmt **)realloc(prog->stmts, ((prog->nstmts + 1) * sizeof(Stmt *)));

  stmt->id = prog->nstmts;

  prog->nvar = PLMAX(prog->nvar, (int)stmt->dim);
  prog->stmts[prog->nstmts] = stmt;
  prog->nstmts++;

  pluto_pad_stmt_transformations(prog);
}

/* Create a statement and add it to the program
 * iterators: domain iterators
 * trans: schedule/transformation
 * domain: domain
 * text: statement text
 */
void pluto_add_stmt(PlutoProg *prog, const PlutoConstraints *domain,
                    const PlutoMatrix *trans, char **iterators,
                    const char *text, PlutoStmtType type) {
  int nstmts;

  assert(trans != NULL);
  assert(trans->ncols == domain->ncols);

  nstmts = prog->nstmts;

  prog->stmts = (Stmt **)realloc(prog->stmts, ((nstmts + 1) * sizeof(Stmt *)));

  Stmt *stmt = pluto_create_stmt(domain->ncols - prog->npar - 1, domain, trans,
                                 iterators, text, type);
  stmt->id = nstmts;

  /* Initialize intra statement deps to Null. Will be updated when fcg is built
   */
  stmt->intra_stmt_dep_cst = NULL;

  prog->nvar = PLMAX(prog->nvar, (int)stmt->dim);

  prog->stmts[nstmts] = stmt;
  prog->nstmts++;

  pluto_pad_stmt_transformations(prog);
}

Dep *pluto_dep_alloc() {
  Dep *dep = (Dep *)malloc(sizeof(Dep));

  dep->id = -1;
  dep->satvec = NULL;
  dep->dpolytope = NULL;
  dep->bounding_poly = NULL;
  dep->depsat_poly = NULL;
  dep->satisfied = false;
  dep->satisfaction_level = -1;
  dep->dirvec = NULL;
  dep->src_acc = NULL;
  dep->dest_acc = NULL;
  dep->cst = NULL;
  dep->bounding_cst = NULL;
  dep->src_unique_dpolytope = NULL;

  return dep;
}

Dep *pluto_dep_dup(Dep *d) {
  Dep *dep = (Dep *)malloc(sizeof(Dep));

  dep->id = d->id;
  dep->src = d->src;
  dep->dest = d->dest;
  dep->src_acc = d->src_acc;
  dep->dest_acc = d->dest_acc;
  dep->dpolytope = pluto_constraints_dup(d->dpolytope);
  dep->bounding_poly = pluto_constraints_dup(d->bounding_poly);

  dep->src_unique_dpolytope =
      d->src_unique_dpolytope ? pluto_constraints_dup(d->src_unique_dpolytope)
                              : NULL;

  dep->depsat_poly =
      d->depsat_poly ? pluto_constraints_dup(d->depsat_poly) : NULL;
  dep->satvec = NULL; // TODO
  dep->type = d->type;
  dep->satisfied = d->satisfied;
  dep->satisfaction_level = d->satisfaction_level;
  dep->dirvec = NULL; // TODO
  dep->cst = d->cst ? pluto_constraints_dup(d->cst) : NULL;
  dep->bounding_cst =
      d->bounding_cst ? pluto_constraints_dup(d->bounding_cst) : NULL;

  return dep;
}

/*
 * Only very essential information is needed to allocate; rest can be
 * populated as needed
 */
Stmt *pluto_stmt_alloc(unsigned dim, const PlutoConstraints *domain,
                       const PlutoMatrix *trans) {
  /* Have to provide a transformation */
  assert(trans != NULL);

  Stmt *stmt = (Stmt *)malloc(sizeof(Stmt));

  /* id will be assigned when added to PlutoProg */
  stmt->id = -1;
  stmt->dim = dim;
  stmt->dim_orig = dim;
  if (domain != NULL) {
    stmt->domain = pluto_constraints_dup(domain);
  } else {
    stmt->domain = NULL;
  }

  stmt->trans = pluto_matrix_dup(trans);

  stmt->hyp_types =
      (PlutoHypType *)malloc(stmt->trans->nrows * sizeof(PlutoHypType));
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    stmt->hyp_types[i] = H_LOOP;
  }

  stmt->text = NULL;
  stmt->tile = 1;
  stmt->num_tiled_loops = 0;
  stmt->reads = NULL;
  stmt->writes = NULL;
  stmt->nreads = 0;
  stmt->nwrites = 0;

  /* For diamond tiling */
  stmt->evicted_hyp = NULL;
  stmt->evicted_hyp_pos = -1;

  stmt->first_tile_dim = 0;
  stmt->last_tile_dim = -1;

  stmt->type = STMT_UNKNOWN;
  stmt->ploop_id = -1;

  if (dim >= 1) {
    stmt->is_orig_loop = (bool *)malloc(dim * sizeof(bool));
    stmt->iterators = (char **)malloc(sizeof(char *) * dim);
    for (unsigned i = 0; i < stmt->dim; i++) {
      stmt->iterators[i] = NULL;
    }
  } else {
    stmt->is_orig_loop = NULL;
    stmt->iterators = NULL;
  }

  return stmt;
}

PlutoAccess *pluto_access_dup(const PlutoAccess *acc) {
  assert(acc);

  PlutoAccess *nacc = (PlutoAccess *)malloc(sizeof(PlutoAccess));
  nacc->mat = pluto_matrix_dup(acc->mat);
  nacc->name = strdup(acc->name);
  nacc->sym_id = acc->sym_id;

  return nacc;
}

void pluto_access_free(PlutoAccess *acc) {
  if (acc) {
    pluto_matrix_free(acc->mat);
    free(acc->name);
    free(acc);
  }
}

void pluto_stmt_free(Stmt *stmt) {
  pluto_constraints_free(stmt->domain);

  pluto_matrix_free(stmt->trans);

  free(stmt->hyp_types);
  free(stmt->text);

  for (unsigned j = 0; j < stmt->dim; j++) {
    if (stmt->iterators[j] != NULL) {
      free(stmt->iterators[j]);
    }
  }

  /* If dim is zero, iterators, is_orig_loop are NULL */
  if (stmt->iterators != NULL) {
    free(stmt->iterators);
    free(stmt->is_orig_loop);
  }

  PlutoAccess **writes = stmt->writes;
  PlutoAccess **reads = stmt->reads;

  if (writes != NULL) {
    for (int i = 0; i < stmt->nwrites; i++) {
      pluto_access_free(writes[i]);
    }
    free(writes);
  }
  if (reads != NULL) {
    for (int i = 0; i < stmt->nreads; i++) {
      pluto_access_free(reads[i]);
    }
    free(reads);
  }

  pluto_matrix_free(stmt->evicted_hyp);

  free(stmt);
}

/* Get transformed domain */
PlutoConstraints *pluto_get_new_domain(const Stmt *stmt) {
  PlutoConstraints *sched;

  PlutoConstraints *newdom = pluto_constraints_dup(stmt->domain);
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    pluto_constraints_add_dim(newdom, 0, NULL);
  }

  sched = pluto_stmt_get_schedule(stmt);

  pluto_constraints_intersect(newdom, sched);

  pluto_constraints_project_out(newdom, stmt->trans->nrows, stmt->dim);

  pluto_constraints_free(sched);

  return newdom;
}

/*
 * Checks if the range of the variable at depth 'depth' can be bound by a
 * constant; returns the constant of -1 if it can't be
 *
 * WARNING: If cnst is a list, looks at just the first element
 *
 * TODO: Not general now: difference being constant can be implied through
 * other inequalities
 *
 * */
int get_const_bound_difference(const PlutoConstraints *cnst, int depth) {
  int constdiff, c, _lcm;

  assert(cnst != NULL);
  PlutoConstraints *cst = pluto_constraints_dup(cnst);

  pluto_constraints_project_out(cst, depth + 1, cst->ncols - 1 - depth - 1);
  assert(depth >= 0 && depth <= (int)cst->ncols - 2);

  constdiff = INT_MAX;

  unsigned r;
  for (r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0)
      break;
  }
  /* Variable doesn't appear */
  if (r == cst->nrows)
    return -1;

  /* Scale rows so that the coefficient of depth var is the same */
  _lcm = 1;
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0)
      _lcm = lcm(_lcm, llabs(cst->val[r][depth]));
  }
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->val[r][depth] != 0) {
      for (unsigned c = 0; c < cst->ncols; c++) {
        cst->val[r][c] = cst->val[r][c] * (_lcm / llabs(cst->val[r][depth]));
      }
    }
  }

  /* Equality to a function of parameters/constant implies single point */
  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r] && cst->val[r][depth] != 0) {
      for (c = depth + 1; c < (int)cst->ncols - 1; c++) {
        if (cst->val[r][c] != 0) {
          break;
        }
      }
      if (c == (int)cst->ncols - 1) {
        constdiff = 1;
        // printf("constdiff is 1\n");
      }
    }
  }

  for (unsigned r = 0; r < cst->nrows; r++) {
    if (cst->is_eq[r])
      continue;
    if (cst->val[r][depth] <= -1) {
      /* Find a lower bound with constant difference */
      for (unsigned r1 = 0; r1 < cst->nrows; r1++) {
        if (cst->is_eq[r1])
          continue;
        if (cst->val[r1][depth] >= 1) {
          for (c = 0; c < (int)cst->ncols - 1; c++) {
            if (cst->val[r1][c] + cst->val[r][c] != 0) {
              break;
            }
          }
          if (c == (int)cst->ncols - 1) {
            constdiff = PLMIN(
                constdiff,
                floorf(cst->val[r][c] / (float)-cst->val[r][depth]) +
                    ceilf(cst->val[r1][c] / (float)cst->val[r1][depth]) + 1);
          }
        }
      }
    }
  }
  pluto_constraints_free(cst);

  if (constdiff == INT_MAX) {
    return -1;
  }

  /* Sometimes empty sets imply negative difference */
  /* It basically means zero points */
  if (constdiff <= -1)
    constdiff = 0;

  return constdiff;
}

#define MINF 0
#define MAXF 1

/* Get expression for pos^{th} constraint in cst;
 * Returned string should be freed with 'free' */
char *get_expr(PlutoConstraints *cst, int pos, const char **params,
               int bound_type) {
  int c, sum;

  char *expr = (char *)malloc(512);
  strcpy(expr, "");

  // printf("Get expr\n");
  // pluto_constraints_print(stdout, cst);

  if (bound_type == MINF)
    assert(cst->val[pos][0] <= -1);
  else
    assert(cst->val[pos][0] >= 1);

  sum = 0;
  for (c = 1; c < (int)cst->ncols - 1; c++) {
    sum += llabs(cst->val[pos][c]);
  }

  if (sum == 0) {
    /* constant */
    if (bound_type == MINF) {
      sprintf(expr + strlen(expr), "%d",
              (int)floorf(cst->val[pos][cst->ncols - 1] /
                          -(float)cst->val[pos][0]));
    } else {
      sprintf(
          expr + strlen(expr), "%d",
          (int)ceilf(-cst->val[pos][cst->ncols - 1] / (float)cst->val[pos][0]));
    }
  } else {
    /* if it's being divided by 1, make it better by not putting
     * floor/ceil */
    if (llabs(cst->val[pos][0]) != 1) {
      if (bound_type == MINF) {
        sprintf(expr + strlen(expr), "floorf((");
      } else {
        sprintf(expr + strlen(expr), "ceilf((");
      }
    }

    for (c = 1; c < (int)cst->ncols - 1; c++) {
      if (cst->val[pos][c] != 0) {
        if (bound_type == MINF) {
          sprintf(expr + strlen(expr),
                  (cst->val[pos][c] >= 1) ? "+%ld*%s" : "%ld*%s",
                  cst->val[pos][c], params[c - 1]);
        } else {
          sprintf(expr + strlen(expr),
                  (cst->val[pos][c] <= -1) ? "+%ld*%s" : "%ld*%s",
                  -cst->val[pos][c], params[c - 1]);
        }
      }
    }

    if (cst->val[pos][c] != 0) {
      if (bound_type == MINF) {
        sprintf(expr + strlen(expr), (cst->val[pos][c] >= 1) ? "+%ld" : "%ld",
                cst->val[pos][c]);
      } else {
        sprintf(expr + strlen(expr), (cst->val[pos][c] <= -1) ? "+%ld" : "%ld",
                -cst->val[pos][c]);
      }
    }

    /* If it's being divided by 1, make it better by not putting
     * floor/ceil. */
    if (llabs(cst->val[pos][0]) != 1) {
      sprintf(expr + strlen(expr), ")/(float)%ld)",
              (bound_type == MINF) ? -cst->val[pos][0] : cst->val[pos][0]);
    }
  }

  return expr;
}

/*
 * Get min or max of all upper or lower bounds (resp).
 * Returned string should be freed with free
 */
char *get_func_of_expr(PlutoConstraints *cst, int offset, int bound_type,
                       const char **params) {
  char *expr, *expr1;

  char *fexpr = (char *)malloc(512);

  strcpy(fexpr, "");

  char func[5];
  if (bound_type == MINF)
    strcpy(func, "min(");
  else
    strcpy(func, "max(");

  if (cst->nrows - offset == 1) {
    expr = get_expr(cst, offset, params, bound_type);
    strcat(fexpr, expr);
  } else {
    /* cst->nrows >= 2 */
    expr = get_expr(cst, offset, params, bound_type);
    strcat(fexpr, func);
    strcat(fexpr, expr);
    expr1 = get_func_of_expr(cst, offset + 1, bound_type, params);
    strcat(fexpr, ",");
    strcat(fexpr, expr1);
    strcat(fexpr, ")");
    free(expr1);
  }
  free(expr);

  return fexpr;
}

/* Return the size of the parametric bounding box for a (contiguous)
 * block of dimensions
 * start: position of start of block
 * num: number of dimensions in block
 * npar: number of parameters in terms of which expression will be computed;
 * these are assumed to be the last 'npar' variables of cst
 * parmas: strings for 'npar' parameters
 * Return: expression describing the maximum number of points 'block'
 * vars traverse for any value of '0..start-1' variables
 *
 * This function is constant-aware, i.e., if possible, it will exploit the
 * fact that the range of a variable is bounded by a constant. The underlying
 * call to get_parametric_extent_const for each of the 'num' dimensions
 * achieves this.
 */
char *get_parametric_bounding_box(const PlutoConstraints *cst, int start,
                                  int num, int npar, const char **params) {
  int k;
  char *buf_size;

  buf_size = (char *)malloc(2048 * 8);
  strcpy(buf_size, "(");

  const PlutoConstraints *cst_tmp = cst;
  while (cst_tmp != NULL) {
    sprintf(buf_size + strlen(buf_size), "+1");
    for (k = 0; k < num; k++) {
      char *extent;
      get_parametric_extent_const(cst_tmp, start + k, npar, params, &extent,
                                  NULL);
      sprintf(buf_size + strlen(buf_size), "*(%s)", extent);
      free(extent);
    }
    cst_tmp = cst_tmp->next;
  }
  sprintf(buf_size + strlen(buf_size), ")");

  return buf_size;
}

/*  Parametric extent of the pos^th variable in cst
 *  Extent computation is constant-aware, i.e., look when it can be
 *  bounded by a constant; if not, just a difference of max and min
 *  expressions of parameters is returned;  last 'npar'  ones are
 *  treated as parameters; *extent should be freed by 'free'
 */
void get_parametric_extent_const(const PlutoConstraints *cst, int pos, int npar,
                                 const char **params, char **extent,
                                 char **p_lbexpr) {
  int constdiff = get_const_bound_difference(cst, pos);

  if ((p_lbexpr == NULL) && (constdiff != -1)) {
    *extent = (char *)malloc(sizeof(int) * 8);
    sprintf(*extent, "%d", constdiff);
  } else {
    get_parametric_extent(cst, pos, npar, params, extent, p_lbexpr);
  }
}

/* Get lower and upper bound expression as a function of parameters for pos^th
 * variable; last npar in cst are treated as parameters
 * lbexpr and ubexpr should be freed with free
 * */
void get_lb_ub_expr(const PlutoConstraints *cst, int pos, int npar,
                    const char **params, char **lbexpr, char **ubexpr) {
  PlutoConstraints *lb, *ub, *lbs, *ubs;
  char *lbe, *ube;

  PlutoConstraints *dup = pluto_constraints_dup(cst);

  pluto_constraints_project_out(dup, 0, pos);
  pluto_constraints_project_out(dup, 1, dup->ncols - npar - 1 - 1);

  lbs = pluto_constraints_alloc(dup->nrows, dup->ncols, cst->context);
  ubs = pluto_constraints_alloc(dup->nrows, dup->ncols, cst->context);

  for (unsigned i = 0; i < dup->nrows; i++) {
    if (dup->is_eq[i] && dup->val[i][0] != 0) {
      lb = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(lbs, lb);
      pluto_constraints_free(lb);

      ub = pluto_constraints_select_row(dup, i);
      pluto_constraints_negate_row(ub, 0);
      pluto_constraints_add(ubs, ub);
      pluto_constraints_free(ub);
    }
    if (dup->val[i][0] >= 1) {
      /* Lower bound */
      lb = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(lbs, lb);
      pluto_constraints_free(lb);
    } else if (dup->val[i][0] <= -1) {
      /* Upper bound */
      ub = pluto_constraints_select_row(dup, i);
      pluto_constraints_add(ubs, ub);
      pluto_constraints_free(ub);
    }
  }

  assert(lbs->nrows >= 1);
  assert(ubs->nrows >= 1);
  pluto_constraints_free(dup);

  lbe = get_func_of_expr(lbs, 0, MAXF, params);
  ube = get_func_of_expr(ubs, 0, MINF, params);

  *lbexpr = lbe;
  *ubexpr = ube;

  pluto_constraints_free(lbs);
  pluto_constraints_free(ubs);
}

/*
 * Get expression for difference of upper and lower bound of pos^th variable
 * in cst in terms of parameters;  last 'npar' dimensions of cst are treated
 * as parameters; *extent should be freed by 'free'
 */
void get_parametric_extent(const PlutoConstraints *cst, int pos, int npar,
                           const char **params, char **extent,
                           char **p_lbexpr) {
  char *lbexpr, *ubexpr;

  get_lb_ub_expr(cst, pos, npar, params, &lbexpr, &ubexpr);

  if (!strcmp(lbexpr, ubexpr)) {
    *extent = strdup("1");
  } else {
    *extent =
        (char *)malloc(strlen(lbexpr) + strlen(ubexpr) + strlen(" -  + 1") + 1);
    sprintf(*extent, "%s - %s + 1", ubexpr, lbexpr);
  }
  if (p_lbexpr != NULL) {
    *p_lbexpr = (char *)malloc(strlen(lbexpr) + 1);
    strcpy(*p_lbexpr, lbexpr);
  }
  free(lbexpr);
  free(ubexpr);
}

/* Get Alpha matrix (A matrix - INRIA transformation representation */
PlutoMatrix *get_alpha(const Stmt *stmt, const PlutoProg *prog) {
  PlutoMatrix *a;
  a = pluto_matrix_alloc(stmt->dim, stmt->dim, prog->context);

  unsigned r = 0;
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    if (stmt->hyp_types[i] == H_LOOP ||
        stmt->hyp_types[i] == H_TILE_SPACE_LOOP) {
      for (unsigned c = 0; c < stmt->dim; c++) {
        a->val[r][c] = stmt->trans->val[i][c];
      }
      r++;
      if (r == stmt->dim)
        break;
    }
  }

  assert(r == stmt->dim);

  return a;
}

bool pluto_is_hyperplane_scalar(const Stmt *stmt, int level) {
  assert(level <= (int)stmt->trans->nrows - 1);

  for (unsigned j = 0; j < stmt->dim; j++) {
    if (stmt->trans->val[level][j] != 0)
      return false;
  }

  return true;
}

int pluto_is_hyperplane_loop(const Stmt *stmt, int level) {
  return !pluto_is_hyperplane_scalar(stmt, level);
}

/* Get the remapping matrix: maps time iterators back to the domain
 * iterators; divs: divisors for the rows */
PlutoMatrix *pluto_stmt_get_remapping(const Stmt *stmt, int **divs) {
  PlutoMatrix *trans = stmt->trans;
  PlutoMatrix *remap = pluto_matrix_dup(trans);

  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  *divs = (int *)malloc(sizeof(int) * (stmt->dim + npar + 1));

  for (unsigned i = 0; i < remap->nrows; i++) {
    pluto_matrix_negate_row(remap, remap->nrows - 1 - i);
    pluto_matrix_add_col(remap, 0);
    remap->val[trans->nrows - 1 - i][0] = 1;
  }

  /* Bring the stmt iterators to the left */
  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_move_col(remap, remap->nrows + i, i);
  }

  assert(stmt->dim <= remap->nrows);

  for (unsigned i = 0; i < stmt->dim; i++) {
    if (remap->val[i][i] == 0) {
      unsigned k;
      for (k = i + 1; k < remap->nrows; k++) {
        if (remap->val[k][i] != 0)
          break;
      }
      if (k < remap->nrows) {
        pluto_matrix_interchange_rows(remap, i, k);
      } else {
        /* Can't associate domain iterator with time iterator */
        /* Shouldn't happen with a full-ranked transformation */
        printf("Can't associate domain iterator #%d with time iterators\n",
               i + 1);
        pluto_matrix_print(stdout, remap);
        assert(0);
      }
    }
    assert(remap->val[i][i] != 0);
    for (unsigned k = i + 1; k < remap->nrows; k++) {
      if (remap->val[k][i] == 0)
        continue;
      int64_t _lcm = lcm(remap->val[k][i], remap->val[i][i]);
      int64_t factor1 = _lcm / remap->val[k][i];
      for (unsigned j = i; j < remap->ncols; j++) {
        remap->val[k][j] = remap->val[k][j] * factor1 -
                           remap->val[i][j] * (_lcm / remap->val[i][i]);
      }
    }
  }

  /* Solve upper triangular system now */
  for (int i = stmt->dim - 1; i >= 0; i--) {
    assert(remap->val[i][i] != 0);
    for (int k = i - 1; k >= 0; k--) {
      if (remap->val[k][i] == 0)
        continue;
      int64_t _lcm = lcm(remap->val[k][i], remap->val[i][i]);
      int64_t factor1 = _lcm / remap->val[k][i];
      for (unsigned j = 0; j < remap->ncols; j++) {
        remap->val[k][j] = remap->val[k][j] * (factor1)-remap->val[i][j] *
                           (_lcm / remap->val[i][i]);
      }
    }
  }

  assert(remap->nrows >= stmt->dim);
  for (int i = remap->nrows - 1; i >= (int)stmt->dim; i--) {
    pluto_matrix_remove_row(remap, remap->nrows - 1);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    assert(remap->val[i][i] != 0);
    if (remap->val[i][i] <= -1) {
      pluto_matrix_negate_row(remap, i);
    }
    (*divs)[i] = llabs(remap->val[i][i]);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_remove_col(remap, 0);
  }

  for (unsigned i = 0; i < stmt->dim; i++) {
    pluto_matrix_negate_row(remap, i);
  }

  /* Identity for the parameter and constant part */
  for (int i = 0; i < npar + 1; i++) {
    pluto_matrix_add_row(remap, remap->nrows);
    remap->val[remap->nrows - 1][remap->ncols - npar - 1 + i] = 1;
    (*divs)[stmt->dim + i] = 1;
  }

  assert(remap->nrows == stmt->dim + npar + 1 && "Invalid remapping matrix");

  return remap;
}

void pluto_prog_params_print(const PlutoProg *prog) {
  for (int i = 0; i < prog->npar; i++) {
    printf("%s\n", prog->params[i]);
  }
}

/// Constructs the new access function for `acc' resulting from the
/// transformation attached to 'stmt'.
PlutoMatrix *pluto_get_new_access_func(const PlutoMatrix *acc, const Stmt *stmt,
                                       int **divs) {
  int npar = stmt->domain->ncols - (int)stmt->dim - 1;

  assert(acc->ncols == stmt->dim + npar + 1 &&
         "Invalid access function matrix");

  *divs = (int *)malloc(sizeof(int) * acc->nrows);

  int *remap_divs;
  PlutoMatrix *remap = pluto_stmt_get_remapping(stmt, &remap_divs);
  assert(remap->nrows == stmt->dim + npar + 1);

  int _lcm = 1;
  for (unsigned r = 0; r < remap->nrows; r++) {
    assert(remap_divs[r] != 0);
    _lcm = lcm(_lcm, remap_divs[r]);
  }
  for (unsigned r = 0; r < remap->nrows; r++) {
    for (unsigned c = 0; c < remap->ncols; c++) {
      remap->val[r][c] = (remap->val[r][c] * _lcm) / remap_divs[r];
    }
  }

  PlutoMatrix *newacc = pluto_matrix_product(acc, remap);

  for (unsigned r = 0; r < newacc->nrows; r++) {
    (*divs)[r] = _lcm;
  }

  assert(newacc->ncols == stmt->trans->nrows + npar + 1);

  pluto_matrix_free(remap);
  free(remap_divs);

  return newacc;
}

/* Separates a list of statements at level 'level' */
void pluto_separate_stmts(PlutoProg *prog, Stmt **stmts, int num, int level,
                          int offset) {
  int i, nstmts, k;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
  }
  for (k = 0; k < num; k++) {
    stmts[k]->trans->val[level][stmts[k]->trans->ncols - 1] = offset + 1 + k;
  }

  pluto_prog_add_hyperplane(prog, level, H_SCALAR);
  prog->hProps[level].dep_prop = SEQ;
}

/* Separates a statement from the rest (places it later) at level 'level';
 * this is done by inserting a scalar dimension separating them */
void pluto_separate_stmt(PlutoProg *prog, const Stmt *stmt, int level) {
  int i, nstmts;

  nstmts = prog->nstmts;

  for (i = 0; i < nstmts; i++) {
    pluto_stmt_add_hyperplane(prog->stmts[i], H_SCALAR, level);
  }
  stmt->trans->val[level][stmt->trans->ncols - 1] = 1;

  pluto_prog_add_hyperplane(prog, level, H_SCALAR);
  prog->hProps[level].dep_prop = SEQ;
}

int pluto_stmt_is_member_of(int stmt_id, Stmt **slist, int len) {
  int i;
  for (i = 0; i < len; i++) {
    if (stmt_id == slist[i]->id)
      return 1;
  }
  return 0;
}

int pluto_stmt_is_subset_of(Stmt **s1, int n1, Stmt **s2, int n2) {
  int i;

  for (i = 0; i < n1; i++) {
    if (!pluto_stmt_is_member_of(s1[i]->id, s2, n2))
      return 0;
  }

  return 1;
}

/* Add new to accs if it's an access to a variable not already contained in
 * accs */
void add_if_new_var(PlutoAccess ***accs, int *num, PlutoAccess *new_acc) {
  int i;

  for (i = 0; i < *num; i++) {
    if (!strcmp((*accs)[i]->name, new_acc->name)) {
      break;
    }
  }

  if (i == *num) {
    *accs = (PlutoAccess **)realloc(*accs, (*num + 1) * sizeof(PlutoAccess *));
    (*accs)[*num] = new_acc;
    (*num)++;
  }
}

/* Get all write accesses in the program */
PlutoAccess **pluto_get_all_waccs(const PlutoProg *prog, int *num) {
  PlutoAccess **accs = NULL;
  *num = 0;

  for (unsigned i = 0; i < prog->nstmts; i++) {
    assert(prog->stmts[i]->nwrites == 1);
    add_if_new_var(&accs, num, prog->stmts[i]->writes[0]);
  }
  return accs;
}

int pluto_get_max_ind_hyps_non_scalar(const PlutoProg *prog) {
  int max = 0;

  for (unsigned i = 0; i < prog->nstmts; i++) {
    max = PLMAX(max, pluto_stmt_get_num_ind_hyps_non_scalar(prog->stmts[i]));
  }

  return max;
}

/*
 * The maximum number of linearly independent hyperplanes across all
 * statements
 */
int pluto_get_max_ind_hyps(const PlutoProg *prog) {
  unsigned max = 0;

  for (unsigned i = 0; i < prog->nstmts; i++) {
    max = PLMAX(max, pluto_stmt_get_num_ind_hyps(prog->stmts[i]));
  }

  return max;
}

int pluto_stmt_get_num_ind_hyps_non_scalar(const Stmt *stmt) {
  PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

  /* Ignore padding dimensions, params, and constant part */
  for (unsigned i = stmt->dim_orig; i < stmt->trans->ncols; i++) {
    pluto_matrix_remove_col(tprime, stmt->dim_orig);
  }
  unsigned j = 0;
  for (unsigned i = 0; i < stmt->trans->nrows; i++) {
    if (stmt->hyp_types[i] == H_SCALAR) {
      pluto_matrix_remove_row(tprime, i - j);
      j++;
    }
  }

  unsigned isols = pluto_matrix_get_rank(tprime);
  pluto_matrix_free(tprime);

  return isols;
}

unsigned pluto_stmt_get_num_ind_hyps(const Stmt *stmt) {
  PlutoMatrix *tprime = pluto_matrix_dup(stmt->trans);

  /* Ignore padding dimensions, params, and constant part */
  for (unsigned i = stmt->dim_orig; i < stmt->trans->ncols; i++) {
    pluto_matrix_remove_col(tprime, stmt->dim_orig);
  }

  unsigned isols = pluto_matrix_get_rank(tprime);
  pluto_matrix_free(tprime);

  return isols;
}

/*
 * Are all transformations full column-ranked?
 */
int pluto_transformations_full_ranked(PlutoProg *prog) {
  for (unsigned i = 0; i < prog->nstmts; i++) {
    if (pluto_stmt_get_num_ind_hyps(prog->stmts[i]) <
        prog->stmts[i]->dim_orig) {
      return 0;
    }
  }

  return 1;
}

struct acc_info {
  char *prefix;
  int acc_num;
  isl_union_map **new_maps;
  isl_union_map **schedule;
  isl_map *base_schedule;
};

/*
 * Return clone of a statement
 */
Stmt *pluto_stmt_dup(const Stmt *stmt) {
  Stmt *nstmt = pluto_stmt_alloc(stmt->dim, stmt->domain, stmt->trans);

  nstmt->dim_orig = stmt->dim_orig;
  nstmt->type = stmt->type;

  for (unsigned i = 0; i < stmt->dim; i++) {
    nstmt->iterators[i] = strdup(stmt->iterators[i]);
    nstmt->is_orig_loop[i] = stmt->is_orig_loop[i];
  }
  if (stmt->text)
    nstmt->text = strdup(stmt->text);

  nstmt->nreads = stmt->nreads;
  nstmt->nwrites = stmt->nwrites;

  nstmt->reads = (PlutoAccess **)malloc(nstmt->nreads * sizeof(PlutoAccess *));
  nstmt->writes =
      (PlutoAccess **)malloc(nstmt->nwrites * sizeof(PlutoAccess *));

  for (int i = 0; i < stmt->nreads; i++) {
    nstmt->reads[i] = pluto_access_dup(stmt->reads[i]);
  }

  for (int i = 0; i < stmt->nwrites; i++) {
    nstmt->writes[i] = pluto_access_dup(stmt->writes[i]);
  }

  return nstmt;
}

static void decrement_stmt_id(PlutoProg *prog, int id) {
  int i;

  prog->stmts[id]->id--;

  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    if (dep->src == id) {
      dep->src--;
    }
    if (dep->dest == id) {
      dep->dest--;
    }
  }
}

/* Add statement to program; can't reuse arg stmt pointer any more */
void pluto_remove_stmt(PlutoProg *prog, unsigned stmt_id) {
  assert(prog->nstmts > 0 && "no stmts in program");

  pluto_stmt_free(prog->stmts[stmt_id]);

  for (unsigned i = stmt_id; i < prog->nstmts - 1; i++) {
    prog->stmts[i] = prog->stmts[i + 1];
    decrement_stmt_id(prog, prog->stmts[i]->id);
  }

  prog->nstmts--;

  prog->stmts =
      (Stmt **)realloc(prog->stmts, ((prog->nstmts) * sizeof(Stmt *)));

  for (unsigned i = 0; i < prog->nstmts; i++) {
    prog->nvar = PLMAX(prog->nvar, (int)prog->stmts[i]->dim);
  }
}

void pluto_transformations_pretty_print(const PlutoProg *prog) {
  for (unsigned i = 0; i < prog->nstmts; i++) {
    pluto_stmt_transformation_print(prog->stmts[i]);
  }
}

void pluto_transformation_print_level(const PlutoProg *prog, int level) {
  unsigned nstmts = prog->nstmts;

  for (unsigned i = 0; i < nstmts; i++) {
    fprintf(stdout, "h(S%d) = ", i + 1);
    pluto_stmt_print_hyperplane(stdout, prog->stmts[i], level);
    if (i < nstmts - 1)
      fprintf(stdout, ", ");
  }
  printf("\n");
}

/* List properties of newly found hyperplanes */
void pluto_print_hyperplane_properties(const PlutoProg *prog) {
  int j, numH;
  HyperplaneProperties *hProps;

  hProps = prog->hProps;
  numH = prog->num_hyperplanes;

  if (numH == 0) {
    fprintf(stdout, "No hyperplanes\n");
  }

  /* Note that loop properties are calculated for each dimension in the
   * transformed space (common for all statements) */
  for (j = 0; j < numH; j++) {
    fprintf(stdout, "t%d --> ", j + 1);
    switch (hProps[j].dep_prop) {
    case PARALLEL:
      fprintf(stdout, "parallel ");
      break;
    case SEQ:
      fprintf(stdout, "serial   ");
      break;
    case PIPE_PARALLEL:
      fprintf(stdout, "fwd_dep  ");
      break;
    default:
      fprintf(stdout, "unknown  ");
      break;
    }
    switch (hProps[j].type) {
    case H_LOOP:
      fprintf(stdout, "loop  ");
      break;
    case H_SCALAR:
      fprintf(stdout, "scalar");
      break;
    case H_TILE_SPACE_LOOP:
      fprintf(stdout, "tLoop ");
      break;
    default:
      fprintf(stdout, "unknown  ");
      break;
    }
    fprintf(stdout, " (band %d)", hProps[j].band_num);
    fprintf(stdout, hProps[j].unroll ? "ujam" : "no-ujam");
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
}

void pluto_transformations_print(const PlutoProg *prog) {
  for (unsigned i = 0; i < prog->nstmts; i++) {
    printf("T_(S%d) \n", prog->stmts[i]->id + 1);
    pluto_matrix_print(stdout, prog->stmts[i]->trans);
  }
}

void pluto_stmt_transformation_print(const Stmt *stmt) {
  fprintf(stdout, "T(S%d): ", stmt->id + 1);
  printf("(");
  for (unsigned level = 0; level < stmt->trans->nrows; level++) {
    pluto_stmt_print_hyperplane(stdout, stmt, level);
    if ((int)level <= (int)stmt->trans->nrows - 2)
      printf(", ");
  }
  printf(")\n");

  printf("loop types (");
  for (unsigned level = 0; level < stmt->trans->nrows; level++) {
    if (level > 0)
      printf(", ");
    if (stmt->hyp_types[level] == H_SCALAR)
      printf("scalar");
    else if (stmt->hyp_types[level] == H_LOOP)
      printf("loop");
    else if (stmt->hyp_types[level] == H_TILE_SPACE_LOOP)
      printf("tloop");
    else
      printf("unknown");
  }
  printf(")\n\n");
}

///
PlutoAccess *pluto_create_new_access(int sym_id, const char *name,
                                     PlutoMatrix *acc_mat) {
  PlutoAccess *newacc = (PlutoAccess *)malloc(sizeof(PlutoAccess));
  newacc->sym_id = sym_id;
  newacc->name = strdup(name);
  newacc->mat = acc_mat;
  return newacc;
}

/// Checks if two PlutoAceesess are the same.
static bool are_pluto_accesses_same(PlutoAccess *acc1, PlutoAccess *acc2) {
  // Just compare names in PlutoAccess. sym_id may not be initialized in all
  // cases.
  if (strcmp(acc1->name, acc2->name))
    return false;
  return are_pluto_matrices_equal(acc1->mat, acc2->mat);
}

/// Checks if an access is present in a set of accesses.
static bool is_access_present(std::vector<PlutoAccess *> acc_set,
                              PlutoAccess *acc) {
  for (auto acc_itr = acc_set.begin(); acc_itr != acc_set.end(); acc_itr++) {
    if (are_pluto_accesses_same(*acc_itr, acc))
      return true;
  }
  return false;
}

/// Inserts an access newacc to a vector of accesses if it is not already
/// present.
static void insert_access_if_unique(std::vector<PlutoAccess *> *acc_set,
                                    PlutoAccess *newacc) {
  if (acc_set->empty()) {
    acc_set->push_back(newacc);
  } else if (is_access_present(*acc_set, newacc)) {
    pluto_access_free(newacc);
  } else {
    acc_set->push_back(newacc);
  }
}

/// Returns a vector of unique accesses in the list of input statements. Note
/// that the new access function after the transformation is computed in this
/// function.
static std::vector<PlutoAccess *>
get_unique_accesses_in_stmts(std::vector<Stmt *> stmts, const PlutoProg *prog) {
  std::vector<PlutoAccess *> unique_accesses;
  for (auto stmt_itr = stmts.begin(); stmt_itr != stmts.end(); stmt_itr++) {
    /* All read accesses. */
    Stmt *stmt = *stmt_itr;
    for (int j = 0; j < stmt->nreads; j++) {
      int *divs;
      PlutoAccess *acc = stmt->reads[j];
      PlutoMatrix *new_acc_func =
          pluto_get_new_access_func(acc->mat, stmt, &divs);
      PlutoAccess *newacc =
          pluto_create_new_access(acc->sym_id, acc->name, new_acc_func);
      insert_access_if_unique(&unique_accesses, newacc);
      free(divs);
    }
    /* All write accesses. */
    for (int j = 0; j < stmt->nwrites; j++) {
      int *divs;
      PlutoAccess *acc = stmt->writes[j];
      PlutoMatrix *new_acc_func =
          pluto_get_new_access_func(acc->mat, stmt, &divs);
      PlutoAccess *newacc =
          pluto_create_new_access(acc->sym_id, acc->name, new_acc_func);
      insert_access_if_unique(&unique_accesses, newacc);
      free(divs);
    }
  }
  return unique_accesses;
}

/// Returns the number of unique accesses in a set of statements.
unsigned get_num_unique_accesses_in_stmts(Stmt **stmts, unsigned nstmts,
                                          const PlutoProg *prog) {
  std::vector<Stmt *> stmts_vec;
  for (unsigned i = 0; i < nstmts; i++) {
    stmts_vec.push_back(stmts[i]);
  }
  std::vector<PlutoAccess *> unique_accesses =
      get_unique_accesses_in_stmts(stmts_vec, prog);
  unsigned num_unique_accesses = unique_accesses.size();
  for (auto acc_itr = unique_accesses.begin(); acc_itr != unique_accesses.end();
       acc_itr++) {
    pluto_access_free(*acc_itr);
  }
  return num_unique_accesses;
}

/// Checks if an access is invariant at a given depth. This is very similar to
/// the routine is_invariant in src/polyloop.c but it does not compute the new
/// access function. It assumes that the new access function is already
/// computed.
static bool is_access_invariant(PlutoAccess *acc, unsigned depth) {
  for (unsigned i = 0; i < acc->mat->nrows; i++) {
    if (acc->mat->val[i][depth] != 0)
      return false;
  }
  return true;
}

/// Returns the number of access in the input statements that are invariant
/// with respect to a loop at a perticular depth.
unsigned get_num_invariant_accesses_in_stmts(Stmt **stmts, unsigned nstmts,
                                             unsigned depth,
                                             const PlutoProg *prog) {
  std::vector<Stmt *> stmts_vec;
  for (unsigned i = 0; i < nstmts; i++) {
    stmts_vec.push_back(stmts[i]);
  }
  std::vector<PlutoAccess *> unique_accesses =
      get_unique_accesses_in_stmts(stmts_vec, prog);
  unsigned num_invariant_access = 0;
  for (auto itr = unique_accesses.begin(); itr != unique_accesses.end();
       itr++) {
    if (is_access_invariant(*itr, depth)) {
      num_invariant_access++;
    }
    pluto_access_free(*itr);
  }
  return num_invariant_access;
}
