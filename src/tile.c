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
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "pluto.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"

/* Read tile sizes from file tile.sizes */
static int read_tile_sizes(int *tile_sizes, int *l2_tile_size_ratios,
                           int num_tile_dims, Stmt **stmts, int nstmts,
                           int firstLoop) {
  int i, j;

  FILE *tsfile = fopen("tile.sizes", "r");

  if (!tsfile)
    return 0;

  IF_DEBUG(printf("[pluto] Reading %d tile sizes\n", num_tile_dims););

  if (options->ft >= 0 && options->lt >= 0) {
    num_tile_dims = options->lt - options->ft + 1;
  }

  for (i = 0; i < num_tile_dims && !feof(tsfile); i++) {
    for (j = 0; j < nstmts; j++) {
      if (pluto_is_hyperplane_loop(stmts[j], firstLoop + i))
        break;
    }
    int loop = (j < nstmts);
    if (loop) {
      fscanf(tsfile, "%d", &tile_sizes[i]);
    } else {
      /* Size set for scalar dimension doesn't matter */
      tile_sizes[i] = 42;
    }
  }

  if (i < num_tile_dims) {
    printf("WARNING: not enough tile sizes provided\n");
    fclose(tsfile);
    return 0;
  }

  i = 0;
  while (i < num_tile_dims && !feof(tsfile)) {
    fscanf(tsfile, "%d", &l2_tile_size_ratios[i++]);
  }

  if (i < num_tile_dims) {
    if (options->l2tile)
      printf("WARNING: not enough L2 tile sizes provided; using default\n");
    for (i = 0; i < num_tile_dims; i++) {
      l2_tile_size_ratios[i] = 8;
    }
  }

  fclose(tsfile);
  return 1;
}

/*
 * Reschedule a diamond tile
 */
int pluto_diamond_tile_reschedule(PlutoProg *prog) {
  int i, j, tmp, retval;

  retval = 0;

  for (i = 0; i < prog->nstmts; i++) {
    if (prog->stmts[i]->evicted_hyp) {
      int evicted_hyp_pos = prog->stmts[i]->evicted_hyp_pos;
      int fl = prog->stmts[i]->dim - prog->stmts[i]->dim_orig;
      PlutoMatrix *evicted_hyp = prog->stmts[i]->evicted_hyp;
      assert(fl + evicted_hyp->ncols == prog->stmts[i]->trans->ncols);
      for (j = 0; j < evicted_hyp->ncols; j++) {
        tmp = evicted_hyp->val[0][j];
        evicted_hyp->val[0][j] =
            prog->stmts[i]->trans->val[fl + evicted_hyp_pos][fl + j];
        prog->stmts[i]->trans->val[fl + evicted_hyp_pos][fl + j] = tmp;
      }
      retval = 1;
    }
  }

  return retval;
}

Array *get_corrs_array(PlutoAccess *access, PlutoProg *prog) {

  int k = 0;

  // Search the find the array
  for (k = 0; k < prog->narrays; k++) {
    if (!strcmp(prog->arrays[k]->text, access->name)) {
      return prog->arrays[k];
    }
  }

  return NULL;
}

/* pos: position of domain iterator
 * time_pos: position of time iterator; iter: domain iterator; supply -1
 * if you don't want a scattering function row added for it */
void pluto_array_add_dim(Array *arr, PlutoConstraints *domain, int copy_level,
                         int pos, int time_pos, PlutoProg *prog, char *iter) {
  int i;

  PlutoMatrix *trans = arr->trans;

  int dim = arr->dim;

  assert(pos <= dim);
  assert(time_pos <= trans->nrows);

  pluto_constraints_add_dim(domain, copy_level + pos, NULL);

  if (!arr->dim_updated) {

    arr->dim++;
    arr->num_cur_tiled_loops++;
    arr->iterators =
        (char **)realloc(arr->iterators, (arr->dim) * sizeof(char *));
    for (i = arr->dim - 2; i >= pos; i--) {
      arr->iterators[i + 1] = arr->iterators[i];
    }
    arr->iterators[pos] = strdup(iter);

    /* Array should always have a transformation */
    assert(trans != NULL);
    pluto_matrix_add_col(trans, pos);

    if (time_pos != -1) {
      pluto_matrix_add_row(trans, time_pos);
      trans->val[time_pos][pos] = 1;
    }
  }

  return;
}

void pluto_dist_tile_array_dim(Array *arr, PlutoConstraints *domain,
                               int *dom_tiled_loops, int copy_level,
                               int tile_size, int depth, int dim_orig,
                               PlutoProg *prog) {

  PlutoMatrix *trans = arr->trans;

  int j = 0, k = 0;

  /* 1. Specify tiles in the original domain.
   * NOTE: tile shape info comes in here */

  /* 1.1 Add additional dimensions */
  char iter[6];
  sprintf(iter, "zT%d", domain->ncols);
  int npar = prog->npar;

  print_polylib_visual_sets("domain", domain);

  /* 1.2 Specify tile shapes in the original domain */
  int num_domain_supernodes = 0;

  assert(tile_size >= 1);

  /* Domain supernodes aren't added for scalar dimensions */
  pluto_array_add_dim(arr, domain, copy_level, num_domain_supernodes, 0, prog,
                      iter);

  (*dom_tiled_loops)++;
  /* Add relation b/w tile space variable and intra-tile variables like
   * 32*xt <= 2t+i <= 32xt + 31 */
  /* Lower bound */
  pluto_constraints_add_inequality(domain);

  PlutoConstraints *cst = domain;

  int start = num_domain_supernodes + copy_level + *dom_tiled_loops;

  while (cst != NULL) {
    for (j = start, k = 0; k < arr->dim_orig + npar; j++, k++) {
      cst->val[cst->nrows - 1][j] =
          trans->val[dim_orig][arr->num_cur_tiled_loops + k];
    }

    cst->val[cst->nrows - 1][num_domain_supernodes + copy_level] = -tile_size;

    cst->val[cst->nrows - 1][domain->ncols - 1] =
        trans->val[dim_orig][trans->ncols - 1];

    cst = cst->next;
  }

  print_polylib_visual_sets("domain", domain);

  /* Upper bound */
  pluto_constraints_add_inequality(domain);
  cst = domain;

  while (cst != NULL) {
    for (j = start, k = 0; k < arr->dim_orig + npar; j++, k++) {
      cst->val[cst->nrows - 1][j] =
          -trans->val[dim_orig][arr->num_cur_tiled_loops + k];
    }

    cst->val[cst->nrows - 1][num_domain_supernodes + copy_level] = tile_size;

    cst->val[cst->nrows - 1][domain->ncols - 1] =
        -trans->val[dim_orig][trans->ncols - 1] + tile_size - 1;

    cst = cst->next;
  }

  print_polylib_visual_sets("domain", domain);
  num_domain_supernodes++;
}

/* Returns the data tiles required to allocate the data space parametrized on
 * copy level
 */
void pluto_dist_permute_tile_dims(PlutoConstraints *cst, int start, int last) {
  int j = start, k = last;

  assert(j >= 0 && j < cst->ncols);
  assert(k >= 0 && k < cst->ncols);
  assert(j <= k);

  for (j = start, k = last; j < k; ++j, --k) {
    pluto_constraints_interchange_cols(cst, j, k);
  }

  return;
}

void pluto_dist_permute_iterators(char **iterators, int start, int last) {
  int j = start, k = last;

  char *temp;
  assert(j <= k);

  for (j = start, k = last; j < k; ++j, --k) {
    temp = iterators[last];
    iterators[last] = iterators[start];
    iterators[start] = temp;
  }

  return;
}

PlutoConstraints *pluto_dist_get_required_data_tiles(PlutoConstraints *domain,
                                                     int copylevel,
                                                     char *arr_name,
                                                     PlutoProg *prog) {
  int k, j;
  Array *arr = pluto_get_corrs_array(arr_name, prog);
  PlutoConstraints *tiled_domain = pluto_constraints_dup(domain);

  assert(arr != NULL);

  struct TiledHyperplane *t = arr->tiled_hyperplane;

  int trans_row, dom_tiled_loops = 0;

  for (k = 0; k < arr->num_tiled_loops; ++k) {

    if (!arr->dim_updated) {

      /* after each dim tiling,  tiled_dim needs to be updated */
      for (j = 0; j < arr->num_hyperplanes_found; ++j) {

        if (arr->tiled_hyperplane[j].orig_dim >= t[k].firstD)
          arr->tiled_hyperplane[j].orig_dim++;

        if (arr->tiled_hyperplane[j].tiled_dim >= t[k].firstD)
          arr->tiled_hyperplane[j].tiled_dim++;
      }

      arr->tiled_hyperplane[k].tiled_dim = 0;

      arr->first_tile_dim = t[k].firstD;
      arr->last_tile_dim = t[k].firstD + arr->num_tiled_loops - 1;

      arr->copy_level_used = copylevel;
    }

    trans_row = t[k].orig_dim;

    pluto_dist_tile_array_dim(arr, tiled_domain, &dom_tiled_loops, copylevel,
                              t[k].tile_sizes, t[k].depth, trans_row, prog);
  }

  int i = 0;

  if (!arr->dim_updated) {
    arr->dim_updated = 1;
    arr->npar = prog->npar;

    /* Permute the tiled dims and itertors
     *  jj ii i j -> ii jj i j
     */
    pluto_dist_permute_iterators(arr->iterators, arr->first_tile_dim,
                                 arr->last_tile_dim);

    // Update the tilied_dims
    int diff = arr->num_tiled_loops;
    t = arr->tiled_hyperplane;
    for (i = 0; i < arr->dim_orig; ++i) {

      if (t[i].is_tiled)
        t[i].tiled_dim = t[i].orig_dim - diff;
      else
        diff++;
    }
  }

  pluto_dist_permute_tile_dims(tiled_domain, arr->first_tile_dim + copylevel,
                               arr->last_tile_dim + copylevel);

  print_polylib_visual_sets("tiles", tiled_domain);
  /* Project out the intra data tile iterators */
  pluto_constraints_project_out(tiled_domain, copylevel + arr->num_tiled_loops,
                                arr->dim_orig);

  print_polylib_visual_sets("tiles_perm", tiled_domain);

  return tiled_domain;
}

void pluto_dist_compute_all_data_tiles(PlutoProg *prog, int copylevel) {

  Stmt *stmt = prog->stmts[0];
  int i, j, k;

  copylevel = 2;

  for (j = 0; j < prog->narrays; ++j) {

    PlutoConstraints *data = NULL;
    PlutoConstraints *curr = NULL;
    Array *arr = prog->arrays[j];

    char *arr_name = prog->arrays[j]->text;

    for (i = 0; i < stmt->nreads; ++i) {
      PlutoAccess *acc = stmt->reads[i];

      if (strcmp(arr_name, acc->name) != 0)
        continue;

      curr =
          pluto_compute_region_data(stmt, stmt->domain, acc, copylevel, prog);

      print_polylib_visual_sets("x", curr);

      if (data == NULL)
        data = pluto_constraints_dup(curr);
      else
        pluto_constraints_unionize(data, curr);

      pluto_constraints_free(curr);
    }

    for (i = 0; i < stmt->nwrites; ++i) {
      PlutoAccess *acc = stmt->writes[i];

      if (strcmp(arr_name, acc->name) != 0)
        continue;

      curr =
          pluto_compute_region_data(stmt, stmt->domain, acc, copylevel, prog);

      if (data == NULL)
        data = pluto_constraints_dup(curr);
      else
        pluto_constraints_unionize(data, curr);

      pluto_constraints_free(curr);
    }

    pluto_constraints_free(data);

    struct TiledHyperplane *t = arr->tiled_hyperplane;

    int dom_tiled_loops = 0;
    for (k = 0; k < arr->num_tiled_loops; ++k) {
      pluto_dist_tile_array_dim(arr, arr->array_bounds, &dom_tiled_loops,
                                copylevel, t[k].tile_sizes, t[k].depth,
                                t[k].firstD, prog);

      /* add the tiled dim update code here
       * after each dim tiling,  tiled_dim needs to be updated
       */
      for (i = 0; i < arr->num_hyperplanes_found; ++i) {

        if (arr->tiled_hyperplane[i].orig_dim >= t[i].firstD)
          arr->tiled_hyperplane[i].orig_dim++;

        if (arr->tiled_hyperplane[i].tiled_dim >= t[i].firstD)
          arr->tiled_hyperplane[i].tiled_dim++;
      }

      t[k].tiled_dim = t[k].firstD;
      arr->first_tile_dim = t[k].firstD;
      arr->last_tile_dim = t[k].lastD;
    }
    arr->first_tile_dim += copylevel;
    arr->last_tile_dim += copylevel;

    for (i = 0; i < arr->num_tiled_loops; ++i) {
      arr->tiled_hyperplane[i].orig_dim += copylevel;
      arr->tiled_hyperplane[i].tiled_dim += copylevel;
      arr->tiled_hyperplane[i].firstD += copylevel;
      arr->tiled_hyperplane[i].lastD += copylevel;
    }

    print_polylib_visual_sets("d", arr->array_bounds);
  }
}

/* This function returns the size of a data tile
 * data_tile_size = tile_size1 * tile_size2...
 */
char *pluto_dist_get_tile_size(Array *arr) {

  int i = 0;
  char *tile_size = (char *)malloc(256 * sizeof(char));
  tile_size[0] = 0;

  for (i = 0; i < arr->num_tiled_loops; ++i) {
    if (arr->tiled_hyperplane[i].is_tiled) {
      sprintf(tile_size + strlen(tile_size), "%d",
              arr->tiled_hyperplane[i].tile_sizes);
      if (i != arr->num_tiled_loops - 1)
        sprintf(tile_size + strlen(tile_size), "*");
    }
  }

  return tile_size;
}

/* Generate the code for data tile declarations
 */
void pluto_dist_declarations(PlutoProg *prog) {

  int i = 0, j = 0;

  for (i = 0; i < prog->narrays; ++i) {
    Array *arr = prog->arrays[i];

    PlutoConstraints *domain = pluto_constraints_dup(arr->array_bounds);

    print_polylib_visual_sets("set", domain);

    int start = arr->first_tile_dim + arr->num_hyperplanes_found;
    int num_intra_tile_dim = arr->dim_orig;

    print_polylib_visual_sets("set", domain);
    // project out the intra tile dim
    pluto_constraints_project_out(domain, start, num_intra_tile_dim);

    print_polylib_visual_sets("set", domain);

    /* Using a hack here, domains after project out seems to have
     * extra unbounded constraints, so setting next to null
     * TODO: try fixing this or use the array domain extracted from
     * file.
     */

    domain->next = NULL;
    print_polylib_visual_sets("set", domain);
    /* Number of points in the projected out domain gives the
     * total number of data tiles
     */
    assert(domain->ncols == arr->num_tiled_loops + prog->npar + 1);

    char *decl = (char *)malloc(1024 * sizeof(char));
    decl[0] = 0;

    sprintf(decl + strlen(decl), "%s ", arr->data_type);

    sprintf(decl + strlen(decl), "%s_trans", arr->text);

    char *buf_size = get_parametric_bounding_box(
        domain, arr->first_tile_dim + j, arr->num_tiled_loops, prog->npar,
        (const char **)prog->params);

    sprintf(decl + strlen(decl), "[%s]", buf_size);

    free(buf_size);

    sprintf(decl + strlen(decl), ";");
  }
}

void pluto_dist_print_tiles_index_string(char *new_stmt_text, Array *arr,
                                         PlutoProg *prog) {

  int j = 0, k = 0;

  for (j = 0, k = arr->last_tile_dim; j < arr->dim_orig; ++j, k--) {
    if (arr->tiled_hyperplane[j].is_tiled) {

      sprintf(new_stmt_text + strlen(new_stmt_text), "(%s)", arr->iterators[k]);

      if (j != arr->dim_orig - 1)
        sprintf(new_stmt_text + strlen(new_stmt_text), "* %s + ",
                pluto_dist_get_ptr_dim_stride(arr, j, prog));
    }
  }

  return;
}

void pluto_dist_print_tiles_strided_access(char *new_stmt_text, Array *arr,
                                           PlutoProg *prog) {

  int j = 0, k = 0;

  sprintf(new_stmt_text + strlen(new_stmt_text), "[");
  for (j = 0, k = arr->last_tile_dim; j < arr->dim_orig; ++j, k--) {
    if (arr->tiled_hyperplane[j].is_tiled) {

      sprintf(new_stmt_text + strlen(new_stmt_text), "(%s)", arr->iterators[k]);

      if (j != arr->dim_orig - 1)
        sprintf(new_stmt_text + strlen(new_stmt_text), "* %s + ",
                pluto_dist_get_ptr_dim_stride(arr, j, prog));
    } else
      sprintf(new_stmt_text + strlen(new_stmt_text), "[%s]", arr->iterators[k]);
  }
  sprintf(new_stmt_text + strlen(new_stmt_text), "]");
  return;
}

char *pluto_dist_ptr_init_stmt_text(char *arr_name, PlutoProg *prog) {

  Array *arr = pluto_get_corrs_array(arr_name, prog);

  char *new_stmt_text = (char *)malloc(2048 * sizeof(char));
  new_stmt_text[0] = '\0';

  sprintf(new_stmt_text + strlen(new_stmt_text), "%s_trans", arr->text);

  pluto_dist_print_tiles_strided_access(new_stmt_text, arr, prog);

  sprintf(new_stmt_text + strlen(new_stmt_text), " = NULL; \t");

  return new_stmt_text;
}

char *pluto_dist_malloc_stmt_text(char *arr_name, PlutoProg *prog,
                                  int use_strides) {

  int j = 0, k = 0;

  Array *arr = pluto_get_corrs_array(arr_name, prog);

  char *stmt_text = (char *)malloc(1024 * 5 * sizeof(char));
  stmt_text[0] = '\0';

  sprintf(stmt_text + strlen(stmt_text), "if(%s_trans", arr_name);

  if (use_strides) {
    sprintf(stmt_text + strlen(stmt_text), "[");

    for (j = arr->last_tile_dim, k = arr->first_tile_dim;
         j >= arr->first_tile_dim; --j, k++) {

      sprintf(stmt_text + strlen(stmt_text), "(%s)", arr->iterators[j]);

      if (j != arr->first_tile_dim)
        sprintf(stmt_text + strlen(stmt_text), " * (%s) + ",
                pluto_dist_get_ptr_dim_stride_malloc(arr, k, prog));
    }
    sprintf(stmt_text + strlen(stmt_text), "]");
  } else {
    for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
      sprintf(stmt_text + strlen(stmt_text), "[%s]", arr->iterators[j]);
    }
  }

  sprintf(stmt_text + strlen(stmt_text), " == NULL)\t");

  sprintf(stmt_text + strlen(stmt_text), "%s_trans", arr->text);

  if (use_strides) {
    sprintf(stmt_text + strlen(stmt_text), "[");

    pluto_dist_print_tiles_index_string(stmt_text, arr, prog);

    sprintf(stmt_text + strlen(stmt_text), "]");
  } else {
    for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
      sprintf(stmt_text + strlen(stmt_text), "[%s]", arr->iterators[j]);
    }
  }

  sprintf(stmt_text + strlen(stmt_text), " = (%s *) buffer_alloc_atomic(",
          arr->data_type);

  pluto_dist_print_tiles_index_string(stmt_text, arr, prog);

  sprintf(stmt_text + strlen(stmt_text), ", %s * sizeof(%s), buff_mang_%s); ",
          pluto_dist_get_tile_size(arr), arr->data_type, arr_name);

  return stmt_text;
}

void pluto_dist_intial_copy_stmt(PlutoProg *prog) {

  int i = 0, j = 0;

  for (i = 0; i < prog->narrays; i++) {
    Array *arr = prog->arrays[i];

    char *stmt_text = (char *)malloc(1024 * sizeof(char));

    sprintf(stmt_text, "%s_trans", arr->text);

    for (j = arr->last_tile_dim; j >= arr->first_tile_dim; --j) {
      sprintf(stmt_text + strlen(stmt_text), "[%s]", arr->iterators[j]);
    }

    for (j = 0; j < arr->dim - arr->num_tiled_loops; j++) {
      sprintf(stmt_text + strlen(stmt_text), "[%s]",
              arr->iterators[arr->last_tile_dim + 1 + j]);
    }

    sprintf(stmt_text + strlen(stmt_text), " = ");

    sprintf(stmt_text + strlen(stmt_text), "%s", arr->text);

    struct TiledHyperplane h;

    for (j = 0; j < arr->dim_orig; j++) {
      h = arr->tiled_hyperplane[j];
      if (h.is_tiled) {
        sprintf(stmt_text + strlen(stmt_text), "[%s*%d+%s]",
                arr->iterators[h.tiled_dim], h.tile_sizes,
                arr->iterators[h.orig_dim]);
      } else {
        sprintf(stmt_text + strlen(stmt_text), "[%s]",
                arr->iterators[h.orig_dim]);
      }
    }
  }
}

void pluto_dist_array_intialize(PlutoProg *prog) {

  Array *arr;
  int i = 0;

  for (i = 0; i < prog->narrays; i++) {
    arr = prog->arrays[i];

    assert(arr->trans != NULL);
    arr->trans_orig = pluto_matrix_dup(arr->trans);

    if (arr->tiled_hyperplane != NULL)
      continue;

    assert(arr->num_hyperplanes_found > 0);
  }
}
/* Read tile sizes from file tile.sizes */
int pluto_dist_read_tile_sizes(int *tile_sizes, int num_tile_dims) {
  int i = 0;

  FILE *tsfile = fopen("tile.sizes", "r");

  if (tsfile) {
    for (i = 0; i < num_tile_dims && !feof(tsfile); i++) {
      fscanf(tsfile, "%d", &tile_sizes[i]);
    }
  }

  if (i < num_tile_dims) {
    printf("WARNING: not enough tile sizes provided\n");
    for (; i < num_tile_dims; i++)
      tile_sizes[i] = 32;
  }

  if (tsfile)
    fclose(tsfile);

  return 1;
}

void pluto_dist_array_tile(PlutoProg *prog) {

  Array *arr;
  int i = 0, j = 0;

  int *tile_sizes = (int *)malloc(prog->max_array_dim * sizeof(int));

  pluto_dist_read_tile_sizes(tile_sizes, prog->max_array_dim);

  for (i = 0; i < prog->narrays; i++) {
    arr = prog->arrays[i];

    assert(arr->trans != NULL);

    if (arr->trans_orig != NULL)
      arr->trans_orig = pluto_matrix_dup(arr->trans);

    if (arr->tiled_hyperplane != NULL) {
      arr->tiled_hyperplane = (struct TiledHyperplane *)malloc(
          arr->num_hyperplanes_found * sizeof(struct TiledHyperplane));
    }

    for (j = 0; j < arr->num_hyperplanes_found; j++) {
      arr->tiled_hyperplane[j].orig_dim = j;
      arr->tiled_hyperplane[j].is_tiled = 1;
      arr->tiled_hyperplane[j].tiled_dim = -1;
      arr->tiled_hyperplane[j].loop_tiled = 1;
      arr->tiled_hyperplane[j].depth = j;
      arr->tiled_hyperplane[j].firstD = 0;
      arr->tiled_hyperplane[j].lastD = arr->num_hyperplanes_found;
      arr->tiled_hyperplane[j].tile_sizes = tile_sizes[j];
    }

    arr->num_tiled_loops = arr->num_hyperplanes_found;
  }
}

void pluto_dist_reorder_hyperplanes(Array *arr, int firstD, int lastD,
                                    int cur_dim) {

  int i;

  for (i = 0; i < firstD; i++) {
    if (cur_dim + i < arr->dim)
      arr->hyperplane_mapping[cur_dim + i] = i;
  }
}

/* Copy the loop tiling info like tile size and dims into arrays */
void pluto_dist_update_tiling_info(PlutoProg *prog, Band *band,
                                   int *tile_sizes) {
  Array *arr;
  int i = 0, j = 0, s = 0;

  int depth;

  int firstD = band->loop->depth;
  int lastD = band->loop->depth + band->width - 1;

  Stmt *stmt;
  int row = 0;

  for (depth = firstD; depth <= lastD; depth++) {

    for (s = 0; s < band->loop->nstmts; s++) {

      stmt = band->loop->stmts[s];

      for (i = 0; i < stmt->nreads; i++) {

        if (is_access_scalar(stmt->reads[i]))
          continue;

        arr = get_corrs_array(stmt->reads[i], prog);
        assert(arr != NULL);
        row = arr->hyperplane_mapping[depth - firstD];

        if (row < 0) {
          pluto_dist_reorder_hyperplanes(arr, firstD, lastD, depth);

          row = arr->hyperplane_mapping[depth];
          if (row < 0)
            continue;
        }

        if (row >= arr->dim_orig)
          continue;

        if (arr->tiled_hyperplane[row].is_tiled)
          continue;

        arr->tiled_hyperplane[row].is_tiled = 1;
        arr->tiled_hyperplane[row].loop_tiled = 1;
        arr->tiled_hyperplane[row].depth = depth - firstD;
        arr->tiled_hyperplane[row].firstD = 0;
        arr->tiled_hyperplane[row].lastD = lastD - firstD;
        arr->tiled_hyperplane[row].tile_sizes = tile_sizes[depth - firstD];

        arr->num_tiled_loops++;
      }

      /* TODO: duplicate the same for write access also */
      for (i = 0; i < stmt->nwrites; i++) {
        arr = get_corrs_array(stmt->writes[i], prog);
        assert(arr != NULL);
        // correct mapping from comp hyperplane being tiled to array hyperplane
        row = arr->hyperplane_mapping[depth - firstD];

        if (row < 0) {
          pluto_dist_reorder_hyperplanes(arr, firstD, lastD, depth);

          row = arr->hyperplane_mapping[depth];
          if (row < 0)
            continue;
        }

        assert(row < arr->dim_orig);

        if (arr->tiled_hyperplane[row].is_tiled)
          continue;

        arr->tiled_hyperplane[row].is_tiled = 1;
        arr->tiled_hyperplane[row].loop_tiled = 1;
        arr->tiled_hyperplane[row].depth = depth - firstD;
        arr->tiled_hyperplane[row].firstD = 0;
        arr->tiled_hyperplane[row].lastD = lastD - firstD;
        arr->tiled_hyperplane[row].tile_sizes = tile_sizes[depth - firstD];

        arr->num_tiled_loops++;
      }
    }
    row++;
  }

  /* If any arrays tiling info is not updated then initalize to default
   * tilig info
   */

  for (i = 0; i < prog->narrays; ++i) {
    arr = prog->arrays[i];

    int tiled = 0;
    for (j = 0; j < arr->dim_orig; ++j) {
      tiled += arr->tiled_hyperplane[j].is_tiled;
    }

    if (tiled)
      continue;

    for (j = 0; j < arr->dim_orig; ++j) {
      arr->tiled_hyperplane[j].is_tiled = 1;
      arr->tiled_hyperplane[j].loop_tiled = 1;
      arr->tiled_hyperplane[j].tile_sizes = tile_sizes[j];
      arr->tiled_hyperplane[j].depth = 0;
      arr->tiled_hyperplane[j].firstD = 0;
      arr->tiled_hyperplane[j].lastD = arr->dim_orig;

      arr->num_tiled_loops++;
    }
  }
}

/* Manipulates statement domain and transformation to tile scattering
 * dimensions from firstD to lastD */
void pluto_tile_band(PlutoProg *prog, Band *band, int *tile_sizes) {
  int j, s;
  int depth, npar;

  npar = prog->npar;

  int firstD = band->loop->depth;
  int lastD = band->loop->depth + band->width - 1;

  int num_domain_supernodes[band->loop->nstmts];

  for (s = 0; s < band->loop->nstmts; s++) {
    num_domain_supernodes[s] = 0;
  }

  for (depth = firstD; depth <= lastD; depth++) {
    for (s = 0; s < band->loop->nstmts; s++) {
      Stmt *stmt = band->loop->stmts[s];
      /* 1. Specify tiles in the original domain.
       * NOTE: tile shape info comes in here */

      /* 1.1 Add additional dimensions */
      char iter[6];
      sprintf(iter, "zT%d", stmt->dim);

      PlutoHypType hyp_type =
          (stmt->hyp_types[depth + depth - firstD] == H_SCALAR)
              ? H_SCALAR
              : H_TILE_SPACE_LOOP;

      /* 1.2 Specify tile shapes in the original domain */
      if (hyp_type != H_SCALAR) {
        assert(tile_sizes[depth - firstD] >= 1);
        /* Domain supernodes aren't added for scalar dimensions */
        pluto_stmt_add_dim(stmt, num_domain_supernodes[s], depth, iter,
                           hyp_type, prog);
        /* Add relation b/w tile space variable and intra-tile variables like
         * 32*xt <= 2t+i <= 32xt + 31 */
        /* Lower bound */
        pluto_constraints_add_inequality(stmt->domain);

        for (j = num_domain_supernodes[s] + 1; j < stmt->dim + npar; j++) {
          stmt->domain->val[stmt->domain->nrows - 1][j] =
              stmt->trans
                  ->val[firstD + (depth - firstD) + 1 + (depth - firstD)][j];
        }

        stmt->domain->val[stmt->domain->nrows - 1][num_domain_supernodes[s]] =
            -tile_sizes[depth - firstD];

        stmt->domain->val[stmt->domain->nrows - 1][stmt->domain->ncols - 1] =
            stmt->trans
                ->val[(depth - firstD) + 1 + depth][stmt->dim + prog->npar];

        PlutoConstraints *lb =
            pluto_constraints_select_row(stmt->domain, stmt->domain->nrows - 1);
        pluto_update_deps(stmt, lb, prog);
        pluto_constraints_free(lb);

        /* Upper bound */
        pluto_constraints_add_inequality(stmt->domain);
        for (j = num_domain_supernodes[s] + 1; j < stmt->dim + npar; j++) {
          stmt->domain->val[stmt->domain->nrows - 1][j] =
              -stmt->trans
                   ->val[firstD + (depth - firstD) + 1 + (depth - firstD)][j];
        }

        stmt->domain->val[stmt->domain->nrows - 1][num_domain_supernodes[s]] =
            tile_sizes[depth - firstD];

        stmt->domain->val[stmt->domain->nrows - 1][stmt->domain->ncols - 1] =
            -stmt->trans
                 ->val[(depth - firstD) + 1 + depth][stmt->dim + prog->npar] +
            tile_sizes[depth - firstD] - 1;

        PlutoConstraints *ub =
            pluto_constraints_select_row(stmt->domain, stmt->domain->nrows - 1);
        pluto_update_deps(stmt, ub, prog);
        pluto_constraints_free(ub);

        num_domain_supernodes[s]++;

      } else {
        /* Scattering function for tile space iterator is set the
         * same as its associated domain iterator
         * Dimension is not a loop; tile it trivially
         */
        pluto_stmt_add_hyperplane(stmt, H_SCALAR, depth);
        for (j = 0; j < stmt->dim + npar + 1; j++) {
          stmt->trans->val[depth][j] =
              stmt->trans
                  ->val[firstD + (depth - firstD) + 1 + (depth - firstD)][j];
        }
      }
      stmt->num_tiled_loops++;
      stmt->first_tile_dim = firstD;
      stmt->last_tile_dim = lastD;
    } /* all statements */
  }   /* all scats to be tiled */

  if (options->data_dist)
    pluto_dist_update_tiling_info(prog, band, tile_sizes);
}

/*
 * Updates statement domains and transformations to represent the new
 * tiled code. A schedule of tiles is created for parallel execution if
 * --parallel is on. Intra-tile optimization is done as part of this as well.
 */
void pluto_tile(PlutoProg *prog) {
  unsigned nbands, i, j, n_ibands, num_tiled_levels, nloops;
  Band **bands, **ibands;
  bands = pluto_get_outermost_permutable_bands(prog, &nbands);
  ibands = pluto_get_innermost_permutable_bands(prog, &n_ibands);
  IF_DEBUG(printf("[pluto_tile] Outermost tilable bands\n"););
  IF_DEBUG(pluto_bands_print(bands, nbands););
  IF_DEBUG(printf("[pluto_tile] Innermost tilable bands\n"););
  IF_DEBUG(pluto_bands_print(ibands, n_ibands););

  num_tiled_levels = 0;

  /*
   * Create bands for innermost parallel loops to be 1-d tiled
   * if they are not dominated by any band
   */
  Ploop **loops = pluto_get_parallel_loops(prog, &nloops);
  for (i = 0; i < nloops; i++) {
    if (pluto_loop_is_innermost(loops[i], prog)) {
      for (j = 0; j < nbands; j++) {
        if (is_loop_dominated(loops[i], bands[j]->loop, prog))
          break;
      }
      if (j == nbands) {
        bands = (Band **)realloc(bands, (nbands + 1) * sizeof(Band *));
        bands[nbands++] = pluto_band_alloc(loops[i], 1);
      }
    }
  }
  pluto_loops_free(loops, nloops);

  /* Now, we are ready to tile */
  if (options->lt >= 0 && options->ft >= 0) {
    /* User option specified tiling */

    assert(options->ft <= prog->num_hyperplanes - 1);
    assert(options->lt <= prog->num_hyperplanes - 1);
    assert(options->ft <= options->lt);

    /* L1 tiling */
    pluto_tile_scattering_dims(prog, bands, nbands, 0);
    num_tiled_levels++;

    if (options->l2tile) {
      pluto_tile_scattering_dims(prog, bands, nbands, 1);
      num_tiled_levels++;
    }
  } else {
    /* L1 tiling */
    pluto_tile_scattering_dims(prog, bands, nbands, 0);
    num_tiled_levels++;
    if (options->l2tile) {
      /* L2 tiling */
      pluto_tile_scattering_dims(prog, bands, nbands, 1);
      num_tiled_levels++;
    }
  }

  /* Detect properties after tiling */
  pluto_compute_dep_directions(prog);
  pluto_compute_dep_satisfaction(prog);

  if (!options->silent) {
    fprintf(stdout, "[Pluto] After tiling:\n");
    pluto_transformations_pretty_print(prog);
    /* pluto_print_hyperplane_properties(prog); */
  }

  if (options->diamondtile) {
    int retval;
    retval = pluto_diamond_tile_reschedule(prog);

    if (retval) {
      pluto_compute_dep_directions(prog);
      pluto_compute_dep_satisfaction(prog);
      if (!options->silent) {
        printf("[Pluto] After intra_tile reschedule\n");
        pluto_transformations_pretty_print(prog);
      }
    }
  }

  if (options->intratileopt) {
    int retval = 0;
    for (i = 0; i < nbands; i++) {
      retval |=
          pluto_intra_tile_optimize_band(bands[i], num_tiled_levels, prog);
    }
    if (retval) {
      pluto_compute_dep_directions(prog);
      pluto_compute_dep_satisfaction(prog);
      if (!options->silent) {
        printf("[Pluto] After intra-tile optimize\n");
        pluto_transformations_pretty_print(prog);
      }
    }
  }

  if ((options->parallel || options->dynschedule_graph ||
       options->dynschedule) &&
      !options->identity) {

    int retval = pluto_create_tile_schedule(prog, bands, nbands);
    if (retval && !options->silent) {
      printf("[Pluto] After tile scheduling:\n");
      pluto_transformations_pretty_print(prog);
    }
  }
  pluto_bands_free(bands, nbands);
  pluto_bands_free(ibands, n_ibands);
}

/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void pluto_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands,
                                int l2) {
  int i, j, b;
  int depth;
  int tile_sizes[prog->num_hyperplanes];
  int l2_tile_size_ratios[prog->num_hyperplanes];

  Stmt **stmts = prog->stmts;

  for (j = 0; j < prog->num_hyperplanes; j++) {
    tile_sizes[j] = DEFAULT_L1_TILE_SIZE;
    /* L2 cache is around 64 times L1 cache */
    /* assuming 2-d - this tile size has to be eight
     * times the L1 tile size; NOTE: 8 and NOT
     * 8*default_tile_size -- there is a cumulative multiply
     * involved */
    l2_tile_size_ratios[j] = 8;
  }

  for (b = 0; b < nbands; b++) {
    read_tile_sizes(tile_sizes, l2_tile_size_ratios, bands[b]->width,
                    bands[b]->loop->stmts, bands[b]->loop->nstmts,
                    bands[b]->loop->depth);

    if (l2) {
      pluto_tile_band(prog, bands[b], l2_tile_size_ratios);
    } else {
      pluto_tile_band(prog, bands[b], tile_sizes);
    }
  } /* all bands */

  /* Sink everything to the same depth */
  int max = 0, curr;
  for (i = 0; i < prog->nstmts; i++) {
    max = PLMAX(stmts[i]->trans->nrows, max);
  }
  for (i = 0; i < prog->nstmts; i++) {
    curr = stmts[i]->trans->nrows;
    for (j = curr; j < max; j++) {
      pluto_sink_transformation(stmts[i], stmts[i]->trans->nrows, prog);
    }
  }

  curr = prog->num_hyperplanes;
  for (depth = curr; depth < max; depth++) {
    pluto_prog_add_hyperplane(prog, depth, H_UNKNOWN);
  }
  /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
  pluto_detect_hyperplane_types(prog);
  pluto_detect_hyperplane_types_stmtwise(prog);
}

/* Transform a band of dimensions to get a wavefront
 * (a wavefront of tiles typically)
 *
 * Return: true if something was done, false otherwise
 */
bool pluto_create_tile_schedule_band(PlutoProg *prog, Band *band) {
  int i, j, k, depth;

  /* No need to create tile schedule */
  if (pluto_loop_is_parallel(prog, band->loop))
    return false;

  /* A band can have scalar dimensions; it starts from a loop */
  int loop_depths[band->width];
  /* Number of dimensions which are loops for all statements in this
   * band */
  int nloops;

  loop_depths[0] = band->loop->depth;
  nloops = 1;
  for (depth = band->loop->depth + 1; depth < band->loop->depth + band->width;
       depth++) {
    for (j = 0; j < band->loop->nstmts; j++) {
      if (pluto_is_hyperplane_scalar(band->loop->stmts[j], depth))
        break;
    }
    if (j == band->loop->nstmts) {
      /* All of them are loops */
      loop_depths[nloops++] = depth;
    }
  }

  if (nloops <= 1) {
    /* Band doesn't have at least two dimensions for which all
     * statements have loops at those dimensions */
    return false;
  }

  /* Number of inner parallel dims the wavefront will yield */
  int nip_dims;

  /* loop_depths[0...nloops-1] are the depths for which a tile schedule
   * can be created */
  if (prog->options->multipar) {
    /* Full multi-dimensional wavefront */
    nip_dims = nloops - 1;
  } else {
    /* just use the first two to create a tile schedule */
    nip_dims = 1;
  }

  /* Wavefront satisfies all deps, all inner loops in the band will
   * become parallel */
  int first = band->loop->depth;

  if (!options->innerpar) {
    /* Create the wavefront */
    for (i = 0; i < band->loop->nstmts; i++) {
      Stmt *stmt = band->loop->stmts[i];
      for (k = 1; k <= nip_dims; k++) {
        for (j = 0; j < stmt->trans->ncols; j++) {
          stmt->trans->val[first][j] += stmt->trans->val[loop_depths[k]][j];
        }
      }
    }
  }

  IF_DEBUG(printf("[pluto_create_tile_schedule] Created tile schedule for "););
  IF_DEBUG(printf("t%d to t%d\n", first + 1, loop_depths[nip_dims] + 1));

  /* Update dependence satisfaction levels (better to do this instead of
   * a complete complex dep satisfaction check since we know that the tile
   * schedule will satisfy the dependence satisfied by all the dimensions
   * that is a sum of) */
  for (i = 0; i < prog->ndeps; i++) {
    Dep *dep = prog->deps[i];
    /* satvec s should have been computed */
    if (IS_RAR(dep->type))
      continue;
    if (pluto_stmt_is_member_of(prog->stmts[dep->src]->id, band->loop->stmts,
                                band->loop->nstmts) &&
        pluto_stmt_is_member_of(prog->stmts[dep->dest]->id, band->loop->stmts,
                                band->loop->nstmts)) {
      for (k = 1; k <= nip_dims; k++) {
        dep->satvec[first] |= dep->satvec[loop_depths[k]];
        dep->satvec[loop_depths[k]] = 0;
      }
    }
  }
  /* Recompute dep directions ? not needed */

  return true;
}

bool pluto_create_tile_schedule(PlutoProg *prog, Band **bands, int nbands) {
  int i;
  bool retval = 0;

  IF_DEBUG(printf("creating tile schedule for bands: \n"););
  IF_DEBUG(pluto_bands_print(bands, nbands););

  for (i = 0; i < nbands; i++) {
    retval |= pluto_create_tile_schedule_band(prog, bands[i]);
  }

  return retval;
}

/* Find the innermost permutable nest (at least two tilable hyperplanes) */
void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd) {
  int loop, j, lastloop;

  HyperplaneProperties *hProps = prog->hProps;

  lastloop = getDeepestNonScalarLoop(prog);

  if (hProps[lastloop].dep_prop == SEQ) {
    *bandStart = *bandEnd = lastloop;
    return;
  }

  for (loop = prog->num_hyperplanes - 1; loop >= 0; loop--) {
    if (hProps[loop].type == H_SCALAR)
      continue;
    if (hProps[loop].dep_prop == PIPE_PARALLEL ||
        hProps[loop].dep_prop == PARALLEL) {
      j = loop - 1;
      while (j >= 0 &&
             (hProps[j].dep_prop == PIPE_PARALLEL ||
              hProps[j].dep_prop == PARALLEL) &&
             hProps[j].band_num == hProps[loop].band_num &&
             hProps[j].type == H_LOOP) {
        j--;
      }

      if (j <= loop - 2) {
        *bandEnd = loop;
        *bandStart = j + 1;
        return;
      }
    }
  }
  *bandStart = *bandEnd = lastloop;
}
