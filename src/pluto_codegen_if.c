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
 *
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

#include <assert.h>
#include <math.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cloog/cloog.h>

#include "version.h"

#include "ast_transform.h"
#include "constraints.h"
#include "math_support.h"
#include "pluto.h"
#include "program.h"

static int splitLoops(struct clast_stmt *s, int loop_level, int n,
                      struct clast_stmt **split, int stmt_num,
                      CloogState *state, CloogDomain *domain);

static int get_first_point_loop(Stmt *stmt, const PlutoProg *prog) {
  int i, first_point_loop;

  if ((stmt->type != ORIG) && (stmt->type != ORIG_IN_FUNCTION)) {
    for (i = 0; i < prog->num_hyperplanes; i++) {
      if (!pluto_is_hyperplane_scalar(stmt, i)) {
        return i;
      }
    }
    /* No non-scalar hyperplanes */
    return 0;
  }

  for (i = stmt->last_tile_dim + 1; i < stmt->trans->nrows; i++) {
    if (stmt->hyp_types[i] == H_LOOP)
      break;
  }

  if (i < prog->num_hyperplanes) {
    first_point_loop = i;
  } else {
    /* Should come here only if
     * it's a 0-d statement */
    first_point_loop = 0;
  }

  return first_point_loop;
}

/* Generate and print .cloog file from the transformations computed */
void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog) {
  int i, j;

  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;
  int npar = prog->npar;

  IF_DEBUG(printf("[pluto] generating Cloog file...\n"));
  fprintf(fp, "# CLooG script generated automatically by PLUTO %s\n",
          PLUTO_VERSION);
  fprintf(fp, "# language: C\n");
  fprintf(fp, "c\n\n");

  /* Context: setting conditions on parameters */
  PlutoConstraints *ctx = pluto_constraints_dup(prog->context);
  pluto_constraints_intersect_isl(ctx, prog->codegen_context);
  pluto_constraints_print_polylib(fp, ctx);
  pluto_constraints_free(ctx);

  /* Setting parameter names */
  fprintf(fp, "\n1\n");
  for (i = 0; i < npar; i++) {
    fprintf(fp, "%s ", prog->params[i]);
  }
  fprintf(fp, "\n\n");

  fprintf(fp, "# Number of statements\n");
  fprintf(fp, "%d\n\n", nstmts);

  /* Print statement domains */
  for (i = 0; i < nstmts; i++) {
    fprintf(fp, "# S%d (%s)\n", stmts[i]->id + 1, stmts[i]->text);
    pluto_constraints_print_polylib(fp, stmts[i]->domain);
    fprintf(fp, "0 0 0\n\n");
  }

  fprintf(fp, "# we want cloog to set the iterator names\n");
  fprintf(fp, "0\n\n");

  fprintf(fp, "# Number of scattering functions\n");
  if (nstmts >= 1 && stmts[0]->trans != NULL) {
    fprintf(fp, "%d\n\n", nstmts);

    /* Print scattering functions */
    for (i = 0; i < nstmts; i++) {
      fprintf(fp, "# T(S%d)\n", i + 1);
      PlutoConstraints *sched = pluto_stmt_get_schedule(stmts[i]);
      pluto_constraints_print_polylib(fp, sched);
      fprintf(fp, "\n");
      pluto_constraints_free(sched);
    }

    /* Setting target loop names (all stmts have same number of hyperplanes */
    fprintf(fp, "# we will set the scattering dimension names\n");
    fprintf(fp, "%d\n", stmts[0]->trans->nrows);
    for (i = 0; i < stmts[0]->trans->nrows; i++) {
      fprintf(fp, "t%d", i + 1);
      for (j = 0; j < options->scopnum; j++) {
        fprintf(fp, "_");
      }
      fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  } else {
    fprintf(fp, "0\n\n");
  }
}

static void gen_stmt_macro(const Stmt *stmt, FILE *outfp, PlutoProg *prog) {
  int j;

  for (j = 0; j < stmt->dim; j++) {
    if (stmt->iterators[j] == NULL) {
      printf("Iterator name not set for S%d; required \
                    for generating declarations\n",
             stmt->id + 1);
      assert(0);
    }
  }
  fprintf(outfp, "#define S%d", stmt->id + 1);
  fprintf(outfp, "(");
  for (j = 0; j < stmt->dim; j++) {
    if (j != 0)
      fprintf(outfp, ",");
    fprintf(outfp, "%s", stmt->iterators[j]);
  }
  fprintf(outfp, ")\t");

  /* Generate pragmas for Bee/Cl@k */
  if (options->bee) {
    fprintf(outfp, " __bee_schedule");
    for (j = 0; j < stmt->trans->nrows; j++) {
      fprintf(outfp, "[");
      pluto_affine_function_print(outfp, stmt->trans->val[j], stmt->dim,
                                  (const char **)stmt->iterators);
      fprintf(outfp, "]");
    }
    fprintf(outfp, " _NL_DELIMIT_ ");
  }

  if ((stmt->type == ORIG_IN_FUNCTION || stmt->type == ORIG) &&
      stmt->inner_loop_vec) {

    char *decls = (char *)malloc(1024 * sizeof(char));
    char *mod_stmt_text;
    decls[0] = '\0';
    if (options->data_dist) {
      mod_stmt_text = pluto_dist_modify_stmt_text(stmt->text, 1, prog);
      char *new_stmt_text =
          (char *)malloc(strlen(mod_stmt_text) * 500 * sizeof(char));
      new_stmt_text[0] = '\0';

      int ret = pluto_dist_apply_mod_optimization(mod_stmt_text, new_stmt_text,
                                                  decls, stmt, prog);
      if (ret) {
        fprintf(outfp, "%s\n", new_stmt_text);
      } else
        fprintf(outfp, "%s\n", mod_stmt_text);

      fprintf(outfp, "#define decl_S%d", stmt->id + 1);
      fprintf(outfp, "(");
      for (j = 0; j < stmt->dim; j++) {
        if (j != 0)
          fprintf(outfp, ",");
        fprintf(outfp, "%s", stmt->iterators[j]);
      }
      fprintf(outfp, ")\t");
      fprintf(outfp, "%s\n", decls);

      fprintf(outfp, "#define orig_S%d", stmt->id + 1);
      fprintf(outfp, "(");
      for (j = 0; j < stmt->dim; j++) {
        if (j != 0)
          fprintf(outfp, ",");
        fprintf(outfp, "%s", stmt->iterators[j]);
      }
      fprintf(outfp, ")\t");
      fprintf(outfp, "%s\n", mod_stmt_text);

    } else {
      fprintf(outfp, "%s\n", stmt->text);
    }
  } else if (options->data_dist && stmt->type != LW_COPY_IN &&
             stmt->type != DATA_DIST_COPY) {
    int use_stride = 1;
    fprintf(outfp, "%s\n",
            pluto_dist_modify_stmt_text(stmt->text, use_stride, prog));
  } else
    fprintf(outfp, "%s\n", stmt->text);
}

/* scan and copy till "[" is encountered
*/

void pluto_dist_copy_till_array(char *src, char *dest, int *src_index,
                                int *dest_index) {

  assert(src != NULL);
  assert(dest != NULL);

  while (1) {
    if (src[*src_index] == '\0') {
      dest[*dest_index] = '\0';
      break;
    }

    if (src[*src_index] == '[') {
      dest[*dest_index] = '\0';
      break;
    }

    dest[*dest_index] = src[*src_index];
    (*src_index)++;
    (*dest_index)++;
  }

  return;
}

/* get the array name from the stmt text
*/

char *pluto_dist_get_array_string(char *src, int src_index) {

  assert(src[src_index] == '[');
  // consume '['
  src_index--;

  char *arr_name = (char *)malloc(1024 * sizeof(char));
  char *oper = (char *)malloc(1024 * sizeof(char));
  char c;

  while (src_index > 0) {
    c = src[src_index];

    if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
          (c >= '0' && c <= '9') || c == '_'))
      break;

    src_index--;
  }

  if (src_index > 0)
    sscanf(src + src_index, "%[^a-zA-Z0-9]%[a-zA-Z0-9_][", oper, arr_name);
  //		sscanf(src + src_index, "%s%[a-zA-Z0-9][", oper,  arr_name);
  else
    sscanf(src + src_index, "%[a-zA-Z0-9_][", arr_name);

  return arr_name;
}

void test() {
  int i, j;

  char test[128] = "C[i][j]= alpha*A[i][j] + B[i][j];";
  pluto_dist_get_array_string(test, 1);
  pluto_dist_get_array_string(test, 16);
  //	pluto_dist_get_array_string(test,19);
  for (i = 0; i < 20; ++i) {
    for (j = 0; j < 10; ++j) {
      printf("A[%d][%d] = A_trans[%d][%d][%d][%d]\t\n", i, j, (i + j) / 5,
             i / 5, (i + j) % 5, i % 5);
    }
  }
}

/* return the stride length of tile from curr_dim + 1
 * stride0 = tile_size1 * tile_size2*....
 */
// TODO: handle the case when the data dim is not tiled

char *pluto_dist_get_tile_dim_stride(Array *arr, int curr_dim,
                                     PlutoProg *prog) {
  char *stride = (char *)malloc(256 * sizeof(char));
  stride[0] = '\0';

  sprintf(stride + strlen(stride), "(");
  int i, first = 1;
  for (i = curr_dim + 1; i < arr->dim_orig; ++i) {
    if (arr->tiled_hyperplane[i].is_tiled) {
      if (first)
        sprintf(stride + strlen(stride), "%d",
                arr->tiled_hyperplane[i].tile_sizes);
      else
        sprintf(stride + strlen(stride), " * %d",
                arr->tiled_hyperplane[i].tile_sizes);

      first = 0;
    }
  }
  sprintf(stride + strlen(stride), ")");

  return stride;
}

char *pluto_dist_get_ptr_dim_stride_malloc(Array *arr, int curr_dim,
                                           PlutoProg *prog) {
  char *stride = (char *)malloc(256 * sizeof(char));
  stride[0] = '\0';

  sprintf(stride + strlen(stride), "(");
  int i, first = 1;

  for (i = curr_dim + 1; i <= arr->last_tile_dim; ++i) {
    if (first)
      sprintf(stride + strlen(stride), "%s_size_%d", arr->text, i);
    else
      sprintf(stride + strlen(stride), " * %s_size_%d", arr->text, i);

    first = 0;
  }

  sprintf(stride + strlen(stride), ")");

  return stride;
}

char *pluto_dist_get_ptr_dim_stride(Array *arr, int curr_dim, PlutoProg *prog) {
  char *stride = (char *)malloc(256 * sizeof(char));
  stride[0] = '\0';

  sprintf(stride + strlen(stride), "(");
  int i, first = 1;

  for (i = curr_dim + 1; i < arr->dim_orig; ++i) {
    if (arr->tiled_hyperplane[i].is_tiled) {
      if (first)
        sprintf(stride + strlen(stride), " %s_size_%d ", arr->text,
                arr->tiled_hyperplane[i].tiled_dim);
      else
        sprintf(stride + strlen(stride), " * %s_size_%d ", arr->text,
                arr->tiled_hyperplane[i].tiled_dim);

      first = 0;
    }
  }
  sprintf(stride + strlen(stride), ")");

  return stride;
}

int pluto_dist_get_innermost_loop(const Stmt *stmt, const PlutoProg *prog) {

  PlutoMatrix *trans = stmt->trans;
  int i, j;

  //	pluto_matrix_print(stdout, trans);

  int innermost_loop = -1;
  for (i = trans->nrows - 1; i >= 0; i--) {
    for (j = 0; j < trans->ncols - prog->npar - 1; j++) {

      if (trans->val[i][j] != 0) {
        return j;
      }
    }
  }

  return innermost_loop;
}

char *replacestr(char *string, char *sub, char *replace) {
  if (!string || !sub || !replace)
    return NULL;
  char *pos = string;
  int found = 0;
  while ((pos = strstr(pos, sub))) {
    pos += strlen(sub);
    found++;
  }
  if (found == 0)
    return NULL;
  int size =
      ((strlen(string) - (strlen(sub) * found)) + (strlen(replace) * found)) +
      1;
  char *result = (char *)malloc(size);
  result[0] = 0;
  pos = string;
  char *pos1;
  while ((pos1 = strstr(pos, sub))) {
    int len = (pos1 - pos);
    strncat(result, pos, len);
    strcat(result, replace);
    pos = (pos1 + strlen(sub));
  }
  if (pos != (string + strlen(string)))
    strncat(result, pos, (string - pos));
  return result;
}

void pluto_dist_get_array_indcies(char *src, int *src_index, char *div_ind,
                                  char *mod_ind) {

  assert(src != NULL);

  int mod_ind_len = 0;
  int div_ind_len = 0;

  if (src[*src_index] != '[')
    return;
  (*src_index)++;

  while (1) {
    if (src[*src_index] == '\0') {
      div_ind[div_ind_len] = '\0';
      break;
    }

    if (src[*src_index] == ']') {

      (*src_index)++;
      div_ind[div_ind_len] = '\0';
      break;
    }

    div_ind[div_ind_len] = src[*src_index];
    (*src_index)++;
    div_ind_len++;
  }

  if (src[*src_index] != '[')
    return;
  (*src_index)++;

  while (1) {
    if (src[*src_index] == '\0') {
      mod_ind[mod_ind_len] = '\0';
      break;
    }

    if (src[*src_index] == ']') {
      (*src_index)++;
      mod_ind[mod_ind_len] = '\0';
      break;
    }

    mod_ind[mod_ind_len] = src[*src_index];
    (*src_index)++;
    mod_ind_len++;
  }

  return;
}

int pluto_dist_apply_mod_optimization(char *stmt_text, char *new_stmt_text,
                                      char *decls, const Stmt *stmt,
                                      PlutoProg *prog) {

  if (!options->data_tile_opt)
    return 0;

  int inner_loop = pluto_dist_get_innermost_loop(stmt, prog);
  if (inner_loop < 0)
    return 0;

  char *inner_iter = strdup(stmt->iterators[inner_loop]);

  char *incr = (char *)malloc(1024 * sizeof(char));
  incr[0] = 0;

  char *div_ind = (char *)malloc(1024 * sizeof(char));
  div_ind[0] = 0;
  char *mod_ind = (char *)malloc(1024 * sizeof(char));
  mod_ind[0] = 0;

  int src_index = 0, dest_index = 0;

  char i1[50] = "", i2[50] = "", size_1[50] = "", ci1[50], ci2[50];
  int t1, t2, t3, t4, ret = 0;

  char *cur_decl = (char *)malloc(1024 * sizeof(char));
  cur_decl[0] = 0;

  char **list_indices = prog->decl_indices;
  int ind_size = prog->num_decl_indices;
  int max_pos_offset = 0;
  int max_neg_offset = 0;

  while (1) {

    pluto_dist_copy_till_array(stmt_text, new_stmt_text, &src_index,
                               &dest_index);
    if (src_index >= strlen(stmt_text))
      break;

    char *array_name = pluto_dist_get_array_string(stmt_text, src_index);

    array_name[strlen(array_name) - strlen("_trans")] = 0;

    Array *arr = pluto_get_corrs_array(array_name, prog);

    assert(arr != NULL);

    dest_index -= strlen("_trans");
    new_stmt_text[dest_index] = 0;

    pluto_dist_get_array_indcies(stmt_text, &src_index, div_ind, mod_ind);

    i1[0] = 0;
    i2[0] = 0;
    ci1[0] = 0;
    ci2[0] = 0;
    int i1_pos_offset = 0, i2_pos_offset = 0;
    int i1_neg_offset = 0, i2_neg_offset = 0;

    if (arr->dim_orig > 2)
      return 0;

    if (arr->dim_orig == 1) {
      sscanf(div_ind, "((%d)*( %s ))/%d", &t3, i2, &t4);

      cur_decl[0] = 0;

      sprintf(cur_decl + strlen(cur_decl), "%s_%d", i2, stmt->id);
    } else if (arr->dim_orig == 2) {

      sscanf(div_ind, "((%d)*( %s ))/%d* ( %s ) + ((%d)*( %s ))/%d", &t1, ci1,
             &t2, size_1, &t3, ci2, &t4);

      cur_decl[0] = 0;

      // check for offsets
      ret = sscanf(ci1, "%[a-zA-z0-9]-%d", i1, &i1_neg_offset);
      if (ret == 2 && !strcmp(i1, inner_iter)) {
        max_neg_offset =
            i1_neg_offset > max_neg_offset ? i1_neg_offset : max_neg_offset;
      }
      ret = sscanf(ci1, "%[a-zA-z0-9]+%d", i1, &i1_pos_offset);
      if (ret == 2 && !strcmp(i1, inner_iter)) {
        max_pos_offset =
            i1_pos_offset > max_pos_offset ? i1_pos_offset : max_pos_offset;
      }

      ret = sscanf(ci2, "%[a-zA-z0-9]-%d", i2, &i2_neg_offset);
      if (ret == 2 && !strcmp(i2, inner_iter)) {
        max_neg_offset =
            i2_neg_offset > max_neg_offset ? i2_neg_offset : max_neg_offset;
      }
      ret = sscanf(ci2, "%[a-zA-z0-9]+%d", i2, &i2_pos_offset);
      if (ret == 2 && !strcmp(i2, inner_iter)) {
        max_pos_offset =
            i2_pos_offset > max_pos_offset ? i2_pos_offset : max_pos_offset;
      }

      sprintf(cur_decl + strlen(cur_decl), "%s_%s_%d", i1, i2, stmt->id);
    }

    int count = 0;
    if (!strcmp(i2, inner_iter)) {
      if (i2_neg_offset)
        count =
            sprintf(&(new_stmt_text[dest_index]), "_%s[%s_mod + %s - lbv - %d]",
                    cur_decl, cur_decl, inner_iter, i2_neg_offset);
      else if (i2_pos_offset)
        count =
            sprintf(&(new_stmt_text[dest_index]), "_%s[%s_mod + %s - lbv + %d]",
                    cur_decl, cur_decl, inner_iter, i2_pos_offset);
      else
        count = sprintf(&(new_stmt_text[dest_index]), "_%s[%s_mod + %s - lbv]",
                        cur_decl, cur_decl, inner_iter);
    } else if (arr->dim_orig == 2 && !strcmp(i1, inner_iter)) {
      if (i1_neg_offset)
        count = sprintf(&(new_stmt_text[dest_index]),
                        "_%s[%s_mod + (%s - lbv - %d) * %d]", cur_decl,
                        cur_decl, inner_iter, i1_neg_offset, t4);
      else if (i1_pos_offset)
        count = sprintf(&(new_stmt_text[dest_index]),
                        "_%s[%s_mod + (%s - lbv + %d) * %d]", cur_decl,
                        cur_decl, inner_iter, i1_pos_offset, t4);
      else
        count = sprintf(&(new_stmt_text[dest_index]),
                        "_%s[%s_mod + (%s - lbv) * %d]", cur_decl, cur_decl,
                        inner_iter, t4);
    } else
      count = sprintf(&(new_stmt_text[dest_index]), "_%s[%s_mod]", cur_decl,
                      cur_decl);

    dest_index += count;

    int found = 0;
    int i;
    for (i = 0; i < ind_size; i++) {

      if (strcmp(list_indices[i], cur_decl) == 0) {
        found = 1;
        break;
      }
    }

    if (!found) {
      list_indices[ind_size++] = strdup(cur_decl);

      if (!strcmp(i2, inner_iter) || !strcmp(i1, inner_iter)) {
        char s[100] = "";
        sprintf(s, " %s ", inner_iter);
        char *lbv_div = replacestr(div_ind, s, " lbv ");
        char *lbv_mod = replacestr(mod_ind, s, " lbv ");

        sprintf(decls + strlen(decls), "%s_div = %s;", cur_decl, lbv_div);
        sprintf(decls + strlen(decls), "%s_mod = %s;", cur_decl, lbv_mod);

        sprintf(incr + strlen(incr), "%s_mod++;", cur_decl);
      } else {

        sprintf(decls + strlen(decls), "%s_div = %s;", cur_decl, div_ind);
        sprintf(decls + strlen(decls), "%s_mod = %s;", cur_decl, mod_ind);
      }

      sprintf(prog->vec_decls + strlen(prog->vec_decls), "int %s_div = 0;",
              cur_decl);
      sprintf(prog->vec_decls + strlen(prog->vec_decls), "int %s_mod = 0;",
              cur_decl);
    }

    found = 0;
    char arr_decl[100] = "";
    sprintf(arr_decl, "%s_%s", arr->text, cur_decl);
    for (i = 0; i < ind_size; i++) {

      if (strcmp(list_indices[i], arr_decl) == 0) {
        found = 1;
        break;
      }
    }

    if (!found) {
      list_indices[ind_size++] = strdup(arr_decl);

      sprintf(decls + strlen(decls), "%s = %s_trans[%s_div];", arr_decl,
              arr->text, cur_decl);
      sprintf(prog->vec_decls + strlen(prog->vec_decls), "%s *%s = NULL;",
              arr->data_type, arr_decl);
    }

    //		src_index += strlen(div_ind) + strlen(mod_ind) + 2;
  }

  prog->num_decl_indices = ind_size;
  //	stmt->pos_peel_offset = max_pos_offset;
  //	stmt->neg_peel_offset = max_neg_offset;

  return 1;
}

char *pluto_dist_modify_stmt_text(char *stmt_text, int use_strides,
                                  PlutoProg *prog) {
  int j, k;

  char *new_stmt_text =
      (char *)malloc((strlen(stmt_text) + 1) * 500 * sizeof(char));
  new_stmt_text[0] = 0;

  char *array_name;
  int src_index = 0, dest_index = 0;
  int first, trans;

  while (1) {

    pluto_dist_copy_till_array(stmt_text, new_stmt_text, &src_index,
                               &dest_index);
    if (src_index >= strlen(stmt_text))
      break;

    array_name = pluto_dist_get_array_string(stmt_text, src_index);

    Array *arr = pluto_get_corrs_array(array_name, prog);
    /* If array is not found,just copy the stmt text as it is
    */

    if (arr == NULL) {
      new_stmt_text[dest_index++] = stmt_text[src_index++];
      continue;
    }

    strcat(new_stmt_text, "_trans");

    dest_index += strlen("_trans");

    char iter[arr->dim_orig + arr->npar + 1][100];
    for (j = 0; j < arr->dim_orig; ++j) {
      sscanf(stmt_text + src_index, "[%[^]]]", iter[j]);
      src_index += strlen(iter[j]) + 2;
    }

    for (j = 0; j < arr->npar; ++j) {
      strcpy(iter[j + arr->dim_orig], prog->params[j]);
    }

    strcpy(iter[arr->dim_orig + arr->npar], "1");

    //		assert(arr->trans->nrows == arr->dim);

    // change access according to the transfer function
    if (use_strides) {
      sprintf(new_stmt_text + strlen(new_stmt_text), "[");
      for (j = 0; j < arr->dim_orig; ++j) {
        if (arr->tiled_hyperplane[j].is_tiled) {

          sprintf(new_stmt_text + strlen(new_stmt_text), "(");
          first = 1;
          for (k = 0; k < arr->dim_orig + arr->npar + 1; ++k) {
            trans = arr->trans->val[arr->num_tiled_loops + j]
                                   [arr->num_tiled_loops + k];
            if (trans != 0) {
              if (!first)
                sprintf(new_stmt_text + strlen(new_stmt_text), " + ");
              sprintf(new_stmt_text + strlen(new_stmt_text), "(%d)*( %s )",
                      trans, iter[k]);

              first = 0;
            }
          }
          sprintf(new_stmt_text + strlen(new_stmt_text), ")/%d",
                  arr->tiled_hyperplane[j].tile_sizes);
          if (j != arr->dim_orig - 1)
            sprintf(new_stmt_text + strlen(new_stmt_text), "* %s + ",
                    pluto_dist_get_ptr_dim_stride(arr, j, prog));
        } else
          sprintf(new_stmt_text + strlen(new_stmt_text), "[ %s ]", iter[j]);
      }
      sprintf(new_stmt_text + strlen(new_stmt_text), "]");

    } else {
      for (j = 0; j < arr->dim_orig; ++j) {
        if (arr->tiled_hyperplane[j].is_tiled) {

          sprintf(new_stmt_text + strlen(new_stmt_text), "[(");
          int first = 1;
          for (k = 0; k < arr->dim_orig + arr->npar + 1; ++k) {
            int trans = arr->trans->val[arr->num_tiled_loops + j]
                                       [arr->num_tiled_loops + k];
            if (trans != 0) {
              if (!first)
                sprintf(new_stmt_text + strlen(new_stmt_text), " + ");
              sprintf(new_stmt_text + strlen(new_stmt_text), "(%d)*( %s )",
                      trans, iter[k]);

              first = 0;
            }
          }
          sprintf(new_stmt_text + strlen(new_stmt_text), ")/%d]",
                  arr->tiled_hyperplane[j].tile_sizes);
        }
      }
    }

    sprintf(new_stmt_text + strlen(new_stmt_text), "[");
    for (j = 0; j < arr->dim_orig; ++j) {
      if (arr->tiled_hyperplane[j].is_tiled) {

        sprintf(new_stmt_text + strlen(new_stmt_text), "(");
        first = 1;
        for (k = 0; k < arr->dim_orig + arr->npar + 1; ++k) {
          trans = arr->trans->val[arr->num_tiled_loops + j]
                                 [arr->num_tiled_loops + k];
          if (trans != 0) {
            if (!first)
              sprintf(new_stmt_text + strlen(new_stmt_text), " + ");
            sprintf(new_stmt_text + strlen(new_stmt_text), "(%d)*( %s )", trans,
                    iter[k]);

            first = 0;
          }
        }
        sprintf(new_stmt_text + strlen(new_stmt_text), ")%%%d",
                arr->tiled_hyperplane[j].tile_sizes);
        if (j != arr->dim_orig - 1)
          sprintf(new_stmt_text + strlen(new_stmt_text), "* %s + ",
                  pluto_dist_get_tile_dim_stride(arr, j, prog));
      } else
        sprintf(new_stmt_text + strlen(new_stmt_text), "[%s]", iter[j]);
    }
    sprintf(new_stmt_text + strlen(new_stmt_text), "]");

    dest_index = strlen(new_stmt_text);

    free(array_name);
  }

  return new_stmt_text;
}

/* Generate variable declarations and macros */
int generate_declarations(PlutoProg *prog, FILE *outfp) {
  int i, j;

  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  prog->num_decl_indices = 0;
  prog->decl_indices = (char **)malloc(1024 * sizeof(char *));
  prog->vec_decls = malloc(5 * 1024);
  prog->vec_decls[0] = 0;

  /* Generate statement macros */
  for (i = 0; i < nstmts; i++) {
    gen_stmt_macro(stmts[i], outfp, prog);
  }
  fprintf(outfp, "\n");
  fprintf(outfp, "%s", prog->vec_decls);

  free(prog->decl_indices);
  free(prog->vec_decls);
  prog->num_decl_indices = 0;

  char *suffix = (char *)malloc((options->scopnum + 1) * sizeof(char));
  suffix[0] = '\0';
  for (j = 0; j < options->scopnum; j++) {
    strcat(suffix, "_");
  }

  /* Insert "dynsched" specific code
   * Number of parameters to task_computation is fixed as 4 in the static code -
   * scheduler.h
   * So hardcode the parameters to the function as t1, t2, t3, t4
   */
  if (options->dynschedule_graph_old) {
    fprintf(outfp,
            "\n\nvoid task_computation(int t1, int t2, int t3, int t4)\n{\n\n");
    fprintf(outfp, "\tint ");
    for (i = 4; i < stmts[0]->trans->nrows; i++) {
      if (i != 4)
        fprintf(outfp, ", ");
      fprintf(outfp, "t%d", i + 1);
    }
    fprintf(outfp, ";\n\n");
  } else {
    /* Scattering iterators. */
    if (prog->num_hyperplanes >= 1) {
      fprintf(outfp, "\t\tint ");
      int start = (prog->num_parameterized_loops > 0)
                      ? prog->num_parameterized_loops
                      : 0;
      for (i = start; i < prog->num_hyperplanes; i++) {
        if (i != start)
          fprintf(outfp, ", ");
        fprintf(outfp, "t%d%s", i + 1, suffix);
        if (prog->hProps[i].unroll) {
          fprintf(outfp, ", t%dt%s, newlb_t%d%s, newub_t%d%s", i + 1, suffix,
                  i + 1, suffix, i + 1, suffix);
        }
      }
      fprintf(outfp, ";\n\n");
    }
  }

  if (options->parallel || options->distmem) {
    fprintf(outfp, "\t\tint lb%s, ub%s, lbd%s, ubd%s, lb2%s, ub2%s;\n", suffix,
            suffix, suffix, suffix, suffix, suffix);
  } else if (options->dynschedule) {
    fprintf(outfp, "\tint lbp, ubp;\n");
  }

  /* For vectorizable loop bound replacement */
  fprintf(outfp, "\t\tregister int lbv%s, ubv%s;\n\n", suffix, suffix);

  free(suffix);

  return 0;
}

void pluto_mark_vec_stmts(FILE *cloogfp, PlutoProg *prog) {

  CloogInput *input;
  CloogOptions *cloogOptions;
  CloogState *state;
  struct clast_stmt *root;
  int nstmts = prog->nstmts;
  state = cloog_state_malloc();
  cloogOptions = cloog_options_malloc(state);

  cloogOptions->fs = malloc(nstmts * sizeof(int));
  cloogOptions->ls = malloc(nstmts * sizeof(int));
  cloogOptions->fs_ls_size = nstmts;

  int i;
  for (i = 0; i < nstmts; i++) {
    cloogOptions->fs[i] = -1;
    cloogOptions->ls[i] = -1;
  }

  input = cloog_input_read(cloogfp, cloogOptions);

  IF_DEBUG(printf("[Pluto] cloog_clast_create\n"));
  root = cloog_clast_create_from_input(input, cloogOptions);
  if (options->prevector) {
    pluto_mark_vector(root, prog, cloogOptions);
  }

  cloog_clast_free(root);
  cloog_options_free(cloogOptions);
  cloog_state_free(state);
  rewind(cloogfp);
}
/* Call cloog and generate code for the transformed program
 *
 * cloogf, cloogl: set to -1 if you want the function to decide
 *
 * --cloogf, --cloogl overrides everything; next cloogf, cloogl if != -1,
 *  then the function takes care of the rest
 */
int pluto_gen_cloog_code(const PlutoProg *prog, int cloogf, int cloogl,
                         FILE *cloogfp, FILE *outfp) {
  CloogInput *input;
  CloogOptions *cloogOptions;
  CloogState *state;
  int i;

  struct clast_stmt *root;

  Stmt **stmts = prog->stmts;
  int nstmts = prog->nstmts;

  state = cloog_state_malloc();
  cloogOptions = cloog_options_malloc(state);

  cloogOptions->fs = (int *)malloc(nstmts * sizeof(int));
  cloogOptions->ls = (int *)malloc(nstmts * sizeof(int));
  cloogOptions->fs_ls_size = nstmts;

  for (i = 0; i < nstmts; i++) {
    cloogOptions->fs[i] = -1;
    cloogOptions->ls[i] = -1;
  }

  cloogOptions->name = "CLooG file produced by PLUTO";
  cloogOptions->compilable = 0;
  cloogOptions->esp = 1;
  cloogOptions->strides = 1;
  cloogOptions->quiet = !options->debug;

  /* Generates better code in general */
  cloogOptions->backtrack = options->cloogbacktrack;

  if (options->cloogf >= 1 && options->cloogl >= 1) {
    cloogOptions->f = options->cloogf;
    cloogOptions->l = options->cloogl;
  } else {
    if (cloogf >= 1 && cloogl >= 1) {
      cloogOptions->f = cloogf;
      cloogOptions->l = cloogl;
    } else if (options->dynschedule_graph || options->dynschedule ||
               options->distmem) {
      /* We will set f/l statement-wise */
      for (i = 0; i < nstmts; i++) {
        if ((stmts[i]->type == ORIG) || (stmts[i]->type == ORIG_IN_FUNCTION)) {
          cloogOptions->fs[i] = get_first_point_loop(stmts[i], prog) + 1;
          cloogOptions->ls[i] = prog->num_hyperplanes;
        } else if (stmts[i]->type == IN_FUNCTION) {
          cloogOptions->fs[i] = 1;
          cloogOptions->ls[i] = prog->num_hyperplanes;
        } else {
          cloogOptions->fs[i] = prog->num_hyperplanes - 1;
          cloogOptions->ls[i] = prog->num_hyperplanes;
        }
      }
    } else if (options->tile) {
      for (i = 0; i < nstmts; i++) {
        cloogOptions->fs[i] = get_first_point_loop(stmts[i], prog) + 1;
        cloogOptions->ls[i] = prog->num_hyperplanes;
      }
    } else {
      /* Default */
      cloogOptions->f = 1;
      /* last level to optimize: number of hyperplanes;
       * since Pluto provides full-ranked transformations */
      cloogOptions->l = prog->num_hyperplanes;
    }
  }

  if (options->data_dist && options->data_tile_opt)
    cloogOptions->invariant_decl = 1;

  if (!options->silent) {
    if (nstmts >= 1 && cloogOptions->fs[0] >= 1) {
      printf("[pluto] using statement-wise -fs/-ls options: ");
      for (i = 0; i < nstmts; i++) {
        printf("S%d(%d,%d), ", i + 1, cloogOptions->fs[i], cloogOptions->ls[i]);
      }
      printf("\n");
    } else {
      printf("[pluto] using Cloog -f/-l options: %d %d\n", cloogOptions->f,
             cloogOptions->l);
    }
  }

  if (options->cloogsh)
    cloogOptions->sh = 1;

  cloogOptions->name = "PLUTO-produced CLooG file";

  fprintf(outfp, "/* Start of CLooG code */\n");
  /* Get the code from CLooG */
  IF_DEBUG(printf("[pluto] cloog_input_read\n"));
  input = cloog_input_read(cloogfp, cloogOptions);
  IF_DEBUG(printf("[pluto] cloog_clast_create\n"));
  root = cloog_clast_create_from_input(input, cloogOptions);
  if (options->prevector) {
    pluto_mark_vector(root, prog, cloogOptions);
  }
  if (prog->num_parameterized_loops != -1) {
    // statement ID doesn't matter since 'root' isn't printed
    char iter[5];
    sprintf(iter, "t%d", prog->num_parameterized_loops);
    int stmtids[prog->nstmts];
    for (i = 0; i < prog->nstmts; i++) {
      stmtids[i] = prog->stmts[i]->id + 1;
    }
    struct clast_for **parameterized_loops;
    int *stmts;
    int nloops, nstmts;
    ClastFilter filter = { iter, stmtids, prog->nstmts, exact };
    clast_filter(root, filter, &parameterized_loops, &nloops, &stmts, &nstmts);
    assert(nloops == 1 && "[Pluto] parallel poly loop not found in AST\n");
    clast_pprint(outfp, parameterized_loops[0]->body, 0, cloogOptions);
    free(parameterized_loops);
    free(stmts);
  } else {
    if (options->dynschedule) {
      pluto_mark_parallel_dynschedule(root, prog, cloogOptions);
    } else if (options->parallel || options->distmem) {
      pluto_mark_parallel(root, prog, cloogOptions);
    }
    clast_pprint(outfp, root, 0, cloogOptions);
  }
  cloog_clast_free(root);

  fprintf(outfp, "/* End of CLooG code */\n");

  if (options->dynschedule_graph_old) {
    fprintf(outfp, "\n\n}\n\n");
  }

  cloog_options_free(cloogOptions);
  cloog_state_free(state);

  return 0;
}

/* Generate code for a single multicore; the ploog script will insert openmp
 * pragmas later */
int pluto_multicore_codegen(FILE *cloogfp, FILE *outfp, PlutoProg *prog) {
  if (options->parallel) {
    fprintf(outfp, "#include <omp.h>\n\n");
  }

  // pluto_mark_vec_stmts(cloogfp, prog);
  generate_declarations(prog, outfp);

  if (options->multipar) {
    fprintf(outfp, "\tomp_set_nested(1);\n");
    fprintf(outfp, "\tomp_set_num_threads(2);\n");
  }

  pluto_gen_cloog_code(prog, -1, -1, cloogfp, outfp);

  return 0;
}

/* Decides which loops to mark parallel and generates the corresponding OpenMP
 * pragmas and writes them out to a file. They are later read by a script
 * (ploog) and appropriately inserted into the output Cloog code
 *
 * Returns: the number of parallel loops for which OpenMP pragmas were generated
 *
 * Generate the #pragma comment -- will be used by a syntactic scanner
 * to put in place -- should implement this with CLast in future
 **/
int pluto_omp_parallelize(PlutoProg *prog) {
  int i, j;

  FILE *outfp = fopen(".pragmas", "w");

  if (!outfp)
    return 1;

  HyperplaneProperties *hProps = prog->hProps;

  int loop;

  char *suffix = malloc((options->scopnum + 1) * sizeof(char));
  suffix[0] = '\0';
  for (j = 0; j < options->scopnum; j++) {
    strcat(suffix, "_");
  }

  /* IMPORTANT: Note that by the time this function is called, pipelined
   * parallelism has already been converted to inner parallelism in
   * tile space (due to a tile schedule) - so we don't need check any
   * PIPE_PARALLEL properties
   */
  /* Detect the outermost sync-free parallel loop - find upto two of them if
   * the multipar option is set */
  int num_parallel_loops = 0;
  for (loop = 0; loop < prog->num_hyperplanes; loop++) {
    if (hProps[loop].dep_prop == PARALLEL && hProps[loop].type != H_SCALAR) {
      // Remember our loops are 1-indexed (t1, t2, ...)
      fprintf(outfp, "t%d%s #pragma omp parallel for shared(", loop + 1,
              suffix);

      for (i = 0; i < loop; i++) {
        fprintf(outfp, "t%d%s,", i + 1, suffix);
      }

      for (i = 0; i < num_parallel_loops + 1; i++) {
        if (i != 0)
          fprintf(outfp, ",");
        fprintf(outfp, "lb%d%s,ub%d%s", i + 1, suffix, i + 1, suffix);
      }

      fprintf(outfp, ") private(");

      if (options->prevector) {
        fprintf(outfp, "ubv%s,lbv%s,", suffix, suffix);
      }

      /* Lower and upper scalars for parallel loops yet to be marked */
      /* NOTE: we extract up to 2 degrees of parallelism
       */
      if (options->multipar) {
        for (i = num_parallel_loops + 1; i < 2; i++) {
          fprintf(outfp, "lb%d%s,ub%d%s,", i + 1, suffix, i + 1, suffix);
        }
      }

      for (i = loop; i < prog->num_hyperplanes; i++) {
        if (i != loop)
          fprintf(outfp, ",");
        fprintf(outfp, "t%d%s", i + 1, suffix);
      }
      fprintf(outfp, ")\n");

      num_parallel_loops++;

      if (!options->multipar || num_parallel_loops == 2) {
        break;
      }
    }
  }

  IF_DEBUG(fprintf(stdout, "[pluto] marked %d loop(s) parallel\n",
                   num_parallel_loops));

  free(suffix);

  fclose(outfp);

  return num_parallel_loops;
}

// s should be root when called
// !!! only the first occurence of a set of statements at loop level n is split
// !!! what if there are multiple occurrences of sets of statements at loop
// level n? e.g. when the loop is split
static int splitLoops(struct clast_stmt *s, int loop_level, int n,
                      struct clast_stmt **split, int stmt_num,
                      CloogState *state, CloogDomain *domain) {
  struct clast_user_stmt *u = NULL;
  CloogStatement *cs = NULL;

  for (; s; s = s->next) {
    if (CLAST_STMT_IS_A(s, stmt_for)) {
      struct clast_for *for_s = (struct clast_for *)s;
      if (loop_level + 1 == n) {
        cs = cloog_statement_alloc(state, stmt_num);
        u = new_clast_user_stmt(domain, cs,
                                NULL); // add a user defined statement
        cloog_statement_free(cs);
        *split = for_s->body;
        for_s->body = &u->stmt;
        break;
      }
      splitLoops(for_s->body, loop_level + 1, n, split, stmt_num, state,
                 domain);
    } else if (CLAST_STMT_IS_A(s, stmt_guard)) {
      splitLoops(((struct clast_guard *)s)->then, loop_level, n, split,
                 stmt_num, state, domain);
      if (*split != NULL)
        break;
    } else if (CLAST_STMT_IS_A(s, stmt_block)) {
      splitLoops(((struct clast_block *)s)->body, loop_level, n, split,
                 stmt_num, state, domain);
      if (*split != NULL)
        break;
    }
  } // end of while
  return 0;
}

void print_dynsched_file(char *srcFileName, FILE *cloogfp, FILE *outfp,
                         PlutoProg *prog) {
  int i, j, k;
  int srcid, destid, startIteratorIndex, noOfIterators, common_band_num,
      stmt_common_tile_cnt;
  char sysloogFileName[256];
  char tempString[256], variablesString[256], argumentsString[256];
  FILE *sysloogfp;
  PlutoConstraints *transformedDpolytope;
  Stmt *srcStmt, *destStmt;

  CloogInput *input;
  CloogOptions *cloogOptions;
  CloogState *state;
  CloogDomain *domain;
  struct clast_stmt *root, *split;
  split = NULL;

  Stmt **stmts = prog->stmts;
  Dep **deps = prog->deps;
  int nstmts = prog->nstmts;
  int ndeps = prog->ndeps;
  int npar = prog->npar;
  char **params = prog->params;
  HyperplaneProperties *hProps = prog->hProps;
  int num_hyperplanes = prog->num_hyperplanes;

  /* Number of common tiles is taken as number of common iterators (after
   * transformation)
   *  in the first non-sequential band.
   *  However the clast functions have assumed that the common tiles start from
   * the first band */
  int num_common_tiles = prog->nvar;
  for (j = 0; j < num_hyperplanes; j++) {
    if (hProps[j].dep_prop != SEQ)
      break;
  }
  common_band_num = hProps[j].band_num;
  for (i = 0; i < nstmts; i++) {
    stmt_common_tile_cnt = 0;
    j = 0;
    while ((hProps[j].band_num != common_band_num) && (j < num_hyperplanes))
      j++;
    while ((hProps[j].band_num == common_band_num) && (j < num_hyperplanes)) {
      j++;
      stmt_common_tile_cnt++;
    }
    if (stmt_common_tile_cnt < num_common_tiles)
      num_common_tiles = stmt_common_tile_cnt;
  }
  IF_DEBUG(printf("No of variables = %d , No of common tiles = %d\n",
                  prog->nvar, num_common_tiles));

  state = cloog_state_malloc();
  cloogOptions = cloog_options_malloc(state);

  cloogOptions->name = "CLooG file produced by PLUTO for dynamic scheduling";
  cloogOptions->compilable = 0;
  cloogOptions->esp =
      0; // !!!roshan not sure if this needs to be enforced - can 1 be used?
  cloogOptions->quiet = options->silent;
  cloogOptions->backtrack = 1; /* Generates better code in general */

  /* Generate a file that has a function to create DAG
   *
   * Generate tiled code to split the common tiled loops
   * to enumerate all vertices
   *
   * The loop is written twice, once to estimate the number of vertices
   * and once to actually add the vertices
   */

  /* Assume that the common tile iterators would be named 't1, t2...'*/
  /* Temporarily dag_add_vertex and dag_add_edge assumes 4-tuple vertex in
   * scheduler.h */
  // build the string for intra-tile iterators, assuming the character 't' for
  // it
  strcpy(tempString, "");
  for (j = 1; j < num_common_tiles; j++)
    sprintf(tempString, "%s t%d,", tempString, j);
  sprintf(tempString, "%s t%d", tempString, j);
  strcpy(variablesString, tempString);
  strcat(variablesString, ", ");
  // append 0 for the arguments to satisfy the static 4-tuple vertex in
  // scheduler.h
  /* Careful - temporary stuff*/
  for (j = num_common_tiles + 1; j <= 4; j++)
    strcat(tempString, ", 0");
  strcpy(argumentsString, tempString);
  strcat(argumentsString, ", ");

  fprintf(outfp, "#define S%d() {dag_add_vertex(%s, 1, &vId);}\n", 3,
          tempString);

  // build the string for inter-tile iterators, assuming the character 's' for
  // it
  strcpy(tempString, "");
  for (j = 1; j < num_common_tiles; j++)
    sprintf(tempString, "%s s%d,", tempString, j);
  sprintf(tempString, "%s s%d", tempString, j);
  strcat(variablesString, tempString);
  // append 0 for the arguments to satisfy the static 4-tuple vertex in
  // scheduler.h
  /* Careful - temporary stuff*/
  for (j = num_common_tiles + 1; j <= 4; j++)
    strcat(tempString, ", 0");
  strcat(argumentsString, tempString);

  // add_edge statement has ID 1
  // if there is an issue due to this ID overlapping with that of the original
  // code statement,
  // then it can be changed to nstmts+3
  fprintf(outfp, "#define S%d(%s) {dag_add_edge(%s, 0.0);}\n", 1,
          variablesString, argumentsString);
  fprintf(outfp, "#define S%d() {++vertexCount;}\n", 2);
  fprintf(outfp, "\nvoid generate_dag()\n{\n");
  fprintf(outfp, "  int %s, vId;\n", variablesString);
  fprintf(outfp, "  int vertexCount=0;\n\n");

  if (options->ft == -1) {
    if (stmts[0]->num_tiled_loops < 4) {
      cloogOptions->f = stmts[0]->num_tiled_loops + 1;
      cloogOptions->l = stmts[0]->trans->nrows;
    } else {
      cloogOptions->f = stmts[0]->num_tiled_loops + 1;
      cloogOptions->l = stmts[0]->trans->nrows;
    }
  } else {
    cloogOptions->f = stmts[0]->num_tiled_loops + options->ft + 1;
    cloogOptions->l = stmts[0]->trans->nrows;
  }

  input = cloog_input_read(cloogfp, cloogOptions);
  domain = cloog_domain_copy(input->context);
  root = cloog_clast_create_from_input(input, cloogOptions);

  // generate code to estimate the number of vertices
  splitLoops(root, 0, num_common_tiles, &split, 2, state, domain);
  clast_pprint(outfp, root, 2, cloogOptions);

  assert(split != NULL);
  cloog_clast_free(split);
  split = NULL;

  // initialize the DAG with the estimated number of vertices
  fprintf(outfp, "\n  init_dag(vertexCount);\n\n");

  // generate code to add the vertices
  splitLoops(root, 0, num_common_tiles, &split, 3, state, domain);
  clast_pprint(outfp, root, 2, cloogOptions);

  assert(split != NULL);
  cloog_clast_free(split);
  cloog_clast_free(root);

  fprintf(outfp, "\n");

  /* This is the core part
   *  Generate inter-tile dependencies (cloog file)
   *  Create a system to capture the dependencies in the transformed tiled
   * space:
   *          - source equalities/inequaities
   *              1. original domain + tile definition (domain in tiled.cloog)
   *              2. transformations (scattering functions in tiled.cloog)
   *          - destination equalities/inequaities
   *              1. original domain + tile definition (domain in tiled.cloog)
   *              2. transformations (scattering functions in tiled.cloog)
   *          - equalities of original h-transformation of dependence
   */
  for (i = 0; i < ndeps; i++) {
    /* RAW, WAR, WAW deps matter */
    /*if(!IS_RAR(deps[i]->type)) */
    // !!!roshan not sure if WAR and WAW dependences need to be considered.
    // If they are considered, it leads to slower execution time sometimes (due
    // to increase in size of the DAG?),
    // longer compilation time for generated code of some benchmarks like
    // heat-3d
    // and incorrect results for some benchmarks like heat-2d
    if (deps[i]->type == OSL_DEPENDENCE_RAW) {
      strcpy(sysloogFileName, srcFileName);
      sysloogFileName[strlen(srcFileName) - 2] = '\0';
      sprintf(sysloogFileName, "%s.dep%d.dynsched.sysloog", sysloogFileName, i);

      // Generate the cloog file to generate dag code corresponding to the dep
      sysloogfp = fopen(sysloogFileName, "w+");
      if (!sysloogfp) {
        fprintf(stderr, "Can't open file %s for writing\n", sysloogFileName);
        return;
      }

      fprintf(sysloogfp, "# language: C\n");
      fprintf(sysloogfp, "c\n\n");

      fprintf(sysloogfp, "# Context\n");
      fprintf(sysloogfp, "%d %d\n", npar, npar + 2);

      for (j = 0; j < npar; j++) {
        fprintf(sysloogfp, "1 ");
        for (k = 0; k < npar; k++) {
          if (j == k) {
            fprintf(sysloogfp, "1 ");
          } else
            fprintf(sysloogfp, "0 ");
        }
        fprintf(sysloogfp, "%d\n", -options->codegen_context);
      }

      fprintf(sysloogfp, "\n1\n");
      fprintf(sysloogfp, "# Parameter name(s)\n");
      for (j = 0; j < npar; j++) {
        fprintf(sysloogfp, "%s ", params[j]);
      }
      fprintf(sysloogfp, "\n\n");

      fprintf(sysloogfp, "# Number of statements\n");
      fprintf(sysloogfp, "1\n\n");

      fprintf(sysloogfp, "# Iteration Domain\n");
      fprintf(sysloogfp, "1\n\n");

      srcid = deps[i]->src;
      destid = deps[i]->dest;
      srcStmt = stmts[srcid];
      destStmt = stmts[destid];

      // get the transformed dependence polytope, which contains the scattering
      // functions
      transformedDpolytope =
          pluto_get_transformed_dpoly(deps[i], srcStmt, destStmt);

      assert(transformedDpolytope->ncols == (2 * num_hyperplanes + npar + 1));

      // project out the destination statement intra-tile iterators
      startIteratorIndex = num_hyperplanes + num_common_tiles;
      noOfIterators = 2 * num_hyperplanes - startIteratorIndex;
      pluto_constraints_project_out(transformedDpolytope, startIteratorIndex,
                                    noOfIterators);

      // project out the source statement intra-tile iterators
      startIteratorIndex = num_common_tiles;
      noOfIterators = num_hyperplanes - startIteratorIndex;
      pluto_constraints_project_out(transformedDpolytope, startIteratorIndex,
                                    noOfIterators);

      pluto_constraints_print_polylib(sysloogfp, transformedDpolytope);

      fprintf(sysloogfp, "\n0 0 0\n\n");
      fprintf(sysloogfp, "1\n");
      fprintf(sysloogfp, "# Iterator name(s)\n\n");
      // inter-tile iterators use the character 's', according to assumption
      // above while generating the variable declarations
      for (j = 1; j <= num_common_tiles; j++)
        fprintf(sysloogfp, "s%d ", j);
      // intra-tile iterators use the character 't', according to assumption
      // above while generating the variable declarations
      for (j = 1; j <= num_common_tiles; j++)
        fprintf(sysloogfp, "t%d ", j);
      fprintf(sysloogfp, "\n\n");
      fprintf(sysloogfp, "# Scattering functions\n");
      fprintf(sysloogfp, "0\n\n");

      rewind(sysloogfp);

      cloogOptions->f = 1;
      cloogOptions->l = num_common_tiles;

      input = cloog_input_read(sysloogfp, cloogOptions);
      root = cloog_clast_create_from_input(input, cloogOptions);
      // assuming only one (innermost) statement, whose ID will be 1
      // if it is necessary to change the ID to nstmts+3,
      // the clast statements should be traversed to find the only clast user
      // statement
      // assuming statement ID of 1 will work, even though the ID overlaps with
      // one of the original code statements
      clast_pprint(outfp, root, 2, cloogOptions);
      cloog_clast_free(root);
      fclose(sysloogfp);
    }
  }

  fprintf(outfp, "\n  update_number_vertices();");
  fprintf(outfp, "\n\n}\n\n");
  fprintf(outfp, "#undef S1\n");
  fprintf(outfp, "#undef S2\n");
  fprintf(outfp, "#undef S3\n\n\n");

  cloog_options_free(cloogOptions);

  prog->num_parameterized_loops = num_common_tiles;
}
