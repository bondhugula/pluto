/*
 * PLUTO: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007 Uday Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
#ifndef _PROGRAM_H

#include "pluto.h"
#include "constraints.h"

#include "clan/clan.h"
#include "candl/candl.h"

#include "osl/scop.h"

Stmt *pluto_stmt_alloc(int dim, const PlutoConstraints *domain, const PlutoMatrix *mat);
void pluto_stmt_free(Stmt *stmt);
Stmt *pluto_stmt_dup(const Stmt *src);

void pluto_stmts_print(FILE *fp, Stmt **, int);
void pluto_stmt_print(FILE *fp, const Stmt *stmt);
void pluto_prog_print(FILE *fp, PlutoProg *prog);

Dep *pluto_dep_alloc();
void pluto_dep_print(FILE *fp, const Dep *dep);
void pluto_deps_print(FILE *, PlutoProg *prog);

PlutoProg *pluto_prog_alloc();
void pluto_prog_free(PlutoProg *prog);
PlutoProg *scop_to_pluto_prog(osl_scop_p scop, PlutoOptions *options);

int get_coeff_upper_bound(PlutoProg *prog);

void pluto_prog_add_param(PlutoProg *prog, const char *param, int pos);
void pluto_add_stmt(PlutoProg *prog, 
        const PlutoConstraints *domain,
        const PlutoMatrix *trans,
        char ** iterators,
        const char *text,
        PlutoStmtType type
        );

void pluto_add_stmt_to_end(PlutoProg *prog, 
        const PlutoConstraints *domain,
        char **iterators,
        const char *text,
        int level,
        PlutoStmtType type
        );

void pluto_stmt_add_dim(Stmt *stmt, int pos, int time_pos, const char *iter, 
        PlutoHypType type, PlutoProg *prog);
void pluto_stmt_remove_dim(Stmt *stmt, int pos, PlutoProg *prog);
void pluto_prog_add_hyperplane(PlutoProg *prog, int pos, PlutoHypType type);

int get_const_bound_difference(const PlutoConstraints *cst, int depth);
PlutoMatrix *get_alpha(const Stmt *stmt, const PlutoProg *prog);
PlutoMatrix *pluto_stmt_get_remapping(const Stmt *stmt, int **strides);

void get_parametric_extent(const PlutoConstraints *cst, int pos,
        int npar, const char **params, char **extent);

void get_parametric_extent_const(const PlutoConstraints *cst, int pos,
        int npar, const char **params, char **extent);

char *get_parametric_bounding_box(const PlutoConstraints *cst, int start, 
        int num, int npar, const char **params);

void pluto_separate_stmt(PlutoProg *prog, const Stmt *stmt, int level);
void pluto_separate_stmts(PlutoProg *prog, Stmt **stmts, int num, int level, int offset);

int pluto_is_hyperplane_scalar(const Stmt *stmt, int level);
int pluto_stmt_is_member_of(int stmt_id, Stmt **slist, int len);
PlutoAccess **pluto_get_all_waccs(const PlutoProg *prog, int *num);
int pluto_stmt_is_subset_of(Stmt **s1, int n1, Stmt **s2, int n2);
void pluto_stmt_add_hyperplane(Stmt *stmt, PlutoHypType type, int pos);
PlutoMatrix *pluto_get_new_access_func(const Stmt *stmt, const PlutoMatrix *acc, int **divs); 


int extract_deps(Dep **deps, int first, Stmt **stmts,
        __isl_keep isl_union_map *umap, int type);
int isl_map_count(__isl_take isl_map *map, void *user);

int pluto_get_max_ind_hyps(const PlutoProg *prog);
int pluto_get_max_ind_hyps_non_scalar(const PlutoProg *prog);
int pluto_stmt_get_num_ind_hyps(const Stmt *stmt);
int pluto_stmt_get_num_ind_hyps_non_scalar(const Stmt *stmt);
int pluto_transformations_full_ranked(PlutoProg *prog);
void pluto_pad_stmt_transformations(PlutoProg *prog);

void pluto_populate_scop (osl_scop_p scop, PlutoProg *prog,
                           PlutoOptions *options);

void pluto_access_print(FILE *fp, const PlutoAccess *acc, const Stmt *stmt);
void pluto_transformations_print(const PlutoProg *prog);
void pluto_transformations_pretty_print(const PlutoProg *prog);
void pluto_print_hyperplane_properties(const PlutoProg *prog);
void pluto_stmt_transformation_print(const Stmt *stmt);
void pluto_stmt_print_hyperplane(FILE *fp, const Stmt *stmt, int level);
void pluto_transformation_print_level(const PlutoProg *prog, int level);

Stmt *pluto_stmt_dup(const Stmt *stmt);
PlutoAccess *pluto_access_dup(const PlutoAccess *acc);
void pluto_dep_free(Dep *dep);
Dep *pluto_dep_dup(Dep *d);
void pluto_remove_stmt(PlutoProg *prog, int stmt_id);

int pluto_prog_get_largest_const_in_domains(const PlutoProg *prog);

#endif
