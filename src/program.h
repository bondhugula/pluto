/*
 * PLUTO: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
#include "clan/clan.h"
#include "candl/candl.h"
#include "constraints.h"

Stmt *pluto_stmt_alloc(int dim, const PlutoConstraints *domain, const PlutoMatrix *mat);
void pluto_stmt_free(Stmt *stmt);
void pluto_stmts_print(FILE *fp, Stmt **, int);
void pluto_stmt_print(FILE *fp, const Stmt *stmt);
Stmt *pluto_stmt_dup(const Stmt *src);

void pluto_prog_print(PlutoProg *prog);

Dep *pluto_dep_alloc();
void pluto_dep_print(FILE *fp, Dep *dep);
void pluto_deps_print(FILE *, Dep **, int);

PlutoProg *pluto_prog_alloc();
void pluto_prog_free(PlutoProg *prog);
PlutoProg *scop_to_pluto_prog(scoplib_scop_p scop, PlutoOptions *options);

int get_coeff_upper_bound(PlutoProg *prog);

void pluto_prog_add_param(PlutoProg *prog, const char *param, int pos);
void pluto_add_stmt(PlutoProg *prog, 
        const PlutoConstraints *domain,
        const PlutoMatrix *trans,
        char ** iterators,
        const char *text,
        PlutoStmtType type
        );

void pluto_stmt_add_dim(Stmt *stmt, int pos, int time_pos, const char *iter, 
        PlutoHypType type, PlutoProg *prog);
void pluto_stmt_remove_dim(Stmt *stmt, int pos, PlutoProg *prog);
void pluto_prog_add_hyperplane(PlutoProg *prog, int pos, PlutoHypType type);

void pluto_separate_stmt(PlutoProg *prog, const Stmt *stmt, int level);
void pluto_separate_stmts(PlutoProg *prog, Stmt **stmts, int num, int level);

int extract_stmts(__isl_keep isl_union_set *domains, Stmt **stmts);
int pluto_is_hyperplane_scalar(const Stmt *stmt, int level);
int pluto_stmt_is_member_of(Stmt *s, Stmt **slist, int len);
PlutoAccess **pluto_get_all_waccs(PlutoProg *prog, int *num);
int pluto_stmt_is_subset_of(Stmt **s1, int n1, Stmt **s2, int n2);
void pluto_stmt_add_hyperplane(Stmt *stmt, PlutoHypType type, int pos);
PlutoMatrix *pluto_get_new_access_func(const Stmt *stmt, const PlutoMatrix *acc);

int extract_deps(Dep **deps, int first, Stmt **stmts,
        __isl_keep isl_union_map *umap, int type);

#endif
