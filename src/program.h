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

void stmts_print(FILE *fp, Stmt **, int);
void stmt_free(Stmt *stmt);
Stmt *stmt_copy (Stmt *src);

void deps_print (FILE *, Dep *, int);

PlutoProg *scop_to_pluto_prog(scoplib_scop_p scop, PlutoOptions *options);
void pluto_prog_free(PlutoProg *prog);

PlutoOptions *pluto_options_alloc();
void pluto_options_free(PlutoOptions *);
int get_coeff_upper_bound(PlutoProg *prog);

void pluto_add_parameter(PlutoProg *prog, char *param);
void pluto_stmt_add_dim(Stmt *stmt, int pos, int time_pos, char *iter);
void pluto_prog_add_hyperplane(PlutoProg *prog, int pos);

#endif
