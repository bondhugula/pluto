/*
 * PLuTo: An automatic parallelier and locality optimizer
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
#ifndef _CONSTRAINTS_H
#define _CONSTRAINTS_H

PlutoInequalities *constraints_alloc(int nrows, int ncols);
void constraints_free(PlutoInequalities *cst);
void constraints_print(FILE *fp, PlutoInequalities *);
void constraints_pretty_print(FILE *fp, PlutoInequalities *);
PlutoInequalities *constraints_copy(PlutoInequalities *, PlutoInequalities *);

int best_elim_candidate(PlutoInequalities *, int);
void fourier_motzkin_eliminate(PlutoInequalities *, int n);

PlutoInequalities *pluto2pip(PlutoInequalities *, PlutoInequalities *pipcst);

PlutoInequalities *constraints_add(PlutoInequalities *, PlutoInequalities *);
void constraints_simplify(PlutoInequalities *const cst);

int *constraints_solve(PlutoInequalities *);

#endif
