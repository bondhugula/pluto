/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2015 Uday Bondhugula
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
#ifndef _CONSTRAINTS_POLYLIB_H
#define _CONSTRAINTS_POLYLIB_H

#include "polylib/polylib64.h"

typedef struct pluto_constraints PlutoConstraints;
typedef struct plutoContext PlutoContext;
typedef struct pluto_matrix PlutoMatrix;

Polyhedron *pluto_constraints_to_polylib(const PlutoConstraints *cst);
PlutoConstraints *polylib_to_pluto_constraints(Polyhedron *pol,
                                               PlutoContext *context);
PlutoConstraints *polylib_matrix_to_pluto_constraints(Matrix *polymat,
                                                      PlutoContext *context);

Matrix *pluto_matrix_to_polylib(const PlutoMatrix *mat);
PlutoMatrix *polylib_matrix_to_pluto(Matrix *pmat, PlutoContext *context);

#endif
