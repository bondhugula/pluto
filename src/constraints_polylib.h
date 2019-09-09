#ifndef _CONSTRAINTS_POLYLIB_H
#define _CONSTRAINTS_POLYLIB_H

#include "polylib/polylib64.h"

typedef struct pluto_constraints PlutoConstraints;
typedef struct pluto_matrix PlutoMatrix;

Polyhedron *pluto_constraints_to_polylib(const PlutoConstraints *cst);
PlutoConstraints *polylib_to_pluto_constraints(Polyhedron *pol);
PlutoConstraints *polylib_matrix_to_pluto_constraints(Matrix *polymat);

Matrix *pluto_matrix_to_polylib(const PlutoMatrix *mat);
PlutoMatrix *polylib_matrix_to_pluto(Matrix *pmat);

#endif
