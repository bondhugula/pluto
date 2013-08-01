#include "polylib/polylib64.h"
#include "constraints.h"

Polyhedron *pluto_constraints_to_polylib(const PlutoConstraints *cst);
PlutoConstraints *polylib_to_pluto_constraints(Polyhedron *pol);
PlutoConstraints *polylib_matrix_to_pluto_constraints(Matrix *polymat);

Matrix *pluto_matrix_to_polylib(const PlutoMatrix *mat);
PlutoMatrix *polylib_matrix_to_pluto(Matrix *pmat);

