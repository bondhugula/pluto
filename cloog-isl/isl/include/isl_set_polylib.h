#ifndef ISL_SET_POLYLIB_H
#define ISL_SET_POLYLIB_H

#include <isl_set.h>
#include <isl_polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_basic_set *isl_basic_set_new_from_polylib(Polyhedron *P,
			struct isl_dim *dim);
Polyhedron *isl_basic_set_to_polylib(struct isl_basic_set *bset);

struct isl_set *isl_set_new_from_polylib(Polyhedron *D, struct isl_dim *dim);
Polyhedron *isl_set_to_polylib(struct isl_set *set);

#if defined(__cplusplus)
}
#endif

#endif
