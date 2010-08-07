#ifndef ISL_MAP_POLYLIB_H
#define ISL_MAP_POLYLIB_H

#include <isl_dim.h>
#include <isl_map.h>
#include <isl_polylib.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_basic_map *isl_basic_map_new_from_polylib(Polyhedron *P,
			struct isl_dim *dim);
struct isl_map *isl_map_new_from_polylib(Polyhedron *D, struct isl_dim *dim);
Polyhedron *isl_basic_map_to_polylib(struct isl_basic_map *bmap);
Polyhedron *isl_map_to_polylib(struct isl_map *map);

#if defined(__cplusplus)
}
#endif

#endif
