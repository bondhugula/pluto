#ifndef ISL_MAP_PIPLIB_H
#define ISL_MAP_PIPLIB_H

#include <isl_map.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_map *isl_pip_basic_map_lexmax(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty);
struct isl_map *isl_pip_basic_map_lexmin(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty);
struct isl_map *isl_pip_basic_map_compute_divs(struct isl_basic_map *bmap);

#if defined(__cplusplus)
}
#endif

#endif
