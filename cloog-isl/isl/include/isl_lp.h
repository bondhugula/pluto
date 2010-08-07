#ifndef ISL_LP_H
#define ISL_LP_H

#include <isl_map.h>

enum isl_lp_result {
	isl_lp_error = -1,
	isl_lp_ok = 0,
	isl_lp_unbounded,
	isl_lp_empty
};

#if defined(__cplusplus)
extern "C" {
#endif

enum isl_lp_result isl_solve_lp(struct isl_basic_map *bmap, int maximize,
				      isl_int *f, isl_int denom, isl_int *opt,
				      isl_int *opt_denom);

#if defined(__cplusplus)
}
#endif

#endif
