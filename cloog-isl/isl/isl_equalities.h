#ifndef ISL_EQUALITIES_H
#define ISL_EQUALITIES_H

#include <isl_set.h>
#include <isl_mat.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_mat *isl_mat_variable_compression(
			struct isl_mat *B, struct isl_mat **T2);
struct isl_mat *isl_mat_parameter_compression(
			struct isl_mat *B, struct isl_vec *d);
struct isl_basic_set *isl_basic_set_remove_equalities(
	struct isl_basic_set *bset, struct isl_mat **T, struct isl_mat **T2);

#if defined(__cplusplus)
}
#endif

#endif
