#ifndef ISL_BASIS_REDUCTION_H
#define ISL_BASIS_REDUCTION_H

#include "isl_set.h"
#include "isl_mat.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_mat *isl_basic_set_reduced_basis(struct isl_basic_set *bset);

#if defined(__cplusplus)
}
#endif

#endif
