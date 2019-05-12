#ifndef _ISL_SUPPORT_H
#define _ISL_SUPPORT_H

#include "math_support.h"

#if defined(__cplusplus)
extern "C" {
#endif

isl_stat isl_aff_to_pluto_func(__isl_take isl_set *set, __isl_take isl_aff *aff,
                               void *user);
PlutoMatrix *pluto_matrix_from_isl_mat(__isl_keep isl_mat *mat);
long long isl_val_get_num_ll(__isl_keep isl_val *v);

#if defined(__cplusplus)
}
#endif

#endif // _ISL_SUPPORT_H
