#ifndef ISL_SUPPORT_H
#define ISL_SUPPORT_H

#include "isl/ctx.h"
#include "isl/set.h"

typedef struct pluto_matrix PlutoMatrix;
typedef struct plutoContext PlutoContext;

#if defined(__cplusplus)
extern "C" {
#endif

struct pluto_mat_context_info {
  PlutoMatrix **mat;
  PlutoContext *context;
};

isl_stat isl_aff_to_pluto_func(__isl_take isl_set *set, __isl_take isl_aff *aff,
                               void *user);
PlutoMatrix *pluto_matrix_from_isl_mat(__isl_keep isl_mat *mat,
                                       PlutoContext *context);
long long isl_val_get_num_ll(__isl_keep isl_val *v);

#if defined(__cplusplus)
}
#endif

#endif // ISL_SUPPORT_H
