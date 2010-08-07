#ifndef ISL_SAMPLE_PIP_H
#define ISL_SAMPLE_PIP

#include <isl_set.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_vec *isl_pip_basic_set_sample(struct isl_basic_set *bset);

#if defined(__cplusplus)
}
#endif

#endif
