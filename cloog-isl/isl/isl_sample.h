#ifndef ISL_SAMPLE_H
#define ISL_SAMPLE

#include <isl_set.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_vec *isl_basic_set_sample(struct isl_basic_set *bset);

#if defined(__cplusplus)
}
#endif

#endif
