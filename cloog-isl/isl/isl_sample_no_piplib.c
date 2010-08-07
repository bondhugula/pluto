#include "isl_sample_piplib.h"

struct isl_vec *isl_pip_basic_set_sample(struct isl_basic_set *bset)
{
	isl_basic_set_free(bset);
	return NULL;
}
