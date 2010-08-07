#include <assert.h>
#include "isl_sample.h"
#include "isl_vec.h"

int main(int argc, char **argv)
{
	struct isl_ctx *ctx = isl_ctx_alloc();
	struct isl_basic_set *bset;
	struct isl_vec *sample;

	bset = isl_basic_set_read_from_file(ctx, stdin, 0, ISL_FORMAT_POLYLIB);
	sample = isl_basic_set_sample(isl_basic_set_copy(bset));
	isl_vec_dump(sample, stdout, 0);
	assert(sample);
	if (sample->size > 0)
		assert(isl_basic_set_contains(bset, sample));
	isl_basic_set_free(bset);
	isl_vec_free(sample);
	isl_ctx_free(ctx);

	return 0;
}
