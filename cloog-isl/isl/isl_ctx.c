#include "isl_ctx.h"
#include "isl_vec.h"
#ifdef ISL_POLYLIB
#include <polylib/polylibgmp.h>
#endif

struct isl_ctx *isl_ctx_alloc()
{
	struct isl_ctx *ctx = NULL;

	ctx = isl_calloc_type(NULL, struct isl_ctx);
	if (!ctx)
		goto error;

	if (isl_hash_table_init(ctx, &ctx->name_hash, 0))
		goto error;

	ctx->stats = isl_calloc_type(ctx, struct isl_stats);
	if (!ctx->stats)
		goto error;

	ctx->ref = 0;

	isl_int_init(ctx->one);
	isl_int_set_si(ctx->one, 1);

	isl_int_init(ctx->negone);
	isl_int_set_si(ctx->negone, -1);

	ctx->n_cached = 0;

#ifdef ISL_POLYLIB
	ctx->MaxRays = POL_NO_DUAL | POL_INTEGER;
#endif

	ctx->lp_solver = ISL_LP_TAB;
	ctx->ilp_solver = ISL_ILP_GBR;

	return ctx;
error:
	free(ctx);
	return NULL;
}

void isl_ctx_ref(struct isl_ctx *ctx)
{
	ctx->ref++;
}

void isl_ctx_deref(struct isl_ctx *ctx)
{
	isl_assert(ctx, ctx->ref > 0, return);
	ctx->ref--;
}

void isl_ctx_free(struct isl_ctx *ctx)
{
	if (!ctx)
		return;
	isl_assert(ctx, ctx->ref == 0, return);
	isl_hash_table_clear(&ctx->name_hash);
	isl_blk_clear_cache(ctx);
	isl_int_clear(ctx->one);
	isl_int_clear(ctx->negone);
	free(ctx->stats);
	free(ctx);
}
