#include "isl_blk.h"
#include "isl_ctx.h"

struct isl_blk isl_blk_empty()
{
	struct isl_blk block;
	block.size = 0;
	block.data = NULL;
	return block;
}

static int isl_blk_is_empty(struct isl_blk block)
{
	return block.size == 0 && block.data == NULL;
}

static struct isl_blk isl_blk_error()
{
	struct isl_blk block;
	block.size = -1;
	block.data = NULL;
	return block;
}

int isl_blk_is_error(struct isl_blk block)
{
	return block.size == -1 && block.data == NULL;
}

struct isl_blk isl_blk_alloc(struct isl_ctx *ctx, size_t n)
{
	int i;
	struct isl_blk block;

	if (ctx->n_cached) {
		int best = 0;
		for (i = 1; ctx->cache[best].size != n && i < ctx->n_cached; ++i) {
			if (ctx->cache[best].size < n) {
				if (ctx->cache[i].size > ctx->cache[best].size)
					best = i;
			} else if (ctx->cache[i].size >= n &&
				   ctx->cache[i].size < ctx->cache[best].size)
					best = i;
		}
		block = ctx->cache[best];
		if (--ctx->n_cached != best)
			ctx->cache[best] = ctx->cache[ctx->n_cached];
	} else
		block = isl_blk_empty();

	return isl_blk_extend(ctx, block, n);
}

struct isl_blk isl_blk_extend(struct isl_ctx *ctx, struct isl_blk block,
				size_t new_n)
{
	int i;
	isl_int *p;

	if (block.size >= new_n)
		return block;

	p = block.data;
	block.data = isl_realloc_array(ctx, block.data, isl_int, new_n);
	if (!block.data) {
		free(p);
		return isl_blk_error();
	}

	for (i = block.size; i < new_n; ++i)
		isl_int_init(block.data[i]);
	block.size = new_n;

	return block;
}

static void isl_blk_free_force(struct isl_ctx *ctx, struct isl_blk block)
{
	int i;

	for (i = 0; i < block.size; ++i)
		isl_int_clear(block.data[i]);
	free(block.data);
}

void isl_blk_free(struct isl_ctx *ctx, struct isl_blk block)
{
	if (isl_blk_is_empty(block) || isl_blk_is_error(block))
		return;

	if (ctx->n_cached < ISL_BLK_CACHE_SIZE)
		ctx->cache[ctx->n_cached++] = block;
	else
		isl_blk_free_force(ctx, block);
}

void isl_blk_clear_cache(struct isl_ctx *ctx)
{
	int i;

	for (i = 0; i < ctx->n_cached; ++i)
		isl_blk_free_force(ctx, ctx->cache[i]);
	ctx->n_cached = 0;
}
