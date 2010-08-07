#include "isl_list.h"

struct isl_basic_set_list *isl_basic_set_list_alloc(struct isl_ctx *ctx, int n)
{
	struct isl_basic_set_list *list;

	isl_assert(ctx, n >= 0, return NULL);
	list = isl_alloc(ctx, struct isl_basic_set_list,
			 sizeof(struct isl_basic_set_list) +
			 n * sizeof(struct isl_basic_set *));
	if (!list)
		return NULL;

	list->ctx = ctx;
	isl_ctx_ref(ctx);
	list->ref = 1;
	list->size = n;
	list->n = 0;
	return list;
}

struct isl_basic_set_list *isl_basic_set_list_add(
	struct isl_basic_set_list *list,
	struct isl_basic_set *bset)
{
	if (!list || !bset)
		goto error;
	isl_assert(list->ctx, list->n < list->size, goto error);
	list->p[list->n] = bset;
	list->n++;
	return list;
error:
	isl_basic_set_free(bset);
	isl_basic_set_list_free(list);
	return NULL;
}

void isl_basic_set_list_free(struct isl_basic_set_list *list)
{
	int i;

	if (!list)
		return;

	if (--list->ref > 0)
		return;

	isl_ctx_deref(list->ctx);
	for (i = 0; i < list->n; ++i)
		isl_basic_set_free(list->p[i]);
	free(list);
}
