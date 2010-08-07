#ifndef ISL_LIST_H
#define ISL_LIST_H

#include <isl_ctx.h>

struct isl_basic_set;

struct isl_basic_set_list {
	int ref;
	struct isl_ctx *ctx;

	int n;

	size_t size;
	struct isl_basic_set *p[0];
};

struct isl_basic_set_list *isl_basic_set_list_alloc(struct isl_ctx *ctx, int n);
void isl_basic_set_list_free(struct isl_basic_set_list *list);
struct isl_basic_set_list *isl_basic_set_list_add(
	struct isl_basic_set_list *list,
	struct isl_basic_set *bset);

#endif
