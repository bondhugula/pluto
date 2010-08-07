#ifndef ISL_NAME_H
#define ISL_NAME_H

#include <isl_ctx.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_name {
	int ref;

	const char *name;
};

struct isl_name *isl_name_alloc(struct isl_ctx *ctx, const char *name);
struct isl_name *isl_name_get(struct isl_ctx *ctx, const char *name);
struct isl_name *isl_name_copy(struct isl_ctx *ctx, struct isl_name *name);
void isl_name_free(struct isl_ctx *ctx, struct isl_name *name);

#if defined(__cplusplus)
}
#endif

#endif
