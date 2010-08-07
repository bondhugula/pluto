#ifndef ISL_DIV_H
#define ISL_DIV_H

#include "isl_dim.h"
#include "isl_set.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_div {
	int ref;
	struct isl_ctx *ctx;

	struct isl_basic_map	*bmap;
	isl_int			**line;
};

struct isl_div *isl_div_alloc(struct isl_dim *dim);
struct isl_div *isl_basic_map_div(struct isl_basic_map *bmap, int pos);
struct isl_div *isl_basic_set_div(struct isl_basic_set *bset, int pos);
struct isl_div *isl_div_free(struct isl_div *c);

void isl_div_get_constant(struct isl_div *div, isl_int *v);
void isl_div_get_denominator(struct isl_div *div, isl_int *v);
void isl_div_get_coefficient(struct isl_div *div,
	enum isl_dim_type type, int pos, isl_int *v);
void isl_div_set_constant(struct isl_div *div, isl_int v);
void isl_div_set_denominator(struct isl_div *div, isl_int v);
void isl_div_set_coefficient(struct isl_div *div,
	enum isl_dim_type type, int pos, isl_int v);

#if defined(__cplusplus)
}
#endif

#endif
