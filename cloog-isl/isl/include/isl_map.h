#ifndef ISL_MAP_H
#define ISL_MAP_H

#include <stdio.h>

#include <isl_int.h>
#include <isl_ctx.h>
#include <isl_blk.h>
#include <isl_dim.h>

#if defined(__cplusplus)
extern "C" {
#endif

/* General notes:
 *
 * All structures are reference counted to allow reuse without duplication.
 * A *_copy operation will increase the reference count, while a *_free
 * operation will decrease the reference count and only actually release
 * the structures when the reference count drops to zero.
 *
 * Functions that return an isa structure will in general _destroy_
 * all argument isa structures (the obvious execption begin the _copy
 * functions).  A pointer passed to such a function may therefore
 * never be used after the function call.  If you want to keep a
 * reference to the old structure(s), use the appropriate _copy function.
 */

/* A "basic map" is a relation between two sets of variables,
 * called the "in" and "out" variables.
 *
 * It is implemented as a set with two extra fields:
 * n_in is the number of in variables
 * n_out is the number of out variables
 * n_in + n_out should be equal to set.dim
 */
struct isl_vec;
struct isl_dim;
struct isl_basic_map {
	int ref;
#define ISL_BASIC_MAP_FINAL		(1 << 0)
#define ISL_BASIC_MAP_EMPTY		(1 << 1)
#define ISL_BASIC_MAP_NO_IMPLICIT	(1 << 2)
#define ISL_BASIC_MAP_NO_REDUNDANT	(1 << 3)
#define ISL_BASIC_MAP_RATIONAL		(1 << 4)
#define ISL_BASIC_MAP_NORMALIZED	(1 << 5)
#define ISL_BASIC_MAP_NORMALIZED_DIVS	(1 << 6)
#define ISL_BASIC_MAP_ALL_EQUALITIES	(1 << 7)
	unsigned flags;

	struct isl_ctx *ctx;

	struct isl_dim *dim;
	unsigned extra;

	unsigned n_eq;
	unsigned n_ineq;

	size_t c_size;
	isl_int **eq;
	isl_int **ineq;

	unsigned n_div;

	isl_int **div;

	struct isl_vec *sample;

	struct isl_blk block;
	struct isl_blk block2;
};
struct isl_basic_set;

/* A "map" is a (disjoint) union of basic maps.
 *
 * Currently, the isl_set structure is identical to the isl_map structure
 * and the library depends on this correspondence internally.
 * However, users should not depend on this correspondence.
 */
struct isl_map {
	int ref;
#define ISL_MAP_DISJOINT		(1 << 0)
#define ISL_MAP_NORMALIZED		(1 << 1)
	unsigned flags;

	struct isl_ctx *ctx;

	struct isl_dim *dim;

	int n;

	size_t size;
	struct isl_basic_map *p[0];
};
struct isl_set;

unsigned isl_basic_map_n_in(const struct isl_basic_map *bmap);
unsigned isl_basic_map_n_out(const struct isl_basic_map *bmap);
unsigned isl_basic_map_n_param(const struct isl_basic_map *bmap);
unsigned isl_basic_map_n_div(const struct isl_basic_map *bmap);
unsigned isl_basic_map_total_dim(const struct isl_basic_map *bmap);
unsigned isl_basic_map_dim(const struct isl_basic_map *bmap,
				enum isl_dim_type type);

unsigned isl_map_n_in(const struct isl_map *map);
unsigned isl_map_n_out(const struct isl_map *map);
unsigned isl_map_n_param(const struct isl_map *map);
unsigned isl_map_dim(const struct isl_map *map, enum isl_dim_type type);

struct isl_basic_map *isl_basic_map_alloc(struct isl_ctx *ctx,
		unsigned nparam, unsigned in, unsigned out, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_map *isl_basic_map_identity(struct isl_dim *set_dim);
struct isl_basic_map *isl_basic_map_identity_like(struct isl_basic_map *model);
struct isl_basic_map *isl_basic_map_finalize(struct isl_basic_map *bmap);
void isl_basic_map_free(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_copy(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_extend(struct isl_basic_map *base,
		unsigned nparam, unsigned n_in, unsigned n_out, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_map *isl_basic_map_extend_constraints(
		struct isl_basic_map *base, unsigned n_eq, unsigned n_ineq);
struct isl_basic_map *isl_basic_map_equal(
		struct isl_dim *dim, unsigned n_equal);
struct isl_basic_map *isl_basic_map_less_at(struct isl_dim *dim, unsigned pos);
struct isl_basic_map *isl_basic_map_more_at(struct isl_dim *dim, unsigned pos);
struct isl_basic_map *isl_basic_map_empty(struct isl_ctx *ctx,
		unsigned nparam, unsigned in, unsigned out);
struct isl_basic_map *isl_basic_map_empty_like(struct isl_basic_map *model);
struct isl_basic_map *isl_basic_map_empty_like_map(struct isl_map *model);
struct isl_basic_map *isl_basic_map_universe(struct isl_dim *dim);
struct isl_basic_map *isl_basic_map_convex_hull(struct isl_basic_map *bmap);

struct isl_basic_map *isl_basic_map_intersect_domain(
		struct isl_basic_map *bmap,
		struct isl_basic_set *bset);
struct isl_basic_map *isl_basic_map_intersect_range(
		struct isl_basic_map *bmap,
		struct isl_basic_set *bset);
struct isl_basic_map *isl_basic_map_intersect(
		struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);
struct isl_map *isl_basic_map_union(
		struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);
struct isl_basic_map *isl_basic_map_apply_domain(
		struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);
struct isl_basic_map *isl_basic_map_apply_range(
		struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);
struct isl_basic_map *isl_basic_map_reverse(struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_map_domain(struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_map_range(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_remove(struct isl_basic_map *bmap,
	enum isl_dim_type type, unsigned first, unsigned n);
struct isl_basic_map *isl_basic_map_from_basic_set(struct isl_basic_set *bset,
		struct isl_dim *dim);
struct isl_basic_set *isl_basic_set_from_basic_map(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_simplify(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_detect_equalities(
						struct isl_basic_map *bmap);
#define ISL_FORMAT_POLYLIB	1
#define ISL_FORMAT_OMEGA	2
struct isl_basic_map *isl_basic_map_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format);
struct isl_basic_map *isl_basic_map_fix_si(struct isl_basic_map *bmap,
		enum isl_dim_type type, unsigned pos, int value);

struct isl_basic_map *isl_basic_map_sum(
		struct isl_basic_map *bmap1, struct isl_basic_map *bmap2);
struct isl_basic_map *isl_basic_map_neg(struct isl_basic_map *bmap);
struct isl_basic_map *isl_basic_map_floordiv(struct isl_basic_map *bmap,
		isl_int d);

int isl_basic_map_is_equal(
		struct isl_basic_map *bmap1, struct isl_basic_map *bmap2);

struct isl_map *isl_basic_map_lexmax(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty);
struct isl_map *isl_basic_map_lexmin(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty);

void isl_basic_map_dump(struct isl_basic_map *bmap, FILE *out, int indent);

struct isl_basic_map *isl_map_copy_basic_map(struct isl_map *map);
struct isl_map *isl_map_drop_basic_map(struct isl_map *map,
						struct isl_basic_map *bmap);

int isl_basic_map_fast_is_fixed(struct isl_basic_map *bmap,
	enum isl_dim_type type, unsigned pos, isl_int *val);

int isl_basic_map_is_universe(struct isl_basic_map *bmap);
int isl_basic_map_fast_is_empty(struct isl_basic_map *bmap);
int isl_basic_map_is_empty(struct isl_basic_map *bmap);
int isl_basic_map_is_subset(struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);
int isl_basic_map_is_strict_subset(struct isl_basic_map *bmap1,
		struct isl_basic_map *bmap2);

struct isl_map *isl_map_alloc(struct isl_ctx *ctx,
		unsigned nparam, unsigned in, unsigned out, int n,
		unsigned flags);
struct isl_map *isl_map_universe(struct isl_dim *dim);
struct isl_map *isl_map_empty(struct isl_dim *dim);
struct isl_map *isl_map_empty_like(struct isl_map *model);
struct isl_map *isl_map_empty_like_basic_map(struct isl_basic_map *model);
struct isl_map *isl_map_dup(struct isl_map *map);
struct isl_map *isl_map_add(struct isl_map *map, struct isl_basic_map *bmap);
struct isl_map *isl_map_identity(struct isl_dim *set_dim);
struct isl_map *isl_map_identity_like(struct isl_basic_map *model);
struct isl_map *isl_map_finalize(struct isl_map *map);
void isl_map_free(struct isl_map *map);
struct isl_map *isl_map_copy(struct isl_map *map);
struct isl_map *isl_map_extend(struct isl_map *base,
		unsigned nparam, unsigned n_in, unsigned n_out);
struct isl_map *isl_map_reverse(struct isl_map *map);
struct isl_map *isl_map_union(struct isl_map *map1, struct isl_map *map2);
struct isl_map *isl_map_union_disjoint(
			struct isl_map *map1, struct isl_map *map2);
struct isl_map *isl_map_intersect_domain(
		struct isl_map *map,
		struct isl_set *set);
struct isl_map *isl_map_intersect_range(
		struct isl_map *map,
		struct isl_set *set);
struct isl_map *isl_map_apply_domain(
		struct isl_map *map1,
		struct isl_map *map2);
struct isl_map *isl_map_apply_range(
		struct isl_map *map1,
		struct isl_map *map2);
struct isl_map *isl_map_product(struct isl_map *map1, struct isl_map *map2);
struct isl_map *isl_map_intersect(struct isl_map *map1, struct isl_map *map2);
struct isl_map *isl_map_subtract(struct isl_map *map1, struct isl_map *map2);
struct isl_map *isl_map_fix_input_si(struct isl_map *map,
		unsigned input, int value);
struct isl_map *isl_map_fix_si(struct isl_map *map,
		enum isl_dim_type type, unsigned pos, int value);
struct isl_basic_set *isl_basic_map_deltas(struct isl_basic_map *bmap);
struct isl_set *isl_map_deltas(struct isl_map *map);
struct isl_set *isl_map_range(struct isl_map *map);
struct isl_map *isl_map_detect_equalities(struct isl_map *map);
struct isl_basic_map *isl_map_affine_hull(struct isl_map *map);
struct isl_map *isl_map_remove(struct isl_map *map,
	enum isl_dim_type type, unsigned first, unsigned n);
struct isl_map *isl_map_remove_inputs(struct isl_map *map,
	unsigned first, unsigned n);

struct isl_set *isl_map_domain(struct isl_map *bmap);
struct isl_map *isl_map_from_basic_map(struct isl_basic_map *bmap);
struct isl_map *isl_map_from_range(struct isl_set *set);
struct isl_map *isl_map_from_set(struct isl_set *set, struct isl_dim *dim);
struct isl_set *isl_set_from_map(struct isl_map *map);

int isl_map_fast_is_empty(struct isl_map *map);
int isl_map_is_empty(struct isl_map *map);
int isl_map_is_subset(struct isl_map *map1, struct isl_map *map2);
int isl_map_is_equal(struct isl_map *map1, struct isl_map *map2);

void isl_map_dump(struct isl_map *map, FILE *out, int indent);

int isl_map_fast_input_is_fixed(struct isl_map *map,
		unsigned in, isl_int *val);

struct isl_map *isl_map_coalesce(struct isl_map *map);

int isl_map_fast_is_equal(struct isl_map *map1, struct isl_map *map2);

#if defined(__cplusplus)
}
#endif

#endif
