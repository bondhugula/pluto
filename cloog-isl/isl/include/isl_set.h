#ifndef ISL_SET_H
#define ISL_SET_H

#include "isl_map.h"
#include "isl_list.h"

#if defined(__cplusplus)
extern "C" {
#endif

/* A "basic set" is a basic map with a zero-dimensional
 * domain.
 */
struct isl_basic_set {
	int ref;
#define ISL_BASIC_SET_FINAL		(1 << 0)
#define ISL_BASIC_SET_EMPTY		(1 << 1)
#define ISL_BASIC_SET_NO_IMPLICIT	(1 << 2)
#define ISL_BASIC_SET_NO_REDUNDANT	(1 << 3)
#define ISL_BASIC_SET_RATIONAL		(1 << 4)
#define ISL_BASIC_SET_NORMALIZED	(1 << 5)
#define ISL_BASIC_SET_NORMALIZED_DIVS	(1 << 6)
#define ISL_BASIC_SET_ALL_EQUALITIES	(1 << 7)
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

/* A "set" is a (possibly disjoint) union of basic sets.
 *
 * See the documentation of isl_map.
 */
struct isl_set {
	int ref;
#define ISL_SET_DISJOINT		(1 << 0)
#define ISL_SET_NORMALIZED		(1 << 1)
	unsigned flags;

	struct isl_ctx *ctx;

	struct isl_dim *dim;

	int n;

	size_t size;
	struct isl_basic_set *p[0];
};

unsigned isl_basic_set_n_dim(const struct isl_basic_set *bset);
unsigned isl_basic_set_n_param(const struct isl_basic_set *bset);
unsigned isl_basic_set_total_dim(const struct isl_basic_set *bset);
unsigned isl_basic_set_dim(const struct isl_basic_set *bset,
				enum isl_dim_type type);

unsigned isl_set_n_dim(const struct isl_set *set);
unsigned isl_set_n_param(const struct isl_set *set);
unsigned isl_set_dim(const struct isl_set *set, enum isl_dim_type type);

struct isl_dim *isl_basic_set_get_dim(struct isl_basic_set *bset);
struct isl_dim *isl_set_get_dim(struct isl_set *set);

struct isl_basic_set *isl_basic_set_alloc(struct isl_ctx *ctx,
		unsigned nparam, unsigned dim, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_set *isl_basic_set_extend(struct isl_basic_set *base,
		unsigned nparam, unsigned dim, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_set *isl_basic_set_extend_constraints(
		struct isl_basic_set *base, unsigned n_eq, unsigned n_ineq);
struct isl_basic_set *isl_basic_set_finalize(struct isl_basic_set *bset);
void isl_basic_set_free(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_copy(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_dup(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_empty(struct isl_dim *dim);
struct isl_basic_set *isl_basic_set_empty_like(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_universe(struct isl_dim *dim);
struct isl_basic_set *isl_basic_set_universe_like(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_interval(struct isl_ctx *ctx,
	isl_int min, isl_int max);
void isl_basic_set_dump(struct isl_basic_set *bset,
				FILE *out, int indent);
struct isl_basic_set *isl_basic_set_swap_vars(
		struct isl_basic_set *bset, unsigned n);
struct isl_basic_set *isl_basic_set_intersect(
		struct isl_basic_set *bset1,
		struct isl_basic_set *bset2);
struct isl_basic_set *isl_basic_set_apply(
		struct isl_basic_set *bset,
		struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_set_affine_hull(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_simplify(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_product(struct isl_basic_set_list *list);

#define ISL_FORMAT_POLYLIB		1
#define ISL_FORMAT_OMEGA		2
struct isl_basic_set *isl_basic_set_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format);
struct isl_basic_set *isl_basic_set_read_from_str(struct isl_ctx *ctx,
		const char *str, unsigned nparam, unsigned input_format);
struct isl_set *isl_set_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format);
#define ISL_FORMAT_POLYLIB_CONSTRAINTS	3
void isl_basic_set_print(struct isl_basic_set *bset, FILE *out, int indent,
	const char *prefix, const char *suffix, unsigned output_format);
void isl_set_print(struct isl_set *set, FILE *out, int indent,
	unsigned output_format);
struct isl_basic_set *isl_basic_set_fix_si(struct isl_basic_set *bset,
		enum isl_dim_type type, unsigned pos, int value);

struct isl_basic_set *isl_basic_set_from_underlying_set(
	struct isl_basic_set *bset, struct isl_basic_set *like);
struct isl_set *isl_set_from_underlying_set(
	struct isl_set *set, struct isl_basic_set *like);
struct isl_set *isl_set_to_underlying_set(struct isl_set *set);

int isl_basic_set_is_equal(
		struct isl_basic_set *bset1, struct isl_basic_set *bset2);

struct isl_set *isl_basic_set_lexmin(struct isl_basic_set *bset);
struct isl_set *isl_basic_set_union(
		struct isl_basic_set *bset1,
		struct isl_basic_set *bset2);

int isl_basic_set_compare_at(struct isl_basic_set *bset1,
	struct isl_basic_set *bset2, int pos);

int isl_basic_set_is_universe(struct isl_basic_set *bset);
int isl_basic_set_fast_is_empty(struct isl_basic_set *bset);
int isl_basic_set_is_empty(struct isl_basic_set *bset);

struct isl_set *isl_set_alloc(struct isl_ctx *ctx,
		unsigned nparam, unsigned dim, int n, unsigned flags);
struct isl_set *isl_set_extend(struct isl_set *base,
		unsigned nparam, unsigned dim);
struct isl_set *isl_set_empty(struct isl_dim *dim);
struct isl_set *isl_set_empty_like(struct isl_set *set);
struct isl_set *isl_set_universe(struct isl_dim *dim);
struct isl_set *isl_set_add(struct isl_set *set, struct isl_basic_set *bset);
struct isl_set *isl_set_finalize(struct isl_set *set);
struct isl_set *isl_set_copy(struct isl_set *set);
void isl_set_free(struct isl_set *set);
struct isl_set *isl_set_dup(struct isl_set *set);
struct isl_set *isl_set_from_basic_set(struct isl_basic_set *bset);
struct isl_basic_set *isl_set_affine_hull(struct isl_set *set);
struct isl_basic_set *isl_set_convex_hull(struct isl_set *set);
struct isl_basic_set *isl_set_simple_hull(struct isl_set *set);
struct isl_basic_set *isl_set_bounded_simple_hull(struct isl_set *set);

struct isl_set *isl_set_union_disjoint(
			struct isl_set *set1, struct isl_set *set2);
struct isl_set *isl_set_union(struct isl_set *set1, struct isl_set *set2);
struct isl_set *isl_set_product(struct isl_set *set1, struct isl_set *set2);
struct isl_set *isl_set_intersect(struct isl_set *set1, struct isl_set *set2);
struct isl_set *isl_set_subtract(struct isl_set *set1, struct isl_set *set2);
struct isl_set *isl_set_apply(struct isl_set *set, struct isl_map *map);
struct isl_set *isl_set_fix_dim_si(struct isl_set *set,
		unsigned dim, int value);
struct isl_set *isl_set_lower_bound_dim(struct isl_set *set,
		unsigned dim, isl_int value);
struct isl_basic_set *isl_basic_set_remove_dims(struct isl_basic_set *bset,
		unsigned first, unsigned n);
struct isl_basic_set *isl_basic_set_remove_divs(struct isl_basic_set *bset);
struct isl_set *isl_set_eliminate_dims(struct isl_set *set,
		unsigned first, unsigned n);
struct isl_set *isl_set_remove_dims(struct isl_set *set,
		unsigned first, unsigned n);
struct isl_set *isl_set_remove_divs(struct isl_set *set);

void isl_set_dump(struct isl_set *set, FILE *out, int indent);
struct isl_set *isl_set_swap_vars(struct isl_set *set, unsigned n);
int isl_set_is_empty(struct isl_set *set);
int isl_set_is_subset(struct isl_set *set1, struct isl_set *set2);
int isl_set_is_equal(struct isl_set *set1, struct isl_set *set2);

struct isl_set *isl_basic_set_compute_divs(struct isl_basic_set *bset);
struct isl_set *isl_set_compute_divs(struct isl_set *set);

struct isl_basic_set *isl_set_copy_basic_set(struct isl_set *set);
struct isl_set *isl_set_drop_basic_set(struct isl_set *set,
						struct isl_basic_set *bset);

int isl_basic_set_fast_dim_is_fixed(struct isl_basic_set *bset, unsigned dim,
	isl_int *val);

int isl_set_fast_dim_is_fixed(struct isl_set *set, unsigned dim, isl_int *val);
int isl_set_fast_dim_has_fixed_lower_bound(struct isl_set *set,
	unsigned dim, isl_int *val);

struct isl_basic_set *isl_basic_set_gist(struct isl_basic_set *bset,
						struct isl_basic_set *context);
struct isl_set *isl_set_gist(struct isl_set *set,
	struct isl_basic_set *context);
int isl_basic_set_dim_residue_class(struct isl_basic_set *bset,
	int pos, isl_int *modulo, isl_int *residue);

struct isl_set *isl_set_coalesce(struct isl_set *set);

int isl_set_fast_is_equal(struct isl_set *set1, struct isl_set *set2);
int isl_set_fast_is_disjoint(struct isl_set *set1, struct isl_set *set2);

uint32_t isl_set_get_hash(struct isl_set *set);

int isl_set_dim_is_unique(struct isl_set *set, unsigned dim);

#if defined(__cplusplus)
}
#endif

#endif
