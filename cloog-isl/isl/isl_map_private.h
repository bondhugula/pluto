#include "isl_set.h"
#include "isl_map.h"
#include "isl_vec.h"

unsigned isl_basic_map_offset(struct isl_basic_map *bmap,
					enum isl_dim_type type);

int isl_map_compatible_domain(struct isl_map *map, struct isl_set *set);
int isl_basic_map_compatible_domain(struct isl_basic_map *bmap,
		struct isl_basic_set *bset);
int isl_basic_map_compatible_range(struct isl_basic_map *bmap,
		struct isl_basic_set *bset);

struct isl_basic_map *isl_basic_map_extend_dim(struct isl_basic_map *base,
		struct isl_dim *dim, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_set *isl_basic_set_extend_dim(struct isl_basic_set *base,
		struct isl_dim *dim, unsigned extra,
		unsigned n_eq, unsigned n_ineq);
struct isl_basic_set *isl_basic_set_add_constraints(struct isl_basic_set *bset1,
		struct isl_basic_set *bset2, unsigned pos);

int isl_basic_set_contains(struct isl_basic_set *bset, struct isl_vec *vec);

struct isl_basic_set *isl_basic_set_alloc_dim(struct isl_dim *dim,
		unsigned extra, unsigned n_eq, unsigned n_ineq);
struct isl_set *isl_set_alloc_dim(struct isl_dim *dim, int n, unsigned flags);
struct isl_basic_map *isl_basic_map_alloc_dim(struct isl_dim *dim,
		unsigned extra, unsigned n_eq, unsigned n_ineq);
struct isl_map *isl_map_alloc_dim(struct isl_dim *dim, int n, unsigned flags);

unsigned isl_basic_map_total_dim(const struct isl_basic_map *bmap);

int isl_basic_map_alloc_equality(struct isl_basic_map *bmap);
int isl_basic_set_alloc_equality(struct isl_basic_set *bset);
int isl_basic_set_free_inequality(struct isl_basic_set *bset, unsigned n);
int isl_basic_map_free_equality(struct isl_basic_map *bmap, unsigned n);
int isl_basic_set_free_equality(struct isl_basic_set *bset, unsigned n);
int isl_basic_set_alloc_inequality(struct isl_basic_set *bset);
int isl_basic_map_alloc_inequality(struct isl_basic_map *bmap);
int isl_basic_map_free_inequality(struct isl_basic_map *bmap, unsigned n);
int isl_basic_map_alloc_div(struct isl_basic_map *bmap);
int isl_basic_set_alloc_div(struct isl_basic_set *bset);
int isl_basic_map_free_div(struct isl_basic_map *bmap, unsigned n);
void isl_basic_map_inequality_to_equality(
		struct isl_basic_map *bmap, unsigned pos);
int isl_basic_map_drop_equality(struct isl_basic_map *bmap, unsigned pos);
int isl_basic_set_drop_equality(struct isl_basic_set *bset, unsigned pos);
int isl_basic_set_drop_inequality(struct isl_basic_set *bset, unsigned pos);
int isl_basic_map_drop_inequality(struct isl_basic_map *bmap, unsigned pos);

int isl_inequality_negate(struct isl_basic_map *bmap, unsigned pos);

struct isl_basic_set *isl_basic_set_cow(struct isl_basic_set *bset);
struct isl_basic_map *isl_basic_map_cow(struct isl_basic_map *bmap);
struct isl_set *isl_set_cow(struct isl_set *set);
struct isl_map *isl_map_cow(struct isl_map *map);

struct isl_basic_map *isl_basic_map_set_to_empty(struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_set_set_to_empty(struct isl_basic_set *bset);
struct isl_map *isl_basic_map_compute_divs(struct isl_basic_map *bmap);
struct isl_map *isl_map_compute_divs(struct isl_map *map);
struct isl_basic_map *isl_basic_map_align_divs(
		struct isl_basic_map *dst, struct isl_basic_map *src);
struct isl_basic_set *isl_basic_set_align_divs(
		struct isl_basic_set *dst, struct isl_basic_set *src);
struct isl_map *isl_map_align_divs(struct isl_map *map);
struct isl_set *isl_set_align_divs(struct isl_set *set);
struct isl_basic_map *isl_basic_map_gauss(
	struct isl_basic_map *bmap, int *progress);
struct isl_basic_set *isl_basic_set_gauss(
	struct isl_basic_set *bset, int *progress);
struct isl_basic_map *isl_basic_map_normalize_constraints(
	struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_set_normalize_constraints(
	struct isl_basic_set *bset);
struct isl_basic_map *isl_basic_map_implicit_equalities(
						struct isl_basic_map *bmap);
struct isl_basic_set *isl_basic_map_underlying_set(struct isl_basic_map *bmap);
struct isl_set *isl_map_underlying_set(struct isl_map *map);
struct isl_basic_map *isl_basic_map_overlying_set(struct isl_basic_set *bset,
	struct isl_basic_map *like);
struct isl_basic_map *isl_basic_map_drop(struct isl_basic_map *bmap,
	enum isl_dim_type type, unsigned first, unsigned n);
struct isl_basic_set *isl_basic_set_drop_dims(
		struct isl_basic_set *bset, unsigned first, unsigned n);
struct isl_set *isl_set_drop_dims(
		struct isl_set *set, unsigned first, unsigned n);
struct isl_map *isl_map_drop_inputs(
		struct isl_map *map, unsigned first, unsigned n);
struct isl_map *isl_map_drop(struct isl_map *map,
	enum isl_dim_type type, unsigned first, unsigned n);

struct isl_map *isl_map_remove_empty_parts(struct isl_map *map);
struct isl_set *isl_set_remove_empty_parts(struct isl_set *set);

struct isl_set *isl_set_normalize(struct isl_set *set);

struct isl_set *isl_set_drop_vars(
		struct isl_set *set, unsigned first, unsigned n);

struct isl_basic_map *isl_basic_map_eliminate_vars(
	struct isl_basic_map *bmap, unsigned pos, unsigned n);
struct isl_basic_set *isl_basic_set_eliminate_vars(
	struct isl_basic_set *bset, unsigned pos, unsigned n);

int isl_basic_set_constraint_is_redundant(struct isl_basic_set **bset,
	isl_int *c, isl_int *opt_n, isl_int *opt_d);

struct isl_basic_map *isl_basic_map_drop_redundant_divs(
	struct isl_basic_map *bmap);

struct isl_basic_set *isl_basic_set_recession_cone(struct isl_basic_set *bset);
struct isl_basic_set *isl_basic_set_lineality_space(struct isl_basic_set *bset);
