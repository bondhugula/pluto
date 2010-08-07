#include "isl_ctx.h"
#include "isl_seq.h"
#include "isl_set.h"
#include "isl_lp.h"
#include "isl_map.h"
#include "isl_map_private.h"
#include "isl_equalities.h"
#include "isl_sample.h"
#include "isl_tab.h"

struct isl_basic_map *isl_basic_map_implicit_equalities(
						struct isl_basic_map *bmap)
{
	struct isl_tab *tab;

	if (!bmap)
		return bmap;

	bmap = isl_basic_map_gauss(bmap, NULL);
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_EMPTY))
		return bmap;
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_NO_IMPLICIT))
		return bmap;
	if (bmap->n_ineq <= 1)
		return bmap;

	tab = isl_tab_from_basic_map(bmap);
	tab = isl_tab_detect_equalities(tab);
	bmap = isl_basic_map_update_from_tab(bmap, tab);
	isl_tab_free(tab);
	bmap = isl_basic_map_gauss(bmap, NULL);
	ISL_F_SET(bmap, ISL_BASIC_MAP_NO_IMPLICIT);
	return bmap;
}

struct isl_basic_set *isl_basic_set_implicit_equalities(
						struct isl_basic_set *bset)
{
	return (struct isl_basic_set *)
		isl_basic_map_implicit_equalities((struct isl_basic_map*)bset);
}

struct isl_map *isl_map_implicit_equalities(struct isl_map *map)
{
	int i;

	if (!map)
		return map;

	for (i = 0; i < map->n; ++i) {
		map->p[i] = isl_basic_map_implicit_equalities(map->p[i]);
		if (!map->p[i])
			goto error;
	}

	return map;
error:
	isl_map_free(map);
	return NULL;
}

/* Make eq[row][col] of both bmaps equal so we can add the row
 * add the column to the common matrix.
 * Note that because of the echelon form, the columns of row row
 * after column col are zero.
 */
static void set_common_multiple(
	struct isl_basic_set *bset1, struct isl_basic_set *bset2,
	unsigned row, unsigned col)
{
	isl_int m, c;

	if (isl_int_eq(bset1->eq[row][col], bset2->eq[row][col]))
		return;

	isl_int_init(c);
	isl_int_init(m);
	isl_int_lcm(m, bset1->eq[row][col], bset2->eq[row][col]);
	isl_int_divexact(c, m, bset1->eq[row][col]);
	isl_seq_scale(bset1->eq[row], bset1->eq[row], c, col+1);
	isl_int_divexact(c, m, bset2->eq[row][col]);
	isl_seq_scale(bset2->eq[row], bset2->eq[row], c, col+1);
	isl_int_clear(c);
	isl_int_clear(m);
}

/* Delete a given equality, moving all the following equalities one up.
 */
static void delete_row(struct isl_basic_set *bset, unsigned row)
{
	isl_int *t;
	int r;

	t = bset->eq[row];
	bset->n_eq--;
	for (r = row; r < bset->n_eq; ++r)
		bset->eq[r] = bset->eq[r+1];
	bset->eq[bset->n_eq] = t;
}

/* Make first row entries in column col of bset1 identical to
 * those of bset2, using the fact that entry bset1->eq[row][col]=a
 * is non-zero.  Initially, these elements of bset1 are all zero.
 * For each row i < row, we set
 *		A[i] = a * A[i] + B[i][col] * A[row]
 *		B[i] = a * B[i]
 * so that
 *		A[i][col] = B[i][col] = a * old(B[i][col])
 */
static void construct_column(
	struct isl_basic_set *bset1, struct isl_basic_set *bset2,
	unsigned row, unsigned col)
{
	int r;
	isl_int a;
	isl_int b;
	unsigned total;

	isl_int_init(a);
	isl_int_init(b);
	total = 1 + isl_basic_set_n_dim(bset1);
	for (r = 0; r < row; ++r) {
		if (isl_int_is_zero(bset2->eq[r][col]))
			continue;
		isl_int_gcd(b, bset2->eq[r][col], bset1->eq[row][col]);
		isl_int_divexact(a, bset1->eq[row][col], b);
		isl_int_divexact(b, bset2->eq[r][col], b);
		isl_seq_combine(bset1->eq[r], a, bset1->eq[r],
					      b, bset1->eq[row], total);
		isl_seq_scale(bset2->eq[r], bset2->eq[r], a, total);
	}
	isl_int_clear(a);
	isl_int_clear(b);
	delete_row(bset1, row);
}

/* Make first row entries in column col of bset1 identical to
 * those of bset2, using only these entries of the two matrices.
 * Let t be the last row with different entries.
 * For each row i < t, we set
 *	A[i] = (A[t][col]-B[t][col]) * A[i] + (B[i][col]-A[i][col) * A[t]
 *	B[i] = (A[t][col]-B[t][col]) * B[i] + (B[i][col]-A[i][col) * B[t]
 * so that
 *	A[i][col] = B[i][col] = old(A[t][col]*B[i][col]-A[i][col]*B[t][col])
 */
static int transform_column(
	struct isl_basic_set *bset1, struct isl_basic_set *bset2,
	unsigned row, unsigned col)
{
	int i, t;
	isl_int a, b, g;
	unsigned total;

	for (t = row-1; t >= 0; --t)
		if (isl_int_ne(bset1->eq[t][col], bset2->eq[t][col]))
			break;
	if (t < 0)
		return 0;

	total = 1 + isl_basic_set_n_dim(bset1);
	isl_int_init(a);
	isl_int_init(b);
	isl_int_init(g);
	isl_int_sub(b, bset1->eq[t][col], bset2->eq[t][col]);
	for (i = 0; i < t; ++i) {
		isl_int_sub(a, bset2->eq[i][col], bset1->eq[i][col]);
		isl_int_gcd(g, a, b);
		isl_int_divexact(a, a, g);
		isl_int_divexact(g, b, g);
		isl_seq_combine(bset1->eq[i], g, bset1->eq[i], a, bset1->eq[t],
				total);
		isl_seq_combine(bset2->eq[i], g, bset2->eq[i], a, bset2->eq[t],
				total);
	}
	isl_int_clear(a);
	isl_int_clear(b);
	isl_int_clear(g);
	delete_row(bset1, t);
	delete_row(bset2, t);
	return 1;
}

/* The implementation is based on Section 5.2 of Michael Karr,
 * "Affine Relationships Among Variables of a Program",
 * except that the echelon form we use starts from the last column
 * and that we are dealing with integer coefficients.
 */
static struct isl_basic_set *affine_hull(
	struct isl_basic_set *bset1, struct isl_basic_set *bset2)
{
	unsigned total;
	int col;
	int row;

	total = 1 + isl_basic_set_n_dim(bset1);

	row = 0;
	for (col = total-1; col >= 0; --col) {
		int is_zero1 = row >= bset1->n_eq ||
			isl_int_is_zero(bset1->eq[row][col]);
		int is_zero2 = row >= bset2->n_eq ||
			isl_int_is_zero(bset2->eq[row][col]);
		if (!is_zero1 && !is_zero2) {
			set_common_multiple(bset1, bset2, row, col);
			++row;
		} else if (!is_zero1 && is_zero2) {
			construct_column(bset1, bset2, row, col);
		} else if (is_zero1 && !is_zero2) {
			construct_column(bset2, bset1, row, col);
		} else {
			if (transform_column(bset1, bset2, row, col))
				--row;
		}
	}
	isl_basic_set_free(bset2);
	isl_assert(ctx, row == bset1->n_eq, goto error);
	bset1 = isl_basic_set_normalize_constraints(bset1);
	return bset1;
error:
	isl_basic_set_free(bset1);
	return NULL;
}

static struct isl_basic_set *isl_basic_set_from_vec(struct isl_vec *vec)
{
	int i;
	int k;
	struct isl_basic_set *bset = NULL;
	struct isl_ctx *ctx;
	unsigned dim;

	if (!vec)
		return NULL;
	ctx = vec->ctx;
	isl_assert(ctx, vec->size != 0, goto error);

	bset = isl_basic_set_alloc(ctx, 0, vec->size - 1, 0, vec->size - 1, 0);
	if (!bset)
		goto error;
	dim = isl_basic_set_n_dim(bset);
	for (i = dim - 1; i >= 0; --i) {
		k = isl_basic_set_alloc_equality(bset);
		if (k < 0)
			goto error;
		isl_seq_clr(bset->eq[k], 1 + dim);
		isl_int_neg(bset->eq[k][0], vec->el[1 + i]);
		isl_int_set(bset->eq[k][1 + i], vec->el[0]);
	}
	isl_vec_free(vec);

	return bset;
error:
	isl_basic_set_free(bset);
	isl_vec_free(vec);
	return NULL;
}

/* Find an integer point in "bset" that lies outside of the equality
 * "eq" e(x) = 0.
 * If "up" is true, look for a point satisfying e(x) - 1 >= 0.
 * Otherwise, look for a point satisfying -e(x) - 1 >= 0 (i.e., e(x) <= -1).
 * The point, if found, is returned as a singleton set.
 * If no point can be found, the empty set is returned.
 *
 * Before solving an ILP problem, we first check if simply
 * adding the normal of the constraint to one of the known
 * integer points in the basic set yields another point
 * inside the basic set.
 */
static struct isl_basic_set *outside_point(struct isl_ctx *ctx,
	struct isl_basic_set *bset, isl_int *eq, int up)
{
	struct isl_basic_set *slice = NULL;
	struct isl_vec *sample;
	struct isl_basic_set *point;
	unsigned dim;
	int k;

	dim = isl_basic_set_n_dim(bset);
	sample = isl_vec_alloc(ctx, 1 + dim);
	if (!sample)
		return NULL;
	isl_int_set_si(sample->block.data[0], 1);
	isl_seq_combine(sample->block.data + 1,
		ctx->one, bset->sample->block.data + 1,
		up ? ctx->one : ctx->negone, eq + 1, dim);
	if (isl_basic_set_contains(bset, sample))
		return isl_basic_set_from_vec(sample);
	isl_vec_free(sample);
	sample = NULL;

	slice = isl_basic_set_copy(bset);
	if (!slice)
		goto error;
	slice = isl_basic_set_cow(slice);
	slice = isl_basic_set_extend(slice, 0, dim, 0, 0, 1);
	k = isl_basic_set_alloc_inequality(slice);
	if (k < 0)
		goto error;
	if (up)
		isl_seq_cpy(slice->ineq[k], eq, 1 + dim);
	else
		isl_seq_neg(slice->ineq[k], eq, 1 + dim);
	isl_int_sub_ui(slice->ineq[k][0], slice->ineq[k][0], 1);

	sample = isl_basic_set_sample(slice);
	if (!sample)
		goto error;
	if (sample->size == 0) {
		isl_vec_free(sample);
		point = isl_basic_set_empty_like(bset);
	} else
		point = isl_basic_set_from_vec(sample);

	return point;
error:
	isl_basic_set_free(slice);
	return NULL;
}

struct isl_basic_set *isl_basic_set_recession_cone(struct isl_basic_set *bset)
{
	int i;

	bset = isl_basic_set_cow(bset);
	if (!bset)
		return NULL;
	isl_assert(bset->ctx, bset->n_div == 0, goto error);

	for (i = 0; i < bset->n_eq; ++i)
		isl_int_set_si(bset->eq[i][0], 0);

	for (i = 0; i < bset->n_ineq; ++i)
		isl_int_set_si(bset->ineq[i][0], 0);

	ISL_F_CLR(bset, ISL_BASIC_SET_NO_IMPLICIT);
	return isl_basic_set_implicit_equalities(bset);
error:
	isl_basic_set_free(bset);
	return NULL;
}

static struct isl_basic_set *shift(struct isl_basic_set *bset, isl_int *point)
{
	int i;
	unsigned dim;

	bset = isl_basic_set_cow(bset);
	if (!bset)
		return NULL;

	dim = isl_basic_set_n_dim(bset);
	for (i = 0; i < bset->n_eq; ++i) {
		isl_seq_inner_product(bset->eq[i]+1, point+1, dim,
					&bset->eq[i][0]);
		isl_int_neg(bset->eq[i][0], bset->eq[i][0]);
	}

	for (i = 0; i < bset->n_ineq; ++i) {
		isl_seq_inner_product(bset->ineq[i]+1, point+1, dim,
					&bset->ineq[i][0]);
		isl_int_neg(bset->ineq[i][0], bset->ineq[i][0]);
	}

	return bset;
}

/* Look for all equalities satisfied by the integer points in bset,
 * which is assume not to have any explicit equalities.
 *
 * The equalities are obtained by successively looking for
 * a point that is affinely independent of the points found so far.
 * In particular, for each equality satisfied by the points so far,
 * we check if there is any point on a hyperplane parallel to the
 * corresponding hyperplane shifted by at least one (in either direction).
 *
 * Before looking for any outside points, we first remove the equalities
 * that correspond to the affine hull of the recession cone.
 * These equalities will never be equalities over the whols basic set.
 */
static struct isl_basic_set *uset_affine_hull(struct isl_basic_set *bset)
{
	int i, j;
	struct isl_basic_set *hull = NULL;
	struct isl_vec *sample;
	struct isl_ctx *ctx;
	unsigned dim;

	if (isl_basic_set_is_empty(bset))
		return bset;

	ctx = bset->ctx;
	sample = isl_basic_set_sample(isl_basic_set_copy(bset));
	if (!sample)
		goto error;
	if (sample->size == 0) {
		isl_vec_free(sample);
		hull = isl_basic_set_empty_like(bset);
		isl_basic_set_free(bset);
		return hull;
	} else
		hull = isl_basic_set_from_vec(sample);

	if (hull->n_eq > 0) {
		struct isl_basic_set *cone;
		cone = isl_basic_set_recession_cone(isl_basic_set_copy(bset));
		isl_basic_set_free_inequality(cone, cone->n_ineq);
		cone = isl_basic_set_normalize_constraints(cone);
		cone = shift(cone, bset->sample->block.data);
		hull = affine_hull(hull, cone);
	}

	dim = isl_basic_set_n_dim(bset);
	for (i = 0; i < dim; ++i) {
		struct isl_basic_set *point;
		for (j = 0; j < hull->n_eq; ++j) {
			point = outside_point(ctx, bset, hull->eq[j], 1);
			if (!point)
				goto error;
			if (!ISL_F_ISSET(point, ISL_BASIC_SET_EMPTY))
				break;
			isl_basic_set_free(point);
			point = outside_point(ctx, bset, hull->eq[j], 0);
			if (!point)
				goto error;
			if (!ISL_F_ISSET(point, ISL_BASIC_SET_EMPTY))
				break;
			isl_basic_set_free(point);
		}
		if (j == hull->n_eq)
			break;
		hull = affine_hull(hull, point);
	}
	isl_basic_set_free(bset);

	return hull;
error:
	isl_basic_set_free(bset);
	isl_basic_set_free(hull);
	return NULL;
}

/* Look for all equalities satisfied by the integer points in bmap
 * that are independent of the equalities already explicitly available
 * in bmap.
 *
 * We first remove all equalities already explicitly available,
 * then look for additional equalities in the reduced space
 * and then transform the result to the original space.
 * The original equalities are _not_ added to this set.  This is
 * the responsibility of the calling function.
 * The resulting basic set has all meaning about the dimensions removed.
 * In particular, dimensions that correspond to existential variables
 * in bmap and that are found to be fixed are not removed.
 */
static struct isl_basic_set *equalities_in_underlying_set(
						struct isl_basic_map *bmap)
{
	struct isl_mat *T2 = NULL;
	struct isl_basic_set *bset = NULL;
	struct isl_basic_set *hull = NULL;

	bset = isl_basic_map_underlying_set(bmap);
	bset = isl_basic_set_remove_equalities(bset, NULL, &T2);
	if (!bset)
		goto error;

	hull = uset_affine_hull(bset);
	if (T2)
		hull = isl_basic_set_preimage(hull, T2);

	return hull;
error:
	isl_mat_free(T2);
	isl_basic_set_free(bset);
	isl_basic_set_free(hull);
	return NULL;
}

/* Detect and make explicit all equalities satisfied by the (integer)
 * points in bmap.
 */
struct isl_basic_map *isl_basic_map_detect_equalities(
						struct isl_basic_map *bmap)
{
	int i, j;
	struct isl_basic_set *hull = NULL;

	if (!bmap)
		return NULL;
	if (bmap->n_ineq == 0)
		return bmap;
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_EMPTY))
		return bmap;
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_ALL_EQUALITIES))
		return bmap;
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_RATIONAL))
		return isl_basic_map_implicit_equalities(bmap);

	hull = equalities_in_underlying_set(isl_basic_map_copy(bmap));
	if (!hull)
		goto error;
	if (ISL_F_ISSET(hull, ISL_BASIC_SET_EMPTY)) {
		isl_basic_set_free(hull);
		return isl_basic_map_set_to_empty(bmap);
	}
	bmap = isl_basic_map_extend_dim(bmap, isl_dim_copy(bmap->dim), 0,
					hull->n_eq, 0);
	for (i = 0; i < hull->n_eq; ++i) {
		j = isl_basic_map_alloc_equality(bmap);
		if (j < 0)
			goto error;
		isl_seq_cpy(bmap->eq[j], hull->eq[i],
				1 + isl_basic_set_total_dim(hull));
	}
	isl_basic_set_free(hull);
	ISL_F_SET(bmap, ISL_BASIC_MAP_NO_IMPLICIT | ISL_BASIC_MAP_ALL_EQUALITIES);
	bmap = isl_basic_map_simplify(bmap);
	return isl_basic_map_finalize(bmap);
error:
	isl_basic_set_free(hull);
	isl_basic_map_free(bmap);
	return NULL;
}

struct isl_map *isl_map_detect_equalities(struct isl_map *map)
{
	struct isl_basic_map *bmap;
	int i;

	if (!map)
		return NULL;

	for (i = 0; i < map->n; ++i) {
		bmap = isl_basic_map_copy(map->p[i]);
		bmap = isl_basic_map_detect_equalities(bmap);
		if (!bmap)
			goto error;
		isl_basic_map_free(map->p[i]);
		map->p[i] = bmap;
	}

	return map;
error:
	isl_map_free(map);
	return NULL;
}

/* After computing the rational affine hull (by detecting the implicit
 * equalities), we compute the additional equalities satisfied by
 * the integer points (if any) and add the original equalities back in.
 */
struct isl_basic_map *isl_basic_map_affine_hull(struct isl_basic_map *bmap)
{
	struct isl_basic_set *hull = NULL;

	bmap = isl_basic_map_detect_equalities(bmap);
	bmap = isl_basic_map_cow(bmap);
	isl_basic_map_free_inequality(bmap, bmap->n_ineq);
	return bmap;
}

struct isl_basic_set *isl_basic_set_affine_hull(struct isl_basic_set *bset)
{
	return (struct isl_basic_set *)
		isl_basic_map_affine_hull((struct isl_basic_map *)bset);
}

struct isl_basic_map *isl_map_affine_hull(struct isl_map *map)
{
	int i;
	struct isl_basic_map *model = NULL;
	struct isl_basic_map *hull = NULL;
	struct isl_set *set;

	if (!map)
		return NULL;

	if (map->n == 0) {
		hull = isl_basic_map_empty_like_map(map);
		isl_map_free(map);
		return hull;
	}

	map = isl_map_detect_equalities(map);
	map = isl_map_align_divs(map);
	if (!map)
		return NULL;
	model = isl_basic_map_copy(map->p[0]);
	set = isl_map_underlying_set(map);
	set = isl_set_cow(set);
	if (!set)
		goto error;

	for (i = 0; i < set->n; ++i) {
		set->p[i] = isl_basic_set_cow(set->p[i]);
		set->p[i] = isl_basic_set_affine_hull(set->p[i]);
		set->p[i] = isl_basic_set_gauss(set->p[i], NULL);
		if (!set->p[i])
			goto error;
	}
	set = isl_set_remove_empty_parts(set);
	if (set->n == 0) {
		hull = isl_basic_map_empty_like(model);
		isl_basic_map_free(model);
	} else {
		struct isl_basic_set *bset;
		while (set->n > 1) {
			set->p[0] = affine_hull(set->p[0], set->p[--set->n]);
			if (!set->p[0])
				goto error;
		}
		bset = isl_basic_set_copy(set->p[0]);
		hull = isl_basic_map_overlying_set(bset, model);
	}
	isl_set_free(set);
	hull = isl_basic_map_simplify(hull);
	return isl_basic_map_finalize(hull);
error:
	isl_basic_map_free(model);
	isl_set_free(set);
	return NULL;
}

struct isl_basic_set *isl_set_affine_hull(struct isl_set *set)
{
	return (struct isl_basic_set *)
		isl_map_affine_hull((struct isl_map *)set);
}
