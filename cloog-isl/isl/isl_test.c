#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <isl_ctx.h>
#include <isl_set.h>
#include <isl_constraint.h>

static char *srcdir;

void test_read(struct isl_ctx *ctx)
{
	char filename[PATH_MAX];
	FILE *input;
	int n;
	struct isl_basic_set *bset1, *bset2;
	const char *str = "{[y]: Exists ( alpha : 2alpha = y)}";

	n = snprintf(filename, sizeof(filename),
			"%s/test_inputs/set.omega", srcdir);
	assert(n < sizeof(filename));
	input = fopen(filename, "r");
	assert(input);

	bset1 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_OMEGA);
	bset2 = isl_basic_set_read_from_str(ctx, str, 0, ISL_FORMAT_OMEGA);

	assert(isl_basic_set_is_equal(bset1, bset2) == 1);

	isl_basic_set_free(bset1);
	isl_basic_set_free(bset2);

	fclose(input);
}

/* Construct the basic set { [i] : 5 <= i <= N } */
void test_construction(struct isl_ctx *ctx)
{
	isl_int v;
	struct isl_dim *dim;
	struct isl_basic_set *bset;
	struct isl_constraint *c;

	isl_int_init(v);

	dim = isl_dim_set_alloc(ctx, 1, 1);
	bset = isl_basic_set_universe(dim);

	c = isl_inequality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_param, 0, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_inequality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -5);
	isl_constraint_set_constant(c, v);
	bset = isl_basic_set_add_constraint(bset, c);

	isl_basic_set_free(bset);

	isl_int_clear(v);
}

void test_div(struct isl_ctx *ctx)
{
	isl_int v;
	int pos;
	struct isl_dim *dim;
	struct isl_div *div;
	struct isl_basic_set *bset;
	struct isl_constraint *c;

	isl_int_init(v);

	/* test 1 */
	dim = isl_dim_set_alloc(ctx, 0, 1);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 2 */
	dim = isl_dim_set_alloc(ctx, 0, 1);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 3 */
	dim = isl_dim_set_alloc(ctx, 0, 1);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -3);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 4);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 4 */
	dim = isl_dim_set_alloc(ctx, 0, 1);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 2);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_constant(c, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 6);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(isl_basic_set_is_empty(bset));
	isl_basic_set_free(bset);

	/* test 5 */
	dim = isl_dim_set_alloc(ctx, 0, 2);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 1, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 0);
	isl_basic_set_free(bset);

	/* test 6 */
	dim = isl_dim_set_alloc(ctx, 0, 2);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 6);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 1, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 7 */
	/* This test is a bit tricky.  We set up an equality
	 *		a + 3b + 3c = 6 e0
	 * Normalization of divs creates _two_ divs
	 *		a = 3 e0
	 *		c - b - e0 = 2 e1
	 * Afterwards e0 is removed again because it has coefficient -1
	 * and we end up with the original equality and div again.
	 * Perhaps we can avoid the introduction of this temporary div.
	 */
	dim = isl_dim_set_alloc(ctx, 0, 3);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 1, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 2, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 6);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 8 */
	dim = isl_dim_set_alloc(ctx, 0, 4);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 1, v);
	isl_int_set_si(v, -3);
	isl_constraint_set_coefficient(c, isl_dim_set, 3, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 6);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 2, v);
	isl_int_set_si(v, 1);
	isl_constraint_set_constant(c, v);
	bset = isl_basic_set_add_constraint(bset, c);

	assert(bset->n_div == 1);
	isl_basic_set_free(bset);

	/* test 9 */
	dim = isl_dim_set_alloc(ctx, 0, 2);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 1, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, -2);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, -1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, 3);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	isl_int_set_si(v, 2);
	isl_constraint_set_constant(c, v);
	bset = isl_basic_set_add_constraint(bset, c);

	bset = isl_basic_set_fix_si(bset, isl_dim_set, 0, 2);

	assert(!isl_basic_set_is_empty(bset));

	isl_basic_set_free(bset);

	/* test 10 */
	dim = isl_dim_set_alloc(ctx, 0, 2);
	bset = isl_basic_set_universe(dim);

	c = isl_equality_alloc(isl_dim_copy(bset->dim));
	isl_int_set_si(v, 1);
	isl_constraint_set_coefficient(c, isl_dim_set, 0, v);
	div = isl_div_alloc(isl_dim_copy(bset->dim));
	c = isl_constraint_add_div(c, div, &pos);
	isl_int_set_si(v, -2);
	isl_constraint_set_coefficient(c, isl_dim_div, pos, v);
	bset = isl_basic_set_add_constraint(bset, c);

	bset = isl_basic_set_fix_si(bset, isl_dim_set, 0, 2);

	isl_basic_set_free(bset);

	isl_int_clear(v);
}

void test_application_case(struct isl_ctx *ctx, const char *name)
{
	char filename[PATH_MAX];
	FILE *input;
	int n;
	struct isl_basic_set *bset1, *bset2;
	struct isl_basic_map *bmap;

	n = snprintf(filename, sizeof(filename),
			"%s/test_inputs/%s.omega", srcdir, name);
	assert(n < sizeof(filename));
	input = fopen(filename, "r");
	assert(input);

	bset1 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_OMEGA);
	bmap = isl_basic_map_read_from_file(ctx, input, 0, ISL_FORMAT_OMEGA);

	bset1 = isl_basic_set_apply(bset1, bmap);

	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_OMEGA);

	assert(isl_basic_set_is_equal(bset1, bset2) == 1);

	isl_basic_set_free(bset1);
	isl_basic_set_free(bset2);

	fclose(input);
}

void test_application(struct isl_ctx *ctx)
{
	test_application_case(ctx, "application");
	test_application_case(ctx, "application2");
}

void test_affine_hull_case(struct isl_ctx *ctx, const char *name)
{
	char filename[PATH_MAX];
	FILE *input;
	int n;
	struct isl_basic_set *bset1, *bset2;

	n = snprintf(filename, sizeof(filename),
			"%s/test_inputs/%s.polylib", srcdir, name);
	assert(n < sizeof(filename));
	input = fopen(filename, "r");
	assert(input);

	bset1 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);
	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);

	bset1 = isl_basic_set_affine_hull(bset1);

	assert(isl_basic_set_is_equal(bset1, bset2) == 1);

	isl_basic_set_free(bset1);
	isl_basic_set_free(bset2);

	fclose(input);
}

void test_affine_hull(struct isl_ctx *ctx)
{
	test_affine_hull_case(ctx, "affine2");
	test_affine_hull_case(ctx, "affine");
	test_affine_hull_case(ctx, "affine3");
}

void test_convex_hull_case(struct isl_ctx *ctx, const char *name)
{
	char filename[PATH_MAX];
	FILE *input;
	int n;
	struct isl_basic_set *bset1, *bset2;
	struct isl_set *set;

	n = snprintf(filename, sizeof(filename),
			"%s/test_inputs/%s.polylib", srcdir, name);
	assert(n < sizeof(filename));
	input = fopen(filename, "r");
	assert(input);

	bset1 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);
	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);

	set = isl_basic_set_union(bset1, bset2);
	bset1 = isl_set_convex_hull(set);

	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);

	assert(isl_basic_set_is_equal(bset1, bset2) == 1);

	isl_basic_set_free(bset1);
	isl_basic_set_free(bset2);

	fclose(input);
}

void test_convex_hull(struct isl_ctx *ctx)
{
	test_convex_hull_case(ctx, "convex0");
	test_convex_hull_case(ctx, "convex1");
	test_convex_hull_case(ctx, "convex2");
	test_convex_hull_case(ctx, "convex3");
	test_convex_hull_case(ctx, "convex4");
	test_convex_hull_case(ctx, "convex5");
	test_convex_hull_case(ctx, "convex6");
	test_convex_hull_case(ctx, "convex7");
	test_convex_hull_case(ctx, "convex8");
	test_convex_hull_case(ctx, "convex9");
	test_convex_hull_case(ctx, "convex10");
	test_convex_hull_case(ctx, "convex11");
	test_convex_hull_case(ctx, "convex12");
	test_convex_hull_case(ctx, "convex13");
	test_convex_hull_case(ctx, "convex14");
}

void test_gist_case(struct isl_ctx *ctx, const char *name)
{
	char filename[PATH_MAX];
	FILE *input;
	int n;
	struct isl_basic_set *bset1, *bset2;
	struct isl_set *set;

	n = snprintf(filename, sizeof(filename),
			"%s/test_inputs/%s.polylib", srcdir, name);
	assert(n < sizeof(filename));
	input = fopen(filename, "r");
	assert(input);

	bset1 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);
	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);

	bset1 = isl_basic_set_gist(bset1, bset2);

	bset2 = isl_basic_set_read_from_file(ctx, input, 0, ISL_FORMAT_POLYLIB);

	assert(isl_basic_set_is_equal(bset1, bset2) == 1);

	isl_basic_set_free(bset1);
	isl_basic_set_free(bset2);

	fclose(input);
}

void test_gist(struct isl_ctx *ctx)
{
	test_gist_case(ctx, "gist1");
}

static struct isl_set *s_union(struct isl_ctx *ctx,
	const char *s1, const char *s2)
{
	struct isl_basic_set *bset1;
	struct isl_basic_set *bset2;
	struct isl_set *set1, *set2;

	bset1 = isl_basic_set_read_from_str(ctx, s1, 0, ISL_FORMAT_OMEGA);
	bset2 = isl_basic_set_read_from_str(ctx, s2, 0, ISL_FORMAT_OMEGA);
	set1 = isl_set_from_basic_set(bset1);
	set2 = isl_set_from_basic_set(bset2);
	return isl_set_union(set1, set2);
}

void test_coalesce(struct isl_ctx *ctx)
{
	struct isl_set *set;

	set = s_union(ctx, "{[x,y]: x >= 0 & x <= 10 & y >= 0 & y <= 10}",
			   "{[x,y]: y >= x & x >= 2 & 5 >= y}");
	set = isl_set_coalesce(set);
	assert(set->n == 1);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & 2x + y <= 30 & y <= 10 & x >= 0}",
		       "{[x,y]: x + y >= 10 & y <= x & x + y <= 20 & y >= 0}");
	set = isl_set_coalesce(set);
	assert(set->n == 1);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & 2x + y <= 30 & y <= 10 & x >= 0}",
		       "{[x,y]: x + y >= 10 & y <= x & x + y <= 19 & y >= 0}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x >= 6 & x <= 10 & y <= x}");
	set = isl_set_coalesce(set);
	assert(set->n == 1);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x >= 7 & x <= 10 & y <= x}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x >= 6 & x <= 10 & y + 1 <= x}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x = 6 & y <= 6}");
	set = isl_set_coalesce(set);
	assert(set->n == 1);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x = 7 & y <= 6}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x = 6 & y <= 5}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);

	set = s_union(ctx, "{[x,y]: y >= 0 & x <= 5 & y <= x}",
		       "{[x,y]: y >= 0 & x = 6 & y <= 7}");
	set = isl_set_coalesce(set);
	assert(set->n == 2);
	isl_set_free(set);
}

int main()
{
	struct isl_ctx *ctx;

	srcdir = getenv("srcdir");

	ctx = isl_ctx_alloc();
	test_read(ctx);
	test_construction(ctx);
	test_div(ctx);
	test_application(ctx);
	test_affine_hull(ctx);
	test_convex_hull(ctx);
	test_gist(ctx);
	test_coalesce(ctx);
	isl_ctx_free(ctx);
	return 0;
}
