#include <ctype.h>
#include <stdio.h>
#include <isl_set.h>
#include "isl_dim.h"
#include "isl_map_private.h"
#include "isl_input_omega.h"

static char *next_line(FILE *input, char *line, unsigned len)
{
	char *p;

	do {
		if (!(p = fgets(line, len, input)))
			return NULL;
		while (isspace(*p) && *p != '\n')
			++p;
	} while (*p == '#' || *p == '\n');

	return p;
}

static struct isl_basic_set *isl_basic_set_read_from_file_polylib(
		struct isl_ctx *ctx, FILE *input, unsigned nparam)
{
	struct isl_basic_set *bset = NULL;
	int i, j;
	unsigned n_row, n_col;
	unsigned dim;
	char line[1024];
	char val[1024];
	char *p;

	isl_assert(ctx, next_line(input, line, sizeof(line)), return NULL);
	isl_assert(ctx, sscanf(line, "%u %u", &n_row, &n_col) == 2, return NULL);
	isl_assert(ctx, n_col >= 2, return NULL);
	dim = n_col - 2 - nparam;
	bset = isl_basic_set_alloc(ctx, nparam, dim, 0, n_row, n_row);
	if (!bset)
		return NULL;
	for (i = 0; i < n_row; ++i) {
		int type;
		int offset;
		int n;
		int k;
		isl_int *c;

		p = next_line(input, line, sizeof(line));
		isl_assert(ctx, p, goto error);
		n = sscanf(p, "%u%n", &type, &offset);
		isl_assert(ctx, n != 0, goto error);
		p += offset;
		isl_assert(ctx, type == 0 || type == 1, goto error);
		if (type == 0) {
			k = isl_basic_set_alloc_equality(bset);
			c = bset->eq[k];
		} else {
			k = isl_basic_set_alloc_inequality(bset);
			c = bset->ineq[k];
		}
		isl_assert(ctx, k >= 0, goto error);
		for (j = 0; j < dim; ++j) {
			n = sscanf(p, "%s%n", val, &offset);
			isl_assert(ctx, n != 0, goto error);
			isl_int_read(c[1+nparam+j], val);
			p += offset;
		}
		for (j = 0; j < nparam; ++j) {
			n = sscanf(p, "%s%n", val, &offset);
			isl_assert(ctx, n != 0, goto error);
			isl_int_read(c[1+j], val);
			p += offset;
		}
		n = sscanf(p, "%s%n", val, &offset);
		isl_assert(ctx, n != 0, goto error);
		isl_int_read(c[0], val);
	}
	bset = isl_basic_set_simplify(bset);
	bset = isl_basic_set_finalize(bset);
	return bset;
error:
	isl_basic_set_free(bset);
	return NULL;
}

static struct isl_set *isl_set_read_from_file_polylib(
		struct isl_ctx *ctx, FILE *input, unsigned nparam)
{
	struct isl_set *set = NULL;
	char line[1024];
	int i;
	unsigned n;

	isl_assert(ctx, next_line(input, line, sizeof(line)), return NULL);
	isl_assert(ctx, sscanf(line, "%u", &n) == 1, return NULL);

	set = isl_set_alloc(ctx, nparam, 0, n, 0);
	if (!set)
		return NULL;
	if (n == 0)
		return set;
	for (i = 0; i < n; ++i) {
		set->p[i] = isl_basic_set_read_from_file_polylib(ctx, input,
								 nparam);
		if (!set->p[i])
			goto error;
	}
	set->n = n;
	isl_dim_free(set->dim);
	set->dim = isl_dim_copy(set->p[0]->dim);
	for (i = 1; i < n; ++i)
		isl_assert(ctx, isl_dim_equal(set->dim, set->p[i]->dim),
				goto error);
	return set;
error:
	isl_set_free(set);
	return NULL;
}

struct isl_basic_set *isl_basic_set_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format)
{
	if (input_format == ISL_FORMAT_POLYLIB)
		return isl_basic_set_read_from_file_polylib(ctx, input, nparam);
	else if (input_format == ISL_FORMAT_OMEGA) {
		isl_assert(ctx, nparam == 0, return NULL);
		return isl_basic_set_read_from_file_omega(ctx, input);
	} else
		isl_assert(ctx, 0, return NULL);
}

struct isl_basic_set *isl_basic_set_read_from_str(struct isl_ctx *ctx,
		const char *str, unsigned nparam, unsigned input_format)
{
	if (input_format == ISL_FORMAT_OMEGA) {
		isl_assert(ctx, nparam == 0, return NULL);
		return isl_basic_set_read_from_str_omega(ctx, str);
	} else
		isl_assert(ctx, 0, return NULL);
}

struct isl_basic_map *isl_basic_map_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format)
{
	if (input_format == ISL_FORMAT_OMEGA)
		return isl_basic_map_read_from_file_omega(ctx, input);
	else
		isl_assert(ctx, 0, return NULL);
}

struct isl_set *isl_set_read_from_file(struct isl_ctx *ctx,
		FILE *input, unsigned nparam, unsigned input_format)
{
	if (input_format == ISL_FORMAT_POLYLIB)
		return isl_set_read_from_file_polylib(ctx, input, nparam);
	else
		isl_assert(ctx, 0, return NULL);
}
