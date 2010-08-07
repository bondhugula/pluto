#include "isl_dim.h"
#include "isl_name.h"

struct isl_dim *isl_dim_alloc(struct isl_ctx *ctx,
			unsigned nparam, unsigned n_in, unsigned n_out)
{
	struct isl_dim *dim;

	dim = isl_alloc_type(ctx, struct isl_dim);
	if (!dim)
		return NULL;

	dim->ctx = ctx;
	isl_ctx_ref(ctx);
	dim->ref = 1;
	dim->nparam = nparam;
	dim->n_in = n_in;
	dim->n_out = n_out;

	dim->n_name = 0;
	dim->names = NULL;

	return dim;
}

struct isl_dim *isl_dim_set_alloc(struct isl_ctx *ctx,
			unsigned nparam, unsigned dim)
{
	return isl_dim_alloc(ctx, nparam, 0, dim);
}

static unsigned global_pos(struct isl_dim *dim,
				 enum isl_dim_type type, unsigned pos)
{
	struct isl_ctx *ctx = dim->ctx;

	switch (type) {
	case isl_dim_param:
		isl_assert(ctx, pos < dim->nparam, return isl_dim_total(dim));
		return pos;
	case isl_dim_in:
		isl_assert(ctx, pos < dim->n_in, return isl_dim_total(dim));
		return pos + dim->nparam;
	case isl_dim_out:
		isl_assert(ctx, pos < dim->n_out, return isl_dim_total(dim));
		return pos + dim->nparam + dim->n_in;
	default:
		isl_assert(ctx, 0, goto error);
	}
	return isl_dim_total(dim);
}

static struct isl_dim *set_name(struct isl_dim *dim,
				 enum isl_dim_type type, unsigned pos,
				 struct isl_name *name)
{
	struct isl_ctx *ctx = dim->ctx;
	dim = isl_dim_cow(dim);

	if (!dim)
		goto error;

	pos = global_pos(dim, type, pos);
	isl_assert(ctx, pos != isl_dim_total(dim), goto error);

	if (pos >= dim->n_name) {
		if (!name)
			return dim;
		if (!dim->names) {
			dim->names = isl_calloc_array(dim->ctx,
					struct isl_name *, isl_dim_total(dim));
			if (!dim->names)
				goto error;
		} else {
			int i;
			dim->names = isl_realloc_array(dim->ctx, dim->names,
					struct isl_name *, isl_dim_total(dim));
			if (!dim->names)
				goto error;
			for (i = dim->n_name; i < isl_dim_total(dim); ++i)
				dim->names[i] = NULL;
		}
		dim->n_name = isl_dim_total(dim);
	}

	dim->names[pos] = name;

	return dim;
error:
	isl_name_free(ctx, name);
	isl_dim_free(dim);
	return NULL;
}

static struct isl_name *get_name(struct isl_dim *dim,
				 enum isl_dim_type type, unsigned pos)
{
	if (!dim)
		return NULL;

	pos = global_pos(dim, type, pos);
	if (pos == isl_dim_total(dim))
		return NULL;
	if (pos >= dim->n_name)
		return NULL;
	return dim->names[pos];
}

static unsigned offset(struct isl_dim *dim, enum isl_dim_type type)
{
	switch (type) {
	case isl_dim_param:	return 0;
	case isl_dim_in:	return dim->nparam;
	case isl_dim_out:	return dim->nparam + dim->n_in;
	}
}

static unsigned n(struct isl_dim *dim, enum isl_dim_type type)
{
	switch (type) {
	case isl_dim_param:	return dim->nparam;
	case isl_dim_in:	return dim->n_in;
	case isl_dim_out:	return dim->n_out;
	}
}

unsigned isl_dim_size(struct isl_dim *dim, enum isl_dim_type type)
{
	return n(dim, type);
}

static struct isl_dim *copy_names(struct isl_dim *dst,
	enum isl_dim_type dst_type, unsigned offset, struct isl_dim *src,
	enum isl_dim_type src_type)
{
	int i;
	struct isl_name *name;

	for (i = 0; i < n(src, src_type); ++i) {
		name = get_name(src, src_type, i);
		if (!name)
			continue;
		dst = set_name(dst, dst_type, offset + i,
					isl_name_copy(dst->ctx, name));
		if (!dst)
			return NULL;
	}
	return dst;
}

struct isl_dim *isl_dim_dup(struct isl_dim *dim)
{
	struct isl_dim *dup;
	dup = isl_dim_alloc(dim->ctx, dim->nparam, dim->n_in, dim->n_out);
	if (!dim->names)
		return dup;
	dup = copy_names(dup, isl_dim_param, 0, dim, isl_dim_param);
	dup = copy_names(dup, isl_dim_in, 0, dim, isl_dim_in);
	dup = copy_names(dup, isl_dim_out, 0, dim, isl_dim_out);
	return dup;
}

struct isl_dim *isl_dim_cow(struct isl_dim *dim)
{
	if (!dim)
		return NULL;

	if (dim->ref == 1)
		return dim;
	dim->ref--;
	return isl_dim_dup(dim);
}

struct isl_dim *isl_dim_copy(struct isl_dim *dim)
{
	if (!dim)
		return NULL;

	dim->ref++;
	return dim;
}

void isl_dim_free(struct isl_dim *dim)
{
	int i;

	if (!dim)
		return;

	if (--dim->ref > 0)
		return;

	for (i = 0; i < dim->n_name; ++i)
		isl_name_free(dim->ctx, dim->names[i]);
	free(dim->names);
	isl_ctx_deref(dim->ctx);
	
	free(dim);
}

struct isl_dim *isl_dim_set_name(struct isl_dim *dim,
				 enum isl_dim_type type, unsigned pos,
				 const char *s)
{
	struct isl_name *name;
	if (!dim)
		return NULL;
	name = isl_name_get(dim->ctx, s);
	if (!name)
		goto error;
	return set_name(dim, type, pos, name);
error:
	isl_dim_free(dim);
	return NULL;
}

const char *isl_dim_get_name(struct isl_dim *dim,
				 enum isl_dim_type type, unsigned pos)
{
	struct isl_name *name = get_name(dim, type, pos);
	return name ? name->name : NULL;
}

static int match(struct isl_dim *dim1, enum isl_dim_type dim1_type,
		struct isl_dim *dim2, enum isl_dim_type dim2_type)
{
	int i;

	if (n(dim1, dim1_type) != n(dim2, dim2_type))
		return 0;

	if (!dim1->names && !dim2->names)
		return 1;

	for (i = 0; i < n(dim1, dim1_type); ++i) {
		if (get_name(dim1, dim1_type, i) !=
		    get_name(dim2, dim2_type, i))
			return 0;
	}
	return 1;
}

int isl_dim_match(struct isl_dim *dim1, enum isl_dim_type dim1_type,
		struct isl_dim *dim2, enum isl_dim_type dim2_type)
{
	return match(dim1, dim1_type, dim2, dim2_type);
}

static void get_names(struct isl_dim *dim, enum isl_dim_type type,
	unsigned first, unsigned n, struct isl_name **names)
{
	int i;

	for (i = 0; i < n ; ++i)
		names[i] = get_name(dim, type, first+i);
}

struct isl_dim *isl_dim_extend(struct isl_dim *dim,
			unsigned nparam, unsigned n_in, unsigned n_out)
{
	struct isl_name **names = NULL;

	if (!dim)
		return NULL;
	if (dim->nparam == nparam && dim->n_in == n_in && dim->n_out == n_out)
		return dim;

	isl_assert(dim->ctx, dim->nparam <= nparam, goto error);
	isl_assert(dim->ctx, dim->n_in <= n_in, goto error);
	isl_assert(dim->ctx, dim->n_out <= n_out, goto error);

	dim = isl_dim_cow(dim);

	if (dim->names) {
		names = isl_calloc_array(dim->ctx, struct isl_name *,
					 nparam + n_in + n_out);
		if (!names)
			goto error;
		get_names(dim, isl_dim_param, 0, dim->nparam, names);
		get_names(dim, isl_dim_in, 0, dim->n_in, names + nparam);
		get_names(dim, isl_dim_out, 0, dim->n_out,
				names + nparam + n_in);
		free(dim->names);
		dim->names = names;
		dim->n_name = nparam + n_in + n_out;
	}
	dim->nparam = nparam;
	dim->n_in = n_in;
	dim->n_out = n_out;

	return dim;
error:
	free(names);
	isl_dim_free(dim);
	return NULL;
}

struct isl_dim *isl_dim_add(struct isl_dim *dim, enum isl_dim_type type,
	unsigned n)
{
	switch (type) {
	case isl_dim_param:
		return isl_dim_extend(dim,
					dim->nparam + n, dim->n_in, dim->n_out);
	case isl_dim_in:
		return isl_dim_extend(dim,
					dim->nparam, dim->n_in + n, dim->n_out);
	case isl_dim_out:
		return isl_dim_extend(dim,
					dim->nparam, dim->n_in, dim->n_out + n);
	}
	return dim;
}

struct isl_dim *isl_dim_join(struct isl_dim *left, struct isl_dim *right)
{
	struct isl_dim *dim;

	if (!left || !right)
		goto error;

	isl_assert(left->ctx, match(left, isl_dim_param, right, isl_dim_param),
			goto error);
	isl_assert(left->ctx, match(left, isl_dim_out, right, isl_dim_in),
			goto error);

	dim = isl_dim_alloc(left->ctx, left->nparam, left->n_in, right->n_out);
	if (!dim)
		goto error;

	dim = copy_names(dim, isl_dim_param, 0, left, isl_dim_param);
	dim = copy_names(dim, isl_dim_in, 0, left, isl_dim_in);
	dim = copy_names(dim, isl_dim_out, 0, right, isl_dim_out);

	isl_dim_free(left);
	isl_dim_free(right);

	return dim;
error:
	isl_dim_free(left);
	isl_dim_free(right);
	return NULL;
}

struct isl_dim *isl_dim_product(struct isl_dim *left, struct isl_dim *right)
{
	struct isl_dim *dim;

	if (!left || !right)
		goto error;

	isl_assert(left->ctx, match(left, isl_dim_param, right, isl_dim_param),
			goto error);

	dim = isl_dim_alloc(left->ctx, left->nparam,
			left->n_in + right->n_in, left->n_out + right->n_out);
	if (!dim)
		goto error;

	dim = copy_names(dim, isl_dim_param, 0, left, isl_dim_param);
	dim = copy_names(dim, isl_dim_in, 0, left, isl_dim_in);
	dim = copy_names(dim, isl_dim_in, left->n_in, right, isl_dim_in);
	dim = copy_names(dim, isl_dim_out, 0, left, isl_dim_out);
	dim = copy_names(dim, isl_dim_out, left->n_out, right, isl_dim_out);

	isl_dim_free(left);
	isl_dim_free(right);

	return dim;
error:
	isl_dim_free(left);
	isl_dim_free(right);
	return NULL;
}

struct isl_dim *isl_dim_map(struct isl_dim *dim)
{
	struct isl_name **names = NULL;

	if (!dim)
		return NULL;
	isl_assert(dim->ctx, dim->n_in == 0, goto error);
	if (dim->n_out == 0)
		return dim;
	dim = isl_dim_cow(dim);
	if (!dim)
		return NULL;
	if (dim->names) {
		names = isl_calloc_array(dim->ctx, struct isl_name *,
					dim->nparam + dim->n_out + dim->n_out);
		if (!names)
			goto error;
		get_names(dim, isl_dim_param, 0, dim->nparam, names);
		get_names(dim, isl_dim_out, 0, dim->n_out, names + dim->nparam);
	}
	dim->n_in = dim->n_out;
	if (names) {
		copy_names(dim, isl_dim_out, 0, dim, isl_dim_in);
		free(dim->names);
		dim->names = names;
		dim->n_name = dim->nparam + dim->n_out + dim->n_out;
	}
	return dim;
error:
	isl_dim_free(dim);
	return NULL;
}

static struct isl_dim *set_names(struct isl_dim *dim, enum isl_dim_type type,
	unsigned first, unsigned n, struct isl_name **names)
{
	int i;

	for (i = 0; i < n ; ++i)
		dim = set_name(dim, type, first+i, names[i]);

	return dim;
}

struct isl_dim *isl_dim_reverse(struct isl_dim *dim)
{
	unsigned t;
	struct isl_name **names = NULL;

	if (!dim)
		return NULL;
	if (match(dim, isl_dim_in, dim, isl_dim_out))
		return dim;

	dim = isl_dim_cow(dim);
	if (!dim)
		return NULL;

	if (dim->names) {
		names = isl_alloc_array(dim->ctx, struct isl_name *,
					dim->n_in + dim->n_out);
		if (!names)
			goto error;
		get_names(dim, isl_dim_in, 0, dim->n_in, names);
		get_names(dim, isl_dim_out, 0, dim->n_out, names + dim->n_in);
	}

	t = dim->n_in;
	dim->n_in = dim->n_out;
	dim->n_out = t;

	if (dim->names) {
		dim = set_names(dim, isl_dim_out, 0, dim->n_out, names);
		dim = set_names(dim, isl_dim_in, 0, dim->n_in, names + dim->n_out);
		free(names);
	}

	return dim;
error:
	free(names);
	isl_dim_free(dim);
	return NULL;
}

struct isl_dim *isl_dim_drop(struct isl_dim *dim, enum isl_dim_type type,
		unsigned first, unsigned num)
{
	int i;

	if (!dim)
		return NULL;

	if (n == 0)
		return dim;

	isl_assert(dim->ctx, first + num <= n(dim, type), goto error);
	dim = isl_dim_cow(dim);
	if (!dim)
		goto error;
	if (dim->names) {
		for (i = 0; i < num; ++i)
			isl_name_free(dim->ctx, get_name(dim, type, first+i));
		for (i = first+num; i < n(dim, type); ++i)
			set_name(dim, type, i - num, get_name(dim, type, i));
		switch (type) {
		case isl_dim_param:
			get_names(dim, isl_dim_in, 0, dim->n_in,
				dim->names + offset(dim, isl_dim_in) - num);
		case isl_dim_in:
			get_names(dim, isl_dim_out, 0, dim->n_out,
				dim->names + offset(dim, isl_dim_out) - num);
		case isl_dim_out:
			;
		}
	}
	switch (type) {
	case isl_dim_param:	dim->nparam -= num; break;
	case isl_dim_in:	dim->n_in -= num; break;
	case isl_dim_out:	dim->n_out -= num; break;
	}
	return dim;
error:
	isl_dim_free(dim);
	return NULL;
}

struct isl_dim *isl_dim_drop_inputs(struct isl_dim *dim,
		unsigned first, unsigned n)
{
	return isl_dim_drop(dim, isl_dim_in, first, n);
}

struct isl_dim *isl_dim_drop_outputs(struct isl_dim *dim,
		unsigned first, unsigned n)
{
	return isl_dim_drop(dim, isl_dim_out, first, n);
}

struct isl_dim *isl_dim_domain(struct isl_dim *dim)
{
	if (!dim)
		return NULL;
	dim = isl_dim_drop_outputs(dim, 0, dim->n_out);
	return isl_dim_reverse(dim);
}

struct isl_dim *isl_dim_range(struct isl_dim *dim)
{
	if (!dim)
		return NULL;
	return isl_dim_drop_inputs(dim, 0, dim->n_in);
}

struct isl_dim *isl_dim_underlying(struct isl_dim *dim, unsigned n_div)
{
	int i;

	if (!dim)
		return NULL;
	if (n_div == 0 &&
	    dim->nparam == 0 && dim->n_in == 0 && dim->n_name == 0)
		return dim;
	dim = isl_dim_cow(dim);
	if (!dim)
		return NULL;
	dim->n_out += dim->nparam + dim->n_in + n_div;
	dim->nparam = 0;
	dim->n_in = 0;

	for (i = 0; i < dim->n_name; ++i)
		isl_name_free(dim->ctx, get_name(dim, isl_dim_out, i));
	dim->n_name = 0;

	return dim;
}

unsigned isl_dim_total(struct isl_dim *dim)
{
	return dim->nparam + dim->n_in + dim->n_out;
}

int isl_dim_equal(struct isl_dim *dim1, struct isl_dim *dim2)
{
	return match(dim1, isl_dim_param, dim2, isl_dim_param) &&
	       match(dim1, isl_dim_in, dim2, isl_dim_in) &&
	       match(dim1, isl_dim_out, dim2, isl_dim_out);
}

int isl_dim_compatible(struct isl_dim *dim1, struct isl_dim *dim2)
{
	return dim1->nparam == dim2->nparam &&
	       dim1->n_in + dim1->n_out == dim2->n_in + dim2->n_out;
}
