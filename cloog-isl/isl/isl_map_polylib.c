#include "isl_set.h"
#include "isl_map.h"
#include "isl_seq.h"
#include "isl_set_polylib.h"
#include "isl_map_polylib.h"
#include "isl_map_private.h"

static void copy_values_from(isl_int *dst, Value *src, unsigned n)
{
	int i;

	for (i = 0; i < n; ++i)
		value_assign(dst[i], src[i]);
}

static void copy_values_to(Value *dst, isl_int *src, unsigned n)
{
	int i;

	for (i = 0; i < n; ++i)
		value_assign(dst[i], src[i]);
}

static void copy_constraint_from(isl_int *dst, Value *src,
		unsigned nparam, unsigned dim, unsigned extra)
{
	copy_values_from(dst, src+1+dim+extra+nparam, 1);
	copy_values_from(dst+1, src+1+dim+extra, nparam);
	copy_values_from(dst+1+nparam, src+1, dim);
	copy_values_from(dst+1+nparam+dim, src+1+dim, extra);
}

static void copy_constraint_to(Value *dst, isl_int *src,
		unsigned nparam, unsigned dim, unsigned extra)
{
	copy_values_to(dst+1+dim+extra+nparam, src, 1);
	copy_values_to(dst+1+dim+extra, src+1, nparam);
	copy_values_to(dst+1, src+1+nparam, dim);
	copy_values_to(dst+1+dim, src+1+nparam+dim, extra);
}

static int add_equality(struct isl_ctx *ctx, struct isl_basic_map *bmap,
			 Value *constraint)
{
	unsigned nparam;
	unsigned n_in;
	unsigned n_out;
	int i = isl_basic_map_alloc_equality(bmap);
	if (i < 0)
		return -1;
	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);
	copy_constraint_from(bmap->eq[i], constraint, nparam,
			     n_in + n_out, bmap->extra);
	return 0;
}

static int add_inequality(struct isl_ctx *ctx, struct isl_basic_map *bmap,
			 Value *constraint)
{
	unsigned nparam;
	unsigned n_in;
	unsigned n_out;
	int i = isl_basic_map_alloc_inequality(bmap);
	if (i < 0)
		return -1;
	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);
	copy_constraint_from(bmap->ineq[i], constraint, nparam,
			     n_in + n_out, bmap->extra);
	return 0;
}

static struct isl_basic_map *copy_constraints(
			struct isl_ctx *ctx, struct isl_basic_map *bmap,
			Polyhedron *P)
{
	int i;
	unsigned total = isl_basic_map_total_dim(bmap);

	for (i = 0; i < P->NbConstraints; ++i) {
		if (value_zero_p(P->Constraint[i][0])) {
			if (add_equality(ctx, bmap, P->Constraint[i]))
				goto error;
		} else {
			if (add_inequality(ctx, bmap, P->Constraint[i]))
				goto error;
		}
	}
	for (i = 0; i < bmap->extra; ++i) {
		int j = isl_basic_map_alloc_div(bmap);
		if (j == -1)
			goto error;
		isl_seq_clr(bmap->div[j], 1+1+total);
	}
	return bmap;
error:
	isl_basic_map_free(bmap);
	return NULL;
}

struct isl_basic_set *isl_basic_set_new_from_polylib(Polyhedron *P,
			struct isl_dim *dim)
{
	if (!dim)
		return NULL;
	isl_assert(dim->ctx, dim->n_in == 0, return NULL);

	return (struct isl_basic_set *)
		isl_basic_map_new_from_polylib(P, dim);
}

struct isl_basic_map *isl_basic_map_new_from_polylib(Polyhedron *P,
			struct isl_dim *dim)
{
	struct isl_basic_map *bmap;
	unsigned extra;

	if (!dim)
		return NULL;

	isl_assert(dim->ctx, P, goto error);
	isl_assert(dim->ctx, P->Dimension >= isl_dim_total(dim), goto error);

	extra = P->Dimension - isl_dim_total(dim);
	bmap = isl_basic_map_alloc_dim(dim, extra,
					P->NbEq, P->NbConstraints - P->NbEq);
	if (!bmap)
		return NULL;

	bmap = copy_constraints(dim->ctx, bmap, P);
	bmap = isl_basic_map_simplify(bmap);
	return isl_basic_map_finalize(bmap);
error:
	isl_dim_free(dim);
	return NULL;
}

struct isl_set *isl_set_new_from_polylib(Polyhedron *D, struct isl_dim *dim)
{
	struct isl_set *set = NULL;
	Polyhedron *P;
	int n = 0;

	if (!dim)
		return NULL;
	isl_assert(dim->ctx, dim->n_in == 0, goto error);

	for (P = D; P; P = P->next)
		++n;

	set = isl_set_alloc_dim(isl_dim_copy(dim), n, ISL_MAP_DISJOINT);
	if (!set)
		goto error;

	for (P = D; P; P = P->next)
		isl_set_add(set,
		    isl_basic_set_new_from_polylib(P, isl_dim_copy(dim)));
	isl_dim_free(dim);
	set = isl_set_remove_empty_parts(set);
	return set;
error:
	isl_dim_free(dim);
	return NULL;
}

struct isl_map *isl_map_new_from_polylib(Polyhedron *D, struct isl_dim *dim)
{
	struct isl_map *map = NULL;
	Polyhedron *P;
	int n = 0;

	if (!dim)
		return NULL;

	for (P = D; P; P = P->next)
		++n;

	map = isl_map_alloc_dim(isl_dim_copy(dim), n, ISL_MAP_DISJOINT);
	if (!map)
		goto error;

	for (P = D; P; P = P->next)
		isl_map_add(map,
		    isl_basic_map_new_from_polylib(P, isl_dim_copy(dim)));
	isl_dim_free(dim);
	map = isl_map_remove_empty_parts(map);
	return map;
error:
	isl_dim_free(dim);
	return NULL;
}

Polyhedron *isl_basic_map_to_polylib(struct isl_basic_map *bmap)
{
	int i;
	Matrix *M;
	Polyhedron *P;
	unsigned off;
	unsigned nparam;
	unsigned n_in;
	unsigned n_out;

	if (!bmap)
		return NULL;

	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);
	M = Matrix_Alloc(bmap->n_eq + bmap->n_ineq,
		 1 + n_in + n_out + bmap->n_div + nparam + 1);
	for (i = 0; i < bmap->n_eq; ++i) {
		value_set_si(M->p[i][0], 0);
		copy_constraint_to(M->p[i], bmap->eq[i],
			   nparam, n_in + n_out, bmap->n_div);
	}
	off = bmap->n_eq;
	for (i = 0; i < bmap->n_ineq; ++i) {
		value_set_si(M->p[off+i][0], 1);
		copy_constraint_to(M->p[off+i], bmap->ineq[i],
			   nparam, n_in + n_out, bmap->n_div);
	}
	P = Constraints2Polyhedron(M, bmap->ctx->MaxRays);
	Matrix_Free(M);

	return P;
}

Polyhedron *isl_map_to_polylib(struct isl_map *map)
{
	int i;
	Polyhedron *R = NULL;
	Polyhedron **next = &R;

	if (!map)
		return NULL;

	for (i = 0; i < map->n; ++i) {
		*next = isl_basic_map_to_polylib(map->p[i]);
		next = &(*next)->next;
	}

	return R ? R : Empty_Polyhedron(isl_dim_total(map->dim));
}

Polyhedron *isl_basic_set_to_polylib(struct isl_basic_set *bset)
{
	return isl_basic_map_to_polylib((struct isl_basic_map *)bset);
}

Polyhedron *isl_set_to_polylib(struct isl_set *set)
{
	return isl_map_to_polylib((struct isl_map *)set);
}
