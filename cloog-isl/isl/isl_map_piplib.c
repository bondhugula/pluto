#include "isl_set.h"
#include "isl_map.h"
#include "isl_mat.h"
#include "isl_seq.h"
#include "isl_piplib.h"
#include "isl_map_piplib.h"
#include "isl_map_private.h"
#include "isl_equalities.h"

static void copy_values_from(isl_int *dst, Entier *src, unsigned n)
{
	int i;

	for (i = 0; i < n; ++i)
		entier_assign(dst[i], src[i]);
}

static void add_value(isl_int *dst, Entier *src)
{
	mpz_add(*dst, *dst, *src);
}

static void copy_constraint_from(isl_int *dst, PipVector *src,
		unsigned nparam, unsigned n_in, unsigned n_out,
		unsigned extra, int *pos)
{
	int i;

	copy_values_from(dst, src->the_vector+src->nb_elements-1, 1);
	copy_values_from(dst+1, src->the_vector, nparam+n_in);
	isl_seq_clr(dst+1+nparam+n_in, n_out);
	isl_seq_clr(dst+1+nparam+n_in+n_out, extra);
	for (i = 0; i + n_in + nparam < src->nb_elements-1; ++i) {
		int p = pos[i];
		add_value(&dst[1+nparam+n_in+n_out+p],
			  &src->the_vector[n_in+nparam+i]);
	}
}

static int add_inequality(struct isl_ctx *ctx,
		   struct isl_basic_map *bmap, int *pos, PipVector *vec)
{
	unsigned nparam = isl_basic_map_n_param(bmap);
	unsigned n_in = isl_basic_map_n_in(bmap);
	unsigned n_out = isl_basic_map_n_out(bmap);
	unsigned n_div = isl_basic_map_n_div(bmap);
	int i = isl_basic_map_alloc_inequality(bmap);
	if (i < 0)
		return -1;
	copy_constraint_from(bmap->ineq[i], vec,
	    nparam, n_in, n_out, n_div, pos);

	return i;
}

/* For a div d = floor(f/m), add the constraints
 *
 *		f - m d >= 0
 *		-(f-(n-1)) + m d >= 0
 *
 * Note that the second constraint is the negation of
 *
 *		f - m d >= n
 */
static int add_div_constraints(struct isl_ctx *ctx,
	struct isl_basic_map *bmap, int *pos, PipNewparm *p, unsigned div)
{
	int i, j;
	unsigned total = isl_basic_map_total_dim(bmap);
	unsigned div_pos = 1 + total - bmap->n_div + div;

	i = add_inequality(ctx, bmap, pos, p->vector);
	if (i < 0)
		return -1;
	copy_values_from(&bmap->ineq[i][div_pos], &p->deno, 1);
	isl_int_neg(bmap->ineq[i][div_pos], bmap->ineq[i][div_pos]);

	j = isl_basic_map_alloc_inequality(bmap);
	if (j < 0)
		return -1;
	isl_seq_neg(bmap->ineq[j], bmap->ineq[i], 1 + total);
	isl_int_add(bmap->ineq[j][0], bmap->ineq[j][0], bmap->ineq[j][div_pos]);
	isl_int_sub_ui(bmap->ineq[j][0], bmap->ineq[j][0], 1);
	return j;
}

static int add_equality(struct isl_ctx *ctx,
		   struct isl_basic_map *bmap, int *pos,
		   unsigned var, PipVector *vec)
{
	int i;
	unsigned nparam = isl_basic_map_n_param(bmap);
	unsigned n_in = isl_basic_map_n_in(bmap);
	unsigned n_out = isl_basic_map_n_out(bmap);

	isl_assert(ctx, var < n_out, return -1);

	i = isl_basic_map_alloc_equality(bmap);
	if (i < 0)
		return -1;
	copy_constraint_from(bmap->eq[i], vec,
	    nparam, n_in, n_out, bmap->extra, pos);
	isl_int_set_si(bmap->eq[i][1+nparam+n_in+var], -1);

	return i;
}

static int find_div(struct isl_ctx *ctx,
		   struct isl_basic_map *bmap, int *pos, PipNewparm *p)
{
	int i, j;
	unsigned nparam = isl_basic_map_n_param(bmap);
	unsigned n_in = isl_basic_map_n_in(bmap);
	unsigned n_out = isl_basic_map_n_out(bmap);

	i = isl_basic_map_alloc_div(bmap);
	if (i < 0)
		return -1;

	copy_constraint_from(bmap->div[i]+1, p->vector,
	    nparam, n_in, n_out, bmap->extra, pos);

	copy_values_from(bmap->div[i], &p->deno, 1);
	for (j = 0; j < i; ++j)
		if (isl_seq_eq(bmap->div[i], bmap->div[j],
				1+1+isl_basic_map_total_dim(bmap)+j)) {
			isl_basic_map_free_div(bmap, 1);
			return j;
		}

	if (add_div_constraints(ctx, bmap, pos, p, i) < 0)
		return -1;

	return i;
}

/* Count some properties of a quast
 * - maximal number of new parameters
 * - maximal depth
 * - total number of solutions
 * - total number of empty branches
 */
static void quast_count(PipQuast *q, int *maxnew, int depth, int *maxdepth,
		        int *sol, int *nosol)
{
	PipNewparm *p;

	for (p = q->newparm; p; p = p->next)
		if (p->rank > *maxnew)
			*maxnew = p->rank;
	if (q->condition) {
		if (++depth > *maxdepth)
			*maxdepth = depth;
		quast_count(q->next_else, maxnew, depth, maxdepth, sol, nosol);
		quast_count(q->next_then, maxnew, depth, maxdepth, sol, nosol);
	} else {
		if (q->list)
			++(*sol);
		else
			++(*nosol);
	}
}

/*
 * pos: array of length bmap->set.extra, mapping each of the existential
 *		variables PIP proposes to an existential variable in bmap
 * bmap: collects the currently active constraints
 * rest: collects the empty leaves of the quast (if not NULL)
 */
struct scan_data {
	struct isl_ctx 			*ctx;
	struct isl_basic_map 		*bmap;
	struct isl_set			**rest;
	int	   *pos;
};

/*
 * New existentially quantified variables are places after the existing ones.
 */
static struct isl_map *scan_quast_r(struct scan_data *data, PipQuast *q,
				    struct isl_map *map)
{
	PipNewparm *p;
	struct isl_basic_map *bmap = data->bmap;
	unsigned old_n_div = bmap->n_div;
	unsigned nparam = isl_basic_map_n_param(bmap);
	unsigned n_in = isl_basic_map_n_in(bmap);
	unsigned n_out = isl_basic_map_n_out(bmap);

	if (!map)
		goto error;

	for (p = q->newparm; p; p = p->next) {
		int pos;
		unsigned pip_param = nparam + n_in;

		pos = find_div(data->ctx, bmap, data->pos, p);
		if (pos < 0)
			goto error;
		data->pos[p->rank - pip_param] = pos;
	}

	if (q->condition) {
		int pos = add_inequality(data->ctx, bmap, data->pos,
					 q->condition);
		if (pos < 0)
			goto error;
		map = scan_quast_r(data, q->next_then, map);

		if (isl_inequality_negate(bmap, pos))
			goto error;
		map = scan_quast_r(data, q->next_else, map);

		if (isl_basic_map_free_inequality(bmap, 1))
			goto error;
	} else if (q->list) {
		PipList *l;
		int j;
		/* if bmap->n_out is zero, we are only interested in the domains
		 * where a solution exists and not in the actual solution
		 */
		for (j = 0, l = q->list; j < n_out && l; ++j, l = l->next)
			if (add_equality(data->ctx, bmap, data->pos, j,
						l->vector) < 0)
				goto error;
		map = isl_map_add(map, isl_basic_map_copy(bmap));
		if (isl_basic_map_free_equality(bmap, n_out))
			goto error;
	} else if (data->rest) {
		struct isl_basic_set *bset;
		bset = isl_basic_set_from_basic_map(isl_basic_map_copy(bmap));
		bset = isl_basic_set_drop_dims(bset, n_in, n_out);
		if (!bset)
			goto error;
		*data->rest = isl_set_add(*data->rest, bset);
	}

	if (isl_basic_map_free_inequality(bmap, 2*(bmap->n_div - old_n_div)))
		goto error;
	if (isl_basic_map_free_div(bmap, bmap->n_div - old_n_div))
		goto error;
	return map;
error:
	isl_map_free(map);
	return NULL;
}

/*
 * Returns a map of dimension "keep_dim" with "context" as domain and
 * as range the first "isl_dim_size(keep_dim, isl_dim_out)" variables
 * in the quast lists.
 */
static struct isl_map *isl_map_from_quast(struct isl_ctx *ctx, PipQuast *q,
		struct isl_dim *keep_dim,
		struct isl_basic_set *context,
		struct isl_set **rest)
{
	int		pip_param;
	int		nexist;
	int		max_depth;
	int		n_sol, n_nosol;
	struct scan_data	data;
	struct isl_map		*map = NULL;
	struct isl_dim		*dims;
	unsigned		nparam;
	unsigned		dim;
	unsigned		keep;

	data.ctx = ctx;
	data.rest = rest;
	data.bmap = NULL;
	data.pos = NULL;

	if (!context || !keep_dim)
		goto error;

	dim = isl_basic_set_n_dim(context);
	nparam = isl_basic_set_n_param(context);
	keep = isl_dim_size(keep_dim, isl_dim_out);
	pip_param = nparam + dim;

	max_depth = 0;
	n_sol = 0;
	n_nosol = 0;
	nexist = pip_param-1;
	quast_count(q, &nexist, 0, &max_depth, &n_sol, &n_nosol);
	nexist -= pip_param-1;

	if (rest) {
		*rest = isl_set_alloc_dim(isl_dim_copy(context->dim), n_nosol,
					ISL_MAP_DISJOINT);
		if (!*rest)
			goto error;
	}
	map = isl_map_alloc_dim(isl_dim_copy(keep_dim), n_sol,
				ISL_MAP_DISJOINT);
	if (!map)
		goto error;

	dims = isl_dim_reverse(isl_dim_copy(context->dim));
	data.bmap = isl_basic_map_from_basic_set(context, dims);
	data.bmap = isl_basic_map_extend_dim(data.bmap,
		keep_dim, nexist, keep, max_depth+2*nexist);
	if (!data.bmap)
		goto error2;

	if (data.bmap->extra) {
		int i;
		data.pos = isl_alloc_array(ctx, int, data.bmap->extra);
		if (!data.pos)
			goto error;
		for (i = 0; i < data.bmap->n_div; ++i)
			data.pos[i] = i;
	}

	map = scan_quast_r(&data, q, map);
	map = isl_map_finalize(map);
	if (!map)
		goto error2;
	if (rest) {
		*rest = isl_set_finalize(*rest);
		if (!*rest)
			goto error2;
	}
	isl_basic_map_free(data.bmap);
	if (data.pos)
		free(data.pos);
	return map;
error:
	isl_basic_set_free(context);
	isl_dim_free(keep_dim);
error2:
	if (data.pos)
		free(data.pos);
	isl_basic_map_free(data.bmap);
	isl_map_free(map);
	if (rest) {
		isl_set_free(*rest);
		*rest = NULL;
	}
	return NULL;
}

static void copy_values_to(Entier *dst, isl_int *src, unsigned n)
{
	int i;

	for (i = 0; i < n; ++i)
		entier_assign(dst[i], src[i]);
}

static void copy_constraint_to(Entier *dst, isl_int *src,
		unsigned pip_param, unsigned pip_var,
		unsigned extra_front, unsigned extra_back)
{
	copy_values_to(dst+1+extra_front+pip_var+pip_param+extra_back, src, 1);
	copy_values_to(dst+1+extra_front+pip_var, src+1, pip_param);
	copy_values_to(dst+1+extra_front, src+1+pip_param, pip_var);
}

PipMatrix *isl_basic_map_to_pip(struct isl_basic_map *bmap, unsigned pip_param,
			 unsigned extra_front, unsigned extra_back)
{
	int i;
	unsigned nrow;
	unsigned ncol;
	PipMatrix *M;
	unsigned off;
	unsigned pip_var = isl_basic_map_total_dim(bmap) - pip_param;

	nrow = extra_front + bmap->n_eq + bmap->n_ineq;
	ncol = 1 + extra_front + pip_var + pip_param + extra_back + 1;
	M = pip_matrix_alloc(nrow, ncol);
	if (!M)
		return NULL;

	off = extra_front;
	for (i = 0; i < bmap->n_eq; ++i) {
		entier_set_si(M->p[off+i][0], 0);
		copy_constraint_to(M->p[off+i], bmap->eq[i],
				   pip_param, pip_var, extra_front, extra_back);
	}
	off += bmap->n_eq;
	for (i = 0; i < bmap->n_ineq; ++i) {
		entier_set_si(M->p[off+i][0], 1);
		copy_constraint_to(M->p[off+i], bmap->ineq[i],
				   pip_param, pip_var, extra_front, extra_back);
	}
	return M;
}

PipMatrix *isl_basic_set_to_pip(struct isl_basic_set *bset, unsigned pip_param,
			 unsigned extra_front, unsigned extra_back)
{
	return isl_basic_map_to_pip((struct isl_basic_map *)bset,
					pip_param, extra_front, extra_back);
}

static struct isl_map *extremum_on(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty, int max)
{
	PipOptions	*options;
	PipQuast	*sol;
	struct isl_map	*map;
	struct isl_ctx  *ctx;
	PipMatrix *domain = NULL, *context = NULL;
	unsigned	 nparam, n_in, n_out;

	bmap = isl_basic_map_detect_equalities(bmap);
	if (!bmap || !dom)
		goto error;

	ctx = bmap->ctx;
	isl_assert(ctx, isl_basic_map_compatible_domain(bmap, dom), goto error);
	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);

	domain = isl_basic_map_to_pip(bmap, nparam + n_in, 0, dom->n_div);
	if (!domain)
		goto error;
	context = isl_basic_map_to_pip((struct isl_basic_map *)dom, 0, 0, 0);
	if (!context)
		goto error;

	options = pip_options_init();
	options->Simplify = 1;
	options->Maximize = max;
	options->Urs_unknowns = -1;
	options->Urs_parms = -1;
	sol = pip_solve(domain, context, -1, options);

	if (sol) {
		struct isl_basic_set *copy;
		copy = isl_basic_set_copy(dom);
		map = isl_map_from_quast(ctx, sol,
				isl_dim_copy(bmap->dim), copy, empty);
	} else {
		map = isl_map_empty_like_basic_map(bmap);
		if (empty)
			*empty = NULL;
	}
	if (!map)
		goto error;
	if (map->n == 0 && empty) {
		isl_set_free(*empty);
		*empty = isl_set_from_basic_set(dom);
	} else
		isl_basic_set_free(dom);
	isl_basic_map_free(bmap);

	pip_quast_free(sol);
	pip_options_free(options);
	pip_matrix_free(domain);
	pip_matrix_free(context);

	return map;
error:
	if (domain)
		pip_matrix_free(domain);
	if (context)
		pip_matrix_free(context);
	isl_basic_map_free(bmap);
	isl_basic_set_free(dom);
	return NULL;
}

struct isl_map *isl_pip_basic_map_lexmax(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty)
{
	return extremum_on(bmap, dom, empty, 1);
}

struct isl_map *isl_pip_basic_map_lexmin(
		struct isl_basic_map *bmap, struct isl_basic_set *dom,
		struct isl_set **empty)
{
	return extremum_on(bmap, dom, empty, 0);
}

/* Project the given basic set onto its parameter domain, possibly introducing
 * new, explicit, existential variables in the constraints.
 * The input has parameters and output variables.
 * The output has the same parameters, but no output variables, only
 * explicit existentially quantified variables.
 */
static struct isl_set *compute_divs_no_eq(struct isl_basic_set *bset)
{
	PipMatrix *domain = NULL, *context = NULL;
	PipOptions	*options;
	PipQuast	*sol;
	struct isl_ctx  *ctx;
	struct isl_dim	*dim;
	struct isl_map	*map;
	struct isl_set	*set;
	struct isl_basic_set	*dom;
	unsigned	 nparam;

	if (!bset)
		goto error;

	ctx = bset->ctx;
	nparam = isl_basic_set_n_param(bset);

	domain = isl_basic_set_to_pip(bset, nparam, 0, 0);
	if (!domain)
		goto error;
	context = pip_matrix_alloc(0, nparam + 2);
	if (!context)
		goto error;

	options = pip_options_init();
	options->Simplify = 1;
	options->Urs_unknowns = -1;
	options->Urs_parms = -1;
	sol = pip_solve(domain, context, -1, options);

	dom = isl_basic_set_alloc(ctx, nparam, 0, 0, 0, 0);
	map = isl_map_from_quast(ctx, sol, isl_dim_copy(dom->dim), dom, NULL);
	set = (struct isl_set *)map;

	pip_quast_free(sol);
	pip_options_free(options);
	pip_matrix_free(domain);
	pip_matrix_free(context);

	isl_basic_set_free(bset);

	return set;
error:
	if (domain)
		pip_matrix_free(domain);
	if (context)
		pip_matrix_free(context);
	isl_basic_set_free(dom);
	isl_basic_set_free(bset);
	return NULL;
}

static struct isl_map *isl_map_reset_dim(struct isl_map *map,
	struct isl_dim *dim)
{
	int i;

	if (!map || !dim)
		goto error;

	for (i = 0; i < map->n; ++i) {
		isl_dim_free(map->p[i]->dim);
		map->p[i]->dim = isl_dim_copy(dim);
	}
	isl_dim_free(map->dim);
	map->dim = dim;

	return map;
error:
	isl_map_free(map);
	isl_dim_free(dim);
	return NULL;
}

static struct isl_set *isl_set_reset_dim(struct isl_set *set,
	struct isl_dim *dim)
{
	return (struct isl_set *) isl_map_reset_dim((struct isl_map *)set, dim);
}

/* Given a matrix M (mat) and a size n (size), replace mat
 * by the matrix
 *
 *		[ M 0 ]
 *		[ 0 I ]
 *
 * where I is an n x n identity matrix.
 */
static struct isl_mat *append_identity(struct isl_mat *mat, unsigned size)
{
	int i;
	unsigned n_row, n_col;

	n_row = mat->n_row;
	n_col = mat->n_col;
	mat = isl_mat_extend(mat, n_row + size, n_col + size);
	if (!mat)
		return NULL;
	for (i = 0; i < n_row; ++i)
		isl_seq_clr(mat->row[i] + n_col, size);
	for (i = 0; i < size; ++i) {
		isl_seq_clr(mat->row[n_row + i], n_col + size);
		isl_int_set_si(mat->row[n_row + i][n_col + i], 1);
	}
	return mat;
}

/* Apply a preimage specified by "mat" on the parameters of "bset".
 */
static struct isl_basic_set *basic_set_parameter_preimage(
	struct isl_basic_set *bset, struct isl_mat *mat)
{
	unsigned nparam, n_out;

	if (!bset || !mat)
		goto error;

	bset->dim = isl_dim_cow(bset->dim);
	if (!bset->dim)
		goto error;

	nparam = isl_basic_set_dim(bset, isl_dim_param);
	n_out = isl_basic_set_dim(bset, isl_dim_set);

	isl_assert(bset->ctx, mat->n_row == 1 + nparam, goto error);

	mat = append_identity(mat, n_out);
	if (!mat)
		goto error;

	bset->dim->nparam = 0;
	bset->dim->n_out += nparam;
	bset = isl_basic_set_preimage(bset, mat);
	if (bset) {
		bset->dim->nparam = bset->dim->n_out - n_out;
		bset->dim->n_out = n_out;
	}
	return bset;
error:
	isl_mat_free(mat);
	isl_basic_set_free(bset);
	return NULL;
}

/* Apply a preimage specified by "mat" on the parameters of "set".
 */
static struct isl_set *set_parameter_preimage(
	struct isl_set *set, struct isl_mat *mat)
{
	struct isl_dim *dim = NULL;
	unsigned nparam, n_out;

	if (!set || !mat)
		goto error;

	dim = isl_dim_copy(set->dim);
	dim = isl_dim_cow(dim);
	if (!dim)
		goto error;

	nparam = isl_set_dim(set, isl_dim_param);
	n_out = isl_set_dim(set, isl_dim_set);

	isl_assert(set->ctx, mat->n_row == 1 + nparam, goto error);

	mat = append_identity(mat, n_out);
	if (!mat)
		goto error;

	dim->nparam = 0;
	dim->n_out += nparam;
	isl_set_reset_dim(set, dim);
	set = isl_set_preimage(set, mat);
	if (!set)
		goto error2;
	dim = isl_dim_copy(set->dim);
	dim = isl_dim_cow(dim);
	if (!dim)
		goto error2;
	dim->nparam = dim->n_out - n_out;
	dim->n_out = n_out;
	isl_set_reset_dim(set, dim);
	return set;
error:
	isl_dim_free(dim);
	isl_mat_free(mat);
error2:
	isl_set_free(set);
	return NULL;
}

/* Intersect the basic set "bset" with the affine space specified by the
 * equalities in "eq".
 */
static struct isl_basic_set *basic_set_append_equalities(
	struct isl_basic_set *bset, struct isl_mat *eq)
{
	int i, k;
	unsigned len;

	if (!bset || !eq)
		goto error;

	bset = isl_basic_set_extend_dim(bset, isl_dim_copy(bset->dim), 0,
					eq->n_row, 0);
	if (!bset)
		goto error;

	len = 1 + isl_dim_total(bset->dim) + bset->extra;
	for (i = 0; i < eq->n_row; ++i) {
		k = isl_basic_set_alloc_equality(bset);
		if (k < 0)
			goto error;
		isl_seq_cpy(bset->eq[k], eq->row[i], eq->n_col);
		isl_seq_clr(bset->eq[k] + eq->n_col, len - eq->n_col);
	}
	isl_mat_free(eq);

	return bset;
error:
	isl_mat_free(eq);
	isl_basic_set_free(bset);
	return NULL;
}

/* Intersect the set "set" with the affine space specified by the
 * equalities in "eq".
 */
static struct isl_set *set_append_equalities(struct isl_set *set,
	struct isl_mat *eq)
{
	int i;

	if (!set || !eq)
		goto error;

	for (i = 0; i < set->n; ++i) {
		set->p[i] = basic_set_append_equalities(set->p[i],
					isl_mat_copy(eq));
		if (!set->p[i])
			goto error;
	}
	isl_mat_free(eq);
	return set;
error:
	isl_mat_free(eq);
	isl_set_free(set);
	return NULL;
}

/* Project the given basic set onto its parameter domain, possibly introducing
 * new, explicit, existential variables in the constraints.
 * The input has parameters and output variables.
 * The output has the same parameters, but no output variables, only
 * explicit existentially quantified variables.
 *
 * The actual projection is performed by pip, but pip doesn't seem
 * to like equalities very much, so we first remove the equalities
 * among the parameters by performing a variable compression on
 * the parameters.  Afterward, an inverse transformation is performed
 * and the equalities among the parameters are inserted back in.
 */
static struct isl_set *compute_divs(struct isl_basic_set *bset)
{
	int i, j;
	struct isl_mat *eq;
	struct isl_mat *T, *T2;
	struct isl_set *set;
	unsigned nparam, n_out;

	bset = isl_basic_set_cow(bset);
	if (!bset)
		return NULL;

	if (bset->n_eq == 0)
		return compute_divs_no_eq(bset);

	isl_basic_set_gauss(bset, NULL);

	nparam = isl_basic_set_dim(bset, isl_dim_param);
	n_out = isl_basic_set_dim(bset, isl_dim_out);

	for (i = 0, j = n_out - 1; i < bset->n_eq && j >= 0; --j) {
		if (!isl_int_is_zero(bset->eq[i][1 + nparam + j]))
			++i;
	}
	if (i == bset->n_eq)
		return compute_divs_no_eq(bset);

	eq = isl_mat_sub_alloc(bset->ctx, bset->eq, i, bset->n_eq - i,
		0, 1 + nparam);
	eq = isl_mat_cow(eq);
	T = isl_mat_variable_compression(isl_mat_copy(eq), &T2);
	bset = basic_set_parameter_preimage(bset, T);

	set = compute_divs_no_eq(bset);
	set = set_parameter_preimage(set, T2);
	set = set_append_equalities(set, eq);
	return set;
}

/* Compute an explicit representation for all the existentially
 * quantified variables.
 * The input and output dimensions are first turned into parameters
 * and the existential variables into output dimensions.
 * compute_divs then returns a map with the same parameters and
 * no input or output dimensions and the dimension specification
 * is reset to that of the input.
 */
struct isl_map *isl_pip_basic_map_compute_divs(struct isl_basic_map *bmap)
{
	struct isl_basic_set *bset;
	struct isl_set *set;
	struct isl_map *map;
	struct isl_dim *dim, *orig_dim = NULL;
	unsigned	 nparam;
	unsigned	 n_in;
	unsigned	 n_out;

	bmap = isl_basic_map_cow(bmap);
	if (!bmap)
		return NULL;

	nparam = isl_basic_map_dim(bmap, isl_dim_param);
	n_in = isl_basic_map_dim(bmap, isl_dim_in);
	n_out = isl_basic_map_dim(bmap, isl_dim_out);
	dim = isl_dim_set_alloc(bmap->ctx, nparam + n_in + n_out, bmap->n_div);
	if (!dim)
		goto error;

	orig_dim = bmap->dim;
	bmap->dim = dim;
	bmap->extra -= bmap->n_div;
	bmap->n_div = 0;
	bset = (struct isl_basic_set *)bmap;

	set = compute_divs(bset);
	map = (struct isl_map *)set;
	map = isl_map_reset_dim(map, orig_dim);

	return map;
error:
	isl_basic_map_free(bmap);
	return NULL;
}
