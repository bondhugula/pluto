#include <assert.h>
#include "isl_tab.h"

struct tab_lp {
	struct isl_ctx  *ctx;
	struct isl_vec  *row;
	struct isl_tab  *tab;
	struct isl_tab_undo	**stack;
	isl_int		*obj;
	isl_int		 opt;
	isl_int		 opt_denom;
	int	         neq;
	unsigned	 dim;
	int		 n_ineq;
};

static struct tab_lp *init_lp(struct isl_basic_set *bset);
static void set_lp_obj(struct tab_lp *lp, isl_int *row, int dim);
static int solve_lp(struct tab_lp *lp);
static void get_obj_val(struct tab_lp* lp, mpq_t *F);
static void delete_lp(struct tab_lp *lp);
static int add_lp_row(struct tab_lp *lp, isl_int *row, int dim);
static void get_alpha(struct tab_lp* lp, int row, mpq_t *alpha);
static void del_lp_row(struct tab_lp *lp);

#define GBR_LP			    	    struct tab_lp
#define GBR_type		    	    mpq_t
#define GBR_init(v)		    	    mpq_init(v)
#define GBR_clear(v)		    	    mpq_clear(v)
#define GBR_set(a,b)			    mpq_set(a,b)
#define GBR_set_ui(a,b)			    mpq_set_ui(a,b,1)
#define GBR_mul(a,b,c)			    mpq_mul(a,b,c)
#define GBR_lt(a,b)			    (mpq_cmp(a,b) < 0)
#define GBR_floor(a,b)			    mpz_fdiv_q(a,mpq_numref(b),mpq_denref(b))
#define GBR_ceil(a,b)			    mpz_cdiv_q(a,mpq_numref(b),mpq_denref(b))
#define GBR_lp_init(P)		    	    init_lp(P)
#define GBR_lp_set_obj(lp, obj, dim)	    set_lp_obj(lp, obj, dim)
#define GBR_lp_solve(lp)		    solve_lp(lp)
#define GBR_lp_get_obj_val(lp, F)	    get_obj_val(lp, F)
#define GBR_lp_delete(lp)		    delete_lp(lp)
#define GBR_lp_next_row(lp)		    lp->neq
#define GBR_lp_add_row(lp, row, dim)	    add_lp_row(lp, row, dim)
#define GBR_lp_get_alpha(lp, row, alpha)    get_alpha(lp, row, alpha)
#define GBR_lp_del_row(lp)		    del_lp_row(lp);
#include "basis_reduction_templ.c"

/* Set up a tableau for the Cartesian product of bset with itself.
 * This could be optimized by first setting up a tableau for bset
 * and then performing the Cartesian product on the tableau.
 */
static struct isl_tab *gbr_tab(struct isl_basic_set *bset,
	struct isl_vec *row)
{
	int i, j;
	unsigned dim;
	struct isl_tab *tab;

	if (!bset || !row)
		return NULL;

	dim = isl_basic_set_total_dim(bset);
	tab = isl_tab_alloc(bset->ctx, 2 * bset->n_ineq + dim + 1, 2 * dim);

	for (i = 0; i < 2; ++i) {
		isl_seq_clr(row->el + 1 + (1 - i) * dim, dim);
		for (j = 0; j < bset->n_ineq; ++j) {
			isl_int_set(row->el[0], bset->ineq[j][0]);
			isl_seq_cpy(row->el + 1 + i * dim,
				    bset->ineq[j] + 1, dim);
			tab = isl_tab_add_ineq(tab, row->el);
			if (!tab || tab->empty)
				return tab;
		}
	}

	return tab;
}

static struct tab_lp *init_lp(struct isl_basic_set *bset)
{
	struct tab_lp *lp = NULL;

	if (!bset)
		return NULL;

	isl_assert(bset->ctx, bset->n_eq == 0, return NULL);

	lp = isl_calloc_type(bset->ctx, struct tab_lp);
	if (!lp)
		return NULL;

	isl_int_init(lp->opt);
	isl_int_init(lp->opt_denom);

	lp->dim = isl_basic_set_total_dim(bset);
	lp->n_ineq = bset->n_ineq;

	lp->ctx = bset->ctx;
	isl_ctx_ref(lp->ctx);

	lp->stack = isl_alloc_array(lp->ctx, struct isl_tab_undo *, lp->dim);

	lp->row = isl_vec_alloc(lp->ctx, 1 + 2 * lp->dim);
	if (!lp->row)
		goto error;
	lp->tab = gbr_tab(bset, lp->row);
	if (!lp->tab)
		goto error;
	lp->obj = NULL;
	lp->neq = 0;

	return lp;
error:
	delete_lp(lp);
	return NULL;
}

static void set_lp_obj(struct tab_lp *lp, isl_int *row, int dim)
{
	lp->obj = row;
}

static int solve_lp(struct tab_lp *lp)
{
	enum isl_lp_result res;
	unsigned flags = 0;

	isl_int_set_si(lp->row->el[0], 0);
	isl_seq_cpy(lp->row->el + 1, lp->obj, lp->dim);
	isl_seq_neg(lp->row->el + 1 + lp->dim, lp->obj, lp->dim);
	if (lp->neq)
		flags = ISL_TAB_SAVE_DUAL;
	res = isl_tab_min(lp->tab, lp->row->el, lp->ctx->one,
			  &lp->opt, &lp->opt_denom, flags);
	if (res != isl_lp_ok)
		return -1;
	return 0;
}

static void get_obj_val(struct tab_lp* lp, mpq_t *F)
{
	isl_int_neg(mpq_numref(*F), lp->opt);
	isl_int_set(mpq_denref(*F), lp->opt_denom);
}

static void delete_lp(struct tab_lp *lp)
{
	if (!lp)
		return;

	isl_int_clear(lp->opt);
	isl_int_clear(lp->opt_denom);
	isl_vec_free(lp->row);
	free(lp->stack);
	isl_tab_free(lp->tab);
	isl_ctx_deref(lp->ctx);
	free(lp);
}

static int add_lp_row(struct tab_lp *lp, isl_int *row, int dim)
{
	lp->stack[lp->neq] = isl_tab_snap(lp->tab);

	isl_int_set_si(lp->row->el[0], 0);
	isl_seq_cpy(lp->row->el + 1, row, lp->dim);
	isl_seq_neg(lp->row->el + 1 + lp->dim, row, lp->dim);

	lp->tab = isl_tab_add_valid_eq(lp->tab, lp->row->el);

	return lp->neq++;
}

static void get_alpha(struct tab_lp* lp, int row, mpq_t *alpha)
{
	row += 2 * lp->n_ineq;
	isl_int_neg(mpq_numref(*alpha), lp->tab->dual->el[1 + row]);
	isl_int_set(mpq_denref(*alpha), lp->tab->dual->el[0]);
}

static void del_lp_row(struct tab_lp *lp)
{
	lp->neq--;
	isl_tab_rollback(lp->tab, lp->stack[lp->neq]);
}
