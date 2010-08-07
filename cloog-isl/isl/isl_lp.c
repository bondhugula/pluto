#include "isl_ctx.h"
#include "isl_lp.h"
#include "isl_lp_piplib.h"
#include "isl_tab.h"
#include "isl_map_private.h"

enum isl_lp_result isl_tab_solve_lp(struct isl_basic_map *bmap, int maximize,
				      isl_int *f, isl_int denom, isl_int *opt,
				      isl_int *opt_denom)
{
	struct isl_tab *tab;
	enum isl_lp_result res;
	unsigned dim = isl_basic_map_total_dim(bmap);

	if (maximize)
		isl_seq_neg(f, f, 1 + dim);

	bmap = isl_basic_map_gauss(bmap, NULL);
	tab = isl_tab_from_basic_map(bmap);
	res = isl_tab_min(tab, f, bmap->ctx->one, opt, opt_denom, 0);
	isl_tab_free(tab);

	if (maximize)
		isl_seq_neg(f, f, 1 + dim);

	return res;
}

/* Given a basic map "bmap" and an affine combination of the variables "f"
 * with denominator "denom", set *opt/*opt_denom to the minimal
 * (or maximal if "maximize" is true) value attained by f/d over "bmap",
 * assuming the basic map is not empty and the expression cannot attain
 * arbitrarily small (or large) values.
 * If opt_denom is NULL, then *opt is rounded up (or down)
 * to the nearest integer.
 * The return value reflects the nature of the result (empty, unbounded,
 * minmimal or maximal value returned in *opt).
 */
enum isl_lp_result isl_solve_lp(struct isl_basic_map *bmap, int maximize,
				      isl_int *f, isl_int d, isl_int *opt,
				      isl_int *opt_denom)
{
	if (!bmap)
		return isl_lp_error;

	switch (bmap->ctx->lp_solver) {
	case ISL_LP_PIP:
		return isl_pip_solve_lp(bmap, maximize, f, d, opt, opt_denom);
	case ISL_LP_TAB:
		return isl_tab_solve_lp(bmap, maximize, f, d, opt, opt_denom);
	default:
		return isl_lp_error;
	}
}
