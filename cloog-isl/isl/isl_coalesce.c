#include "isl_map_private.h"
#include "isl_tab.h"

#define STATUS_ERROR		-1
#define STATUS_REDUNDANT	 1
#define STATUS_VALID	 	 2
#define STATUS_SEPARATE	 	 3
#define STATUS_CUT	 	 4
#define STATUS_ADJ_EQ	 	 5
#define STATUS_ADJ_INEQ	 	 6

static int status_in(isl_int *ineq, struct isl_tab *tab)
{
	enum isl_ineq_type type = isl_tab_ineq_type(tab, ineq);
	switch (type) {
	case isl_ineq_error:		return STATUS_ERROR;
	case isl_ineq_redundant:	return STATUS_VALID;
	case isl_ineq_separate:		return STATUS_SEPARATE;
	case isl_ineq_cut:		return STATUS_CUT;
	case isl_ineq_adj_eq:		return STATUS_ADJ_EQ;
	case isl_ineq_adj_ineq:		return STATUS_ADJ_INEQ;
	}
}

/* Compute the position of the equalities of basic map "i"
 * with respect to basic map "j".
 * The resulting array has twice as many entries as the number
 * of equalities corresponding to the two inequalties to which
 * each equality corresponds.
 */
static int *eq_status_in(struct isl_map *map, int i, int j,
	struct isl_tab **tabs)
{
	int k, l;
	int *eq = isl_calloc_array(map->ctx, int, 2 * map->p[i]->n_eq);
	unsigned dim;

	dim = isl_basic_map_total_dim(map->p[i]);
	for (k = 0; k < map->p[i]->n_eq; ++k) {
		for (l = 0; l < 2; ++l) {
			isl_seq_neg(map->p[i]->eq[k], map->p[i]->eq[k], 1+dim);
			eq[2 * k + l] = status_in(map->p[i]->eq[k], tabs[j]);
			if (eq[2 * k + l] == STATUS_ERROR)
				goto error;
		}
		if (eq[2 * k] == STATUS_SEPARATE ||
		    eq[2 * k + 1] == STATUS_SEPARATE)
			break;
	}

	return eq;
error:
	free(eq);
	return NULL;
}

/* Compute the position of the inequalities of basic map "i"
 * with respect to basic map "j".
 */
static int *ineq_status_in(struct isl_map *map, int i, int j,
	struct isl_tab **tabs)
{
	int k;
	unsigned n_eq = map->p[i]->n_eq;
	int *ineq = isl_calloc_array(map->ctx, int, map->p[i]->n_ineq);

	for (k = 0; k < map->p[i]->n_ineq; ++k) {
		if (isl_tab_is_redundant(tabs[i], n_eq + k)) {
			ineq[k] = STATUS_REDUNDANT;
			continue;
		}
		ineq[k] = status_in(map->p[i]->ineq[k], tabs[j]);
		if (ineq[k] == STATUS_ERROR)
			goto error;
		if (ineq[k] == STATUS_SEPARATE)
			break;
	}

	return ineq;
error:
	free(ineq);
	return NULL;
}

static int any(int *con, unsigned len, int status)
{
	int i;

	for (i = 0; i < len ; ++i)
		if (con[i] == status)
			return 1;
	return 0;
}

static int count(int *con, unsigned len, int status)
{
	int i;
	int c = 0;

	for (i = 0; i < len ; ++i)
		if (con[i] == status)
			c++;
	return c;
}

static int all(int *con, unsigned len, int status)
{
	int i;

	for (i = 0; i < len ; ++i) {
		if (con[i] == STATUS_REDUNDANT)
			continue;
		if (con[i] != status)
			return 0;
	}
	return 1;
}

static void drop(struct isl_map *map, int i, struct isl_tab **tabs)
{
	isl_basic_map_free(map->p[i]);
	isl_tab_free(tabs[i]);

	if (i != map->n - 1) {
		map->p[i] = map->p[map->n - 1];
		tabs[i] = tabs[map->n - 1];
	}
	tabs[map->n - 1] = NULL;
	map->n--;
}

/* Replace the pair of basic maps i and j but the basic map bounded
 * by the valid constraints in both basic maps.
 */
static int fuse(struct isl_map *map, int i, int j, struct isl_tab **tabs,
	int *ineq_i, int *ineq_j)
{
	int k, l;
	struct isl_basic_map *fused = NULL;
	struct isl_tab *fused_tab = NULL;
	unsigned total = isl_basic_map_total_dim(map->p[i]);

	fused = isl_basic_map_alloc_dim(isl_dim_copy(map->p[i]->dim),
			map->p[i]->n_div,
			map->p[i]->n_eq + map->p[j]->n_eq,
			map->p[i]->n_ineq + map->p[j]->n_ineq);
	if (!fused)
		goto error;

	for (k = 0; k < map->p[i]->n_eq; ++k) {
		int l = isl_basic_map_alloc_equality(fused);
		isl_seq_cpy(fused->eq[l], map->p[i]->eq[k], 1 + total);
	}

	for (k = 0; k < map->p[j]->n_eq; ++k) {
		int l = isl_basic_map_alloc_equality(fused);
		isl_seq_cpy(fused->eq[l], map->p[j]->eq[k], 1 + total);
	}

	for (k = 0; k < map->p[i]->n_ineq; ++k) {
		if (ineq_i[k] != STATUS_VALID)
			continue;
		l = isl_basic_map_alloc_inequality(fused);
		isl_seq_cpy(fused->ineq[l], map->p[i]->ineq[k], 1 + total);
	}

	for (k = 0; k < map->p[j]->n_ineq; ++k) {
		if (ineq_j[k] != STATUS_VALID)
			continue;
		l = isl_basic_map_alloc_inequality(fused);
		isl_seq_cpy(fused->ineq[l], map->p[j]->ineq[k], 1 + total);
	}

	for (k = 0; k < map->p[i]->n_div; ++k) {
		int l = isl_basic_map_alloc_div(fused);
		isl_seq_cpy(fused->div[l], map->p[i]->div[k], 1 + 1 + total);
	}

	fused = isl_basic_map_gauss(fused, NULL);
	ISL_F_SET(fused, ISL_BASIC_MAP_FINAL);
	if (ISL_F_ISSET(map->p[i], ISL_BASIC_MAP_RATIONAL) &&
	    ISL_F_ISSET(map->p[j], ISL_BASIC_MAP_RATIONAL))
		ISL_F_SET(fused, ISL_BASIC_MAP_RATIONAL);

	fused_tab = isl_tab_from_basic_map(fused);
	fused_tab = isl_tab_detect_redundant(fused_tab);
	if (!fused_tab)
		goto error;

	isl_basic_map_free(map->p[i]);
	map->p[i] = fused;
	isl_tab_free(tabs[i]);
	tabs[i] = fused_tab;
	drop(map, j, tabs);

	return 1;
error:
	isl_basic_map_free(fused);
	return -1;
}

/* Given a pair of basic maps i and j such that all constraints are either
 * "valid" or "cut", check if the facets corresponding to the "cut"
 * constraints of i lie entirely within basic map j.
 * If so, replace the pair by the basic map consisting of the valid
 * constraints in both basic maps.
 *
 * To see that we are not introducing any extra points, call the
 * two basic maps A and B and the resulting map U and let x
 * be an element of U \setminus ( A \cup B ).
 * Then there is a pair of cut constraints c_1 and c_2 in A and B such that x
 * violates them.  Let X be the intersection of U with the opposites
 * of these constraints.  Then x \in X.
 * The facet corresponding to c_1 contains the corresponding facet of A.
 * This facet is entirely contained in B, so c_2 is valid on the facet.
 * However, since it is also (part of) a facet of X, -c_2 is also valid
 * on the facet.  This means c_2 is saturated on the facet, so c_1 and
 * c_2 must be opposites of each other, but then x could not violate
 * both of them.
 */
static int check_facets(struct isl_map *map, int i, int j,
	struct isl_tab **tabs, int *ineq_i, int *ineq_j)
{
	int k, l;
	struct isl_tab_undo *snap;
	unsigned n_eq = map->p[i]->n_eq;

	snap = isl_tab_snap(tabs[i]);

	for (k = 0; k < map->p[i]->n_ineq; ++k) {
		if (ineq_i[k] != STATUS_CUT)
			continue;
		tabs[i] = isl_tab_select_facet(tabs[i], n_eq + k);
		for (l = 0; l < map->p[j]->n_ineq; ++l) {
			int stat;
			if (ineq_j[l] != STATUS_CUT)
				continue;
			stat = status_in(map->p[j]->ineq[l], tabs[i]);
			if (stat != STATUS_VALID)
				break;
		}
		isl_tab_rollback(tabs[i], snap);
		if (l < map->p[j]->n_ineq)
			break;
	}

	if (k < map->p[i]->n_ineq)
		/* BAD CUT PAIR */
		return 0;
	return fuse(map, i, j, tabs, ineq_i, ineq_j);
}

/* Both basic maps have at least one inequality with and adjacent
 * (but opposite) inequality in the other basic map.
 * Check that there are no cut constraints and that there is only
 * a single pair of adjacent inequalities.
 * If so, we can replace the pair by a single basic map described
 * by all but the pair of adjacent inequalities.
 * Any additional points introduced lie strictly between the two
 * adjacent hyperplanes and can therefore be integral.
 *
 *        ____			  _____
 *       /    ||\		 /     \
 *      /     || \		/       \
 *      \     ||  \	=>	\        \
 *       \    ||  /		 \       /
 *        \___||_/		  \_____/
 *
 * The test for a single pair of adjancent inequalities is important
 * for avoiding the combination of two basic maps like the following
 *
 *       /|
 *      / |
 *     /__|
 *         _____
 *         |   |
 *         |   |
 *         |___|
 */
static int check_adj_ineq(struct isl_map *map, int i, int j,
	struct isl_tab **tabs, int *ineq_i, int *ineq_j)
{
	int changed = 0;

	if (any(ineq_i, map->p[i]->n_ineq, STATUS_CUT) ||
	    any(ineq_j, map->p[j]->n_ineq, STATUS_CUT))
		/* ADJ INEQ CUT */
		;
	else if (count(ineq_i, map->p[i]->n_ineq, STATUS_ADJ_INEQ) == 1 &&
		 count(ineq_j, map->p[j]->n_ineq, STATUS_ADJ_INEQ) == 1)
		changed = fuse(map, i, j, tabs, ineq_i, ineq_j);
	/* else ADJ INEQ TOO MANY */

	return changed;
}

/* Check if basic map "i" contains the basic map represented
 * by the tableau "tab".
 */
static int contains(struct isl_map *map, int i, int *ineq_i,
	struct isl_tab *tab)
{
	int k, l;
	unsigned dim;

	dim = isl_basic_map_total_dim(map->p[i]);
	for (k = 0; k < map->p[i]->n_eq; ++k) {
		for (l = 0; l < 2; ++l) {
			int stat;
			isl_seq_neg(map->p[i]->eq[k], map->p[i]->eq[k], 1+dim);
			stat = status_in(map->p[i]->eq[k], tab);
			if (stat != STATUS_VALID)
				return 0;
		}
	}

	for (k = 0; k < map->p[i]->n_ineq; ++k) {
		int stat;
		if (ineq_i[k] == STATUS_REDUNDANT)
			continue;
		stat = status_in(map->p[i]->ineq[k], tab);
		if (stat != STATUS_VALID)
			return 0;
	}
	return 1;
}

/* At least one of the basic maps has an equality that is adjacent
 * to inequality.  Make sure that only one of the basic maps has
 * such an equality and that the other basic map has exactly one
 * inequality adjacent to an equality.
 * We call the basic map that has the inequality "i" and the basic
 * map that has the equality "j".
 * If "i" has any "cut" inequality, then relaxing the inequality
 * by one would not result in a basic map that contains the other
 * basic map.
 * Otherwise, we relax the constraint, compute the corresponding
 * facet and check whether it is included in the other basic map.
 * If so, we know that relaxing the constraint extend the basic
 * map with exactly the other basic map (we already know that this
 * other basic map is included in the extension, because there
 * were no "cut" inequalities in "i") and we can replace the
 * two basic maps by thie extension.
 *        ____			  _____
 *       /    || 		 /     |
 *      /     ||  		/      |
 *      \     ||   	=>	\      |
 *       \    ||		 \     |
 *        \___||		  \____|
 */
static int check_adj_eq(struct isl_map *map, int i, int j,
	struct isl_tab **tabs, int *eq_i, int *ineq_i, int *eq_j, int *ineq_j)
{
	int changed = 0;
	int super;
	int k;
	struct isl_tab_undo *snap, *snap2;
	unsigned n_eq = map->p[i]->n_eq;

	if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_ADJ_INEQ) &&
	    any(eq_j, 2 * map->p[j]->n_eq, STATUS_ADJ_INEQ))
		/* ADJ EQ TOO MANY */
		return 0;

	if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_ADJ_INEQ))
		return check_adj_eq(map, j, i, tabs,
					eq_j, ineq_j, eq_i, ineq_i);

	/* j has an equality adjacent to an inequality in i */

	if (any(ineq_i, map->p[i]->n_ineq, STATUS_CUT))
		/* ADJ EQ CUT */
		return 0;
	if (count(eq_j, 2 * map->p[j]->n_eq, STATUS_ADJ_INEQ) != 1 ||
	    count(ineq_i, map->p[i]->n_ineq, STATUS_ADJ_EQ) != 1 ||
	    any(ineq_j, map->p[j]->n_ineq, STATUS_ADJ_EQ) ||
	    any(ineq_i, map->p[i]->n_ineq, STATUS_ADJ_INEQ) ||
	    any(ineq_j, map->p[j]->n_ineq, STATUS_ADJ_INEQ))
		/* ADJ EQ TOO MANY */
		return 0;

	for (k = 0; k < map->p[i]->n_ineq ; ++k)
		if (ineq_i[k] == STATUS_ADJ_EQ)
			break;

	snap = isl_tab_snap(tabs[i]);
	tabs[i] = isl_tab_relax(tabs[i], n_eq + k);
	snap2 = isl_tab_snap(tabs[i]);
	tabs[i] = isl_tab_select_facet(tabs[i], n_eq + k);
	super = contains(map, j, ineq_j, tabs[i]);
	if (super) {
		isl_tab_rollback(tabs[i], snap2);
		map->p[i] = isl_basic_map_cow(map->p[i]);
		if (!map->p[i])
			return -1;
		isl_int_add_ui(map->p[i]->ineq[k][0], map->p[i]->ineq[k][0], 1);
		ISL_F_SET(map->p[i], ISL_BASIC_MAP_FINAL);
		drop(map, j, tabs);
		changed = 1;
	} else
		isl_tab_rollback(tabs[i], snap);

	return changed;
}

/* Check if the union of the given pair of basic maps
 * can be represented by a single basic map.
 * If so, replace the pair by the single basic map and return 1.
 * Otherwise, return 0;
 *
 * We first check the effect of each constraint of one basic map
 * on the other basic map.
 * The constraint may be
 *	redundant	the constraint is redundant in its own
 *			basic map and should be ignore and removed
 *			in the end
 *	valid		all (integer) points of the other basic map
 *			satisfy the constraint
 *	separate	no (integer) point of the other basic map
 *			satisfies the constraint
 *	cut		some but not all points of the other basic map
 *			satisfy the constraint
 *	adj_eq		the given constraint is adjacent (on the outside)
 *			to an equality of the other basic map
 *	adj_ineq	the given constraint is adjacent (on the outside)
 *			to an inequality of the other basic map
 *
 * We consider four cases in which we can replace the pair by a single
 * basic map.  We ignore all "redundant" constraints.
 *
 *	1. all constraints of one basic map are valid
 *		=> the other basic map is a subset and can be removed
 *
 *	2. all constraints of both basic maps are either "valid" or "cut"
 *	   and the facets corresponding to the "cut" constraints
 *	   of one of the basic maps lies entirely inside the other basic map
 *		=> the pair can be replaced by a basic map consisting
 *		   of the valid constraints in both basic maps
 *
 *	3. there is a single pair of adjacent inequalities
 *	   (all other constraints are "valid")
 *		=> the pair can be replaced by a basic map consisting
 *		   of the valid constraints in both basic maps
 *
 *	4. there is a single adjacent pair of an inequality and an equality,
 *	   the other constraints of the basic map containing the inequality are
 *	   "valid".  Moreover, if the inequality the basic map is relaxed
 *	   and then turned into an equality, then resulting facet lies
 *	   entirely inside the other basic map
 *		=> the pair can be replaced by the basic map containing
 *		   the inequality, with the inequality relaxed.
 *
 * Throughout the computation, we maintain a collection of tableaus
 * corresponding to the basic maps.  When the basic maps are dropped
 * or combined, the tableaus are modified accordingly.
 */
static int coalesce_pair(struct isl_map *map, int i, int j,
	struct isl_tab **tabs)
{
	int changed = 0;
	int *eq_i = NULL;
	int *eq_j = NULL;
	int *ineq_i = NULL;
	int *ineq_j = NULL;

	eq_i = eq_status_in(map, i, j, tabs);
	if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_ERROR))
		goto error;
	if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_SEPARATE))
		goto done;

	eq_j = eq_status_in(map, j, i, tabs);
	if (any(eq_j, 2 * map->p[j]->n_eq, STATUS_ERROR))
		goto error;
	if (any(eq_j, 2 * map->p[j]->n_eq, STATUS_SEPARATE))
		goto done;

	ineq_i = ineq_status_in(map, i, j, tabs);
	if (any(ineq_i, map->p[i]->n_ineq, STATUS_ERROR))
		goto error;
	if (any(ineq_i, map->p[i]->n_ineq, STATUS_SEPARATE))
		goto done;

	ineq_j = ineq_status_in(map, j, i, tabs);
	if (any(ineq_j, map->p[j]->n_ineq, STATUS_ERROR))
		goto error;
	if (any(ineq_j, map->p[j]->n_ineq, STATUS_SEPARATE))
		goto done;

	if (all(eq_i, 2 * map->p[i]->n_eq, STATUS_VALID) &&
	    all(ineq_i, map->p[i]->n_ineq, STATUS_VALID)) {
		drop(map, j, tabs);
		changed = 1;
	} else if (all(eq_j, 2 * map->p[j]->n_eq, STATUS_VALID) &&
		   all(ineq_j, map->p[j]->n_ineq, STATUS_VALID)) {
		drop(map, i, tabs);
		changed = 1;
	} else if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_CUT) ||
		   any(eq_j, 2 * map->p[j]->n_eq, STATUS_CUT)) {
		/* BAD CUT */
	} else if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_ADJ_EQ) ||
		   any(eq_j, 2 * map->p[j]->n_eq, STATUS_ADJ_EQ)) {
		/* ADJ EQ PAIR */
	} else if (any(eq_i, 2 * map->p[i]->n_eq, STATUS_ADJ_INEQ) ||
		   any(eq_j, 2 * map->p[j]->n_eq, STATUS_ADJ_INEQ)) {
		changed = check_adj_eq(map, i, j, tabs,
					eq_i, ineq_i, eq_j, ineq_j);
	} else if (any(ineq_i, map->p[i]->n_ineq, STATUS_ADJ_EQ) ||
		   any(ineq_j, map->p[j]->n_ineq, STATUS_ADJ_EQ)) {
		/* Can't happen */
		/* BAD ADJ INEQ */
	} else if (any(ineq_i, map->p[i]->n_ineq, STATUS_ADJ_INEQ) ||
		   any(ineq_j, map->p[j]->n_ineq, STATUS_ADJ_INEQ)) {
		changed = check_adj_ineq(map, i, j, tabs, ineq_i, ineq_j);
	} else
		changed = check_facets(map, i, j, tabs, ineq_i, ineq_j);

done:
	free(eq_i);
	free(eq_j);
	free(ineq_i);
	free(ineq_j);
	return changed;
error:
	free(eq_i);
	free(eq_j);
	free(ineq_i);
	free(ineq_j);
	return -1;
}

static struct isl_map *coalesce(struct isl_map *map, struct isl_tab **tabs)
{
	int i, j;

	for (i = 0; i < map->n - 1; ++i)
		for (j = i + 1; j < map->n; ++j) {
			int changed;
			changed = coalesce_pair(map, i, j, tabs);
			if (changed < 0)
				goto error;
			if (changed)
				return coalesce(map, tabs);
		}
	return map;
error:
	isl_map_free(map);
	return NULL;
}

/* For each pair of basic maps in the map, check if the union of the two
 * can be represented by a single basic map.
 * If so, replace the pair by the single basic map and start over.
 */
struct isl_map *isl_map_coalesce(struct isl_map *map)
{
	int i;
	unsigned n;
	struct isl_tab **tabs = NULL;

	if (!map)
		return NULL;

	if (map->n <= 1)
		return map;

	map = isl_map_align_divs(map);

	tabs = isl_calloc_array(map->ctx, struct isl_tab *, map->n);
	if (!tabs)
		goto error;

	n = map->n;
	for (i = 0; i < map->n; ++i) {
		tabs[i] = isl_tab_from_basic_map(map->p[i]);
		if (!tabs[i])
			goto error;
		if (!ISL_F_ISSET(map->p[i], ISL_BASIC_MAP_NO_IMPLICIT))
			tabs[i] = isl_tab_detect_equalities(tabs[i]);
		if (!ISL_F_ISSET(map->p[i], ISL_BASIC_MAP_NO_REDUNDANT))
			tabs[i] = isl_tab_detect_redundant(tabs[i]);
	}
	for (i = map->n - 1; i >= 0; --i)
		if (tabs[i]->empty)
			drop(map, i, tabs);

	map = coalesce(map, tabs);

	if (map)
		for (i = 0; i < map->n; ++i) {
			map->p[i] = isl_basic_map_update_from_tab(map->p[i],
								    tabs[i]);
			map->p[i] = isl_basic_map_finalize(map->p[i]);
			if (!map->p[i])
				goto error;
			ISL_F_SET(map->p[i], ISL_BASIC_MAP_NO_IMPLICIT);
			ISL_F_SET(map->p[i], ISL_BASIC_MAP_NO_REDUNDANT);
		}

	for (i = 0; i < n; ++i)
		isl_tab_free(tabs[i]);

	free(tabs);

	return map;
error:
	if (tabs)
		for (i = 0; i < n; ++i)
			isl_tab_free(tabs[i]);
	free(tabs);
	return NULL;
}

/* For each pair of basic sets in the set, check if the union of the two
 * can be represented by a single basic set.
 * If so, replace the pair by the single basic set and start over.
 */
struct isl_set *isl_set_coalesce(struct isl_set *set)
{
	(struct isl_set *)isl_map_coalesce((struct isl_map *)set);
}
