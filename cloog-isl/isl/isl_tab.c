#include "isl_map_private.h"
#include "isl_tab.h"

/*
 * The implementation of tableaus in this file was inspired by Section 8
 * of David Detlefs, Greg Nelson and James B. Saxe, "Simplify: a theorem
 * prover for program checking".
 */

struct isl_tab *isl_tab_alloc(struct isl_ctx *ctx,
	unsigned n_row, unsigned n_var)
{
	int i;
	struct isl_tab *tab;

	tab = isl_calloc_type(ctx, struct isl_tab);
	if (!tab)
		return NULL;
	tab->mat = isl_mat_alloc(ctx, n_row, 2 + n_var);
	if (!tab->mat)
		goto error;
	tab->var = isl_alloc_array(ctx, struct isl_tab_var, n_var);
	if (!tab->var)
		goto error;
	tab->con = isl_alloc_array(ctx, struct isl_tab_var, n_row);
	if (!tab->con)
		goto error;
	tab->col_var = isl_alloc_array(ctx, int, n_var);
	if (!tab->col_var)
		goto error;
	tab->row_var = isl_alloc_array(ctx, int, n_row);
	if (!tab->row_var)
		goto error;
	for (i = 0; i < n_var; ++i) {
		tab->var[i].index = i;
		tab->var[i].is_row = 0;
		tab->var[i].is_nonneg = 0;
		tab->var[i].is_zero = 0;
		tab->var[i].is_redundant = 0;
		tab->var[i].frozen = 0;
		tab->col_var[i] = i;
	}
	tab->n_row = 0;
	tab->n_con = 0;
	tab->n_eq = 0;
	tab->max_con = n_row;
	tab->n_col = n_var;
	tab->n_var = n_var;
	tab->n_dead = 0;
	tab->n_redundant = 0;
	tab->need_undo = 0;
	tab->rational = 0;
	tab->empty = 0;
	tab->in_undo = 0;
	tab->bottom.type = isl_tab_undo_bottom;
	tab->bottom.next = NULL;
	tab->top = &tab->bottom;
	return tab;
error:
	isl_tab_free(tab);
	return NULL;
}

static int extend_cons(struct isl_tab *tab, unsigned n_new)
{
	if (tab->max_con < tab->n_con + n_new) {
		struct isl_tab_var *con;

		con = isl_realloc_array(tab->mat->ctx, tab->con,
				    struct isl_tab_var, tab->max_con + n_new);
		if (!con)
			return -1;
		tab->con = con;
		tab->max_con += n_new;
	}
	if (tab->mat->n_row < tab->n_row + n_new) {
		int *row_var;

		tab->mat = isl_mat_extend(tab->mat,
						tab->n_row + n_new, tab->n_col);
		if (!tab->mat)
			return -1;
		row_var = isl_realloc_array(tab->mat->ctx, tab->row_var,
					    int, tab->mat->n_row);
		if (!row_var)
			return -1;
		tab->row_var = row_var;
	}
	return 0;
}

struct isl_tab *isl_tab_extend(struct isl_tab *tab, unsigned n_new)
{
	if (extend_cons(tab, n_new) >= 0)
		return tab;

	isl_tab_free(tab);
	return NULL;
}

static void free_undo(struct isl_tab *tab)
{
	struct isl_tab_undo *undo, *next;

	for (undo = tab->top; undo && undo != &tab->bottom; undo = next) {
		next = undo->next;
		free(undo);
	}
	tab->top = undo;
}

void isl_tab_free(struct isl_tab *tab)
{
	if (!tab)
		return;
	free_undo(tab);
	isl_mat_free(tab->mat);
	isl_vec_free(tab->dual);
	free(tab->var);
	free(tab->con);
	free(tab->row_var);
	free(tab->col_var);
	free(tab);
}

static struct isl_tab_var *var_from_index(struct isl_tab *tab, int i)
{
	if (i >= 0)
		return &tab->var[i];
	else
		return &tab->con[~i];
}

static struct isl_tab_var *var_from_row(struct isl_tab *tab, int i)
{
	return var_from_index(tab, tab->row_var[i]);
}

static struct isl_tab_var *var_from_col(struct isl_tab *tab, int i)
{
	return var_from_index(tab, tab->col_var[i]);
}

/* Check if there are any upper bounds on column variable "var",
 * i.e., non-negative rows where var appears with a negative coefficient.
 * Return 1 if there are no such bounds.
 */
static int max_is_manifestly_unbounded(struct isl_tab *tab,
	struct isl_tab_var *var)
{
	int i;

	if (var->is_row)
		return 0;
	for (i = tab->n_redundant; i < tab->n_row; ++i) {
		if (!isl_int_is_neg(tab->mat->row[i][2 + var->index]))
			continue;
		if (var_from_row(tab, i)->is_nonneg)
			return 0;
	}
	return 1;
}

/* Check if there are any lower bounds on column variable "var",
 * i.e., non-negative rows where var appears with a positive coefficient.
 * Return 1 if there are no such bounds.
 */
static int min_is_manifestly_unbounded(struct isl_tab *tab,
	struct isl_tab_var *var)
{
	int i;

	if (var->is_row)
		return 0;
	for (i = tab->n_redundant; i < tab->n_row; ++i) {
		if (!isl_int_is_pos(tab->mat->row[i][2 + var->index]))
			continue;
		if (var_from_row(tab, i)->is_nonneg)
			return 0;
	}
	return 1;
}

/* Given the index of a column "c", return the index of a row
 * that can be used to pivot the column in, with either an increase
 * (sgn > 0) or a decrease (sgn < 0) of the corresponding variable.
 * If "var" is not NULL, then the row returned will be different from
 * the one associated with "var".
 *
 * Each row in the tableau is of the form
 *
 *	x_r = a_r0 + \sum_i a_ri x_i
 *
 * Only rows with x_r >= 0 and with the sign of a_ri opposite to "sgn"
 * impose any limit on the increase or decrease in the value of x_c
 * and this bound is equal to a_r0 / |a_rc|.  We are therefore looking
 * for the row with the smallest (most stringent) such bound.
 * Note that the common denominator of each row drops out of the fraction.
 * To check if row j has a smaller bound than row r, i.e.,
 * a_j0 / |a_jc| < a_r0 / |a_rc| or a_j0 |a_rc| < a_r0 |a_jc|,
 * we check if -sign(a_jc) (a_j0 a_rc - a_r0 a_jc) < 0,
 * where -sign(a_jc) is equal to "sgn".
 */
static int pivot_row(struct isl_tab *tab,
	struct isl_tab_var *var, int sgn, int c)
{
	int j, r, tsgn;
	isl_int t;

	isl_int_init(t);
	r = -1;
	for (j = tab->n_redundant; j < tab->n_row; ++j) {
		if (var && j == var->index)
			continue;
		if (!var_from_row(tab, j)->is_nonneg)
			continue;
		if (sgn * isl_int_sgn(tab->mat->row[j][2 + c]) >= 0)
			continue;
		if (r < 0) {
			r = j;
			continue;
		}
		isl_int_mul(t, tab->mat->row[r][1], tab->mat->row[j][2 + c]);
		isl_int_submul(t, tab->mat->row[j][1], tab->mat->row[r][2 + c]);
		tsgn = sgn * isl_int_sgn(t);
		if (tsgn < 0 || (tsgn == 0 &&
					    tab->row_var[j] < tab->row_var[r]))
			r = j;
	}
	isl_int_clear(t);
	return r;
}

/* Find a pivot (row and col) that will increase (sgn > 0) or decrease
 * (sgn < 0) the value of row variable var.
 * If not NULL, then skip_var is a row variable that should be ignored
 * while looking for a pivot row.  It is usually equal to var.
 *
 * As the given row in the tableau is of the form
 *
 *	x_r = a_r0 + \sum_i a_ri x_i
 *
 * we need to find a column such that the sign of a_ri is equal to "sgn"
 * (such that an increase in x_i will have the desired effect) or a
 * column with a variable that may attain negative values.
 * If a_ri is positive, then we need to move x_i in the same direction
 * to obtain the desired effect.  Otherwise, x_i has to move in the
 * opposite direction.
 */
static void find_pivot(struct isl_tab *tab,
	struct isl_tab_var *var, struct isl_tab_var *skip_var,
	int sgn, int *row, int *col)
{
	int j, r, c;
	isl_int *tr;

	*row = *col = -1;

	isl_assert(tab->mat->ctx, var->is_row, return);
	tr = tab->mat->row[var->index];

	c = -1;
	for (j = tab->n_dead; j < tab->n_col; ++j) {
		if (isl_int_is_zero(tr[2 + j]))
			continue;
		if (isl_int_sgn(tr[2 + j]) != sgn &&
		    var_from_col(tab, j)->is_nonneg)
			continue;
		if (c < 0 || tab->col_var[j] < tab->col_var[c])
			c = j;
	}
	if (c < 0)
		return;

	sgn *= isl_int_sgn(tr[2 + c]);
	r = pivot_row(tab, skip_var, sgn, c);
	*row = r < 0 ? var->index : r;
	*col = c;
}

/* Return 1 if row "row" represents an obviously redundant inequality.
 * This means
 *	- it represents an inequality or a variable
 *	- that is the sum of a non-negative sample value and a positive
 *	  combination of zero or more non-negative variables.
 */
static int is_redundant(struct isl_tab *tab, int row)
{
	int i;

	if (tab->row_var[row] < 0 && !var_from_row(tab, row)->is_nonneg)
		return 0;

	if (isl_int_is_neg(tab->mat->row[row][1]))
		return 0;

	for (i = tab->n_dead; i < tab->n_col; ++i) {
		if (isl_int_is_zero(tab->mat->row[row][2 + i]))
			continue;
		if (isl_int_is_neg(tab->mat->row[row][2 + i]))
			return 0;
		if (!var_from_col(tab, i)->is_nonneg)
			return 0;
	}
	return 1;
}

static void swap_rows(struct isl_tab *tab, int row1, int row2)
{
	int t;
	t = tab->row_var[row1];
	tab->row_var[row1] = tab->row_var[row2];
	tab->row_var[row2] = t;
	var_from_row(tab, row1)->index = row1;
	var_from_row(tab, row2)->index = row2;
	tab->mat = isl_mat_swap_rows(tab->mat, row1, row2);
}

static void push(struct isl_tab *tab,
	enum isl_tab_undo_type type, struct isl_tab_var *var)
{
	struct isl_tab_undo *undo;

	if (!tab->need_undo)
		return;

	undo = isl_alloc_type(tab->mat->ctx, struct isl_tab_undo);
	if (!undo) {
		free_undo(tab);
		tab->top = NULL;
		return;
	}
	undo->type = type;
	undo->var = var;
	undo->next = tab->top;
	tab->top = undo;
}

/* Mark row with index "row" as being redundant.
 * If we may need to undo the operation or if the row represents
 * a variable of the original problem, the row is kept,
 * but no longer considered when looking for a pivot row.
 * Otherwise, the row is simply removed.
 *
 * The row may be interchanged with some other row.  If it
 * is interchanged with a later row, return 1.  Otherwise return 0.
 * If the rows are checked in order in the calling function,
 * then a return value of 1 means that the row with the given
 * row number may now contain a different row that hasn't been checked yet.
 */
static int mark_redundant(struct isl_tab *tab, int row)
{
	struct isl_tab_var *var = var_from_row(tab, row);
	var->is_redundant = 1;
	isl_assert(tab->mat->ctx, row >= tab->n_redundant, return);
	if (tab->need_undo || tab->row_var[row] >= 0) {
		if (tab->row_var[row] >= 0) {
			var->is_nonneg = 1;
			push(tab, isl_tab_undo_nonneg, var);
		}
		if (row != tab->n_redundant)
			swap_rows(tab, row, tab->n_redundant);
		push(tab, isl_tab_undo_redundant, var);
		tab->n_redundant++;
		return 0;
	} else {
		if (row != tab->n_row - 1)
			swap_rows(tab, row, tab->n_row - 1);
		var_from_row(tab, tab->n_row - 1)->index = -1;
		tab->n_row--;
		return 1;
	}
}

static void mark_empty(struct isl_tab *tab)
{
	if (!tab->empty && tab->need_undo)
		push(tab, isl_tab_undo_empty, NULL);
	tab->empty = 1;
}

/* Given a row number "row" and a column number "col", pivot the tableau
 * such that the associated variables are interchanged.
 * The given row in the tableau expresses
 *
 *	x_r = a_r0 + \sum_i a_ri x_i
 *
 * or
 *
 *	x_c = 1/a_rc x_r - a_r0/a_rc + sum_{i \ne r} -a_ri/a_rc
 *
 * Substituting this equality into the other rows
 *
 *	x_j = a_j0 + \sum_i a_ji x_i
 *
 * with a_jc \ne 0, we obtain
 *
 *	x_j = a_jc/a_rc x_r + a_j0 - a_jc a_r0/a_rc + sum a_ji - a_jc a_ri/a_rc 
 *
 * The tableau
 *
 *	n_rc/d_r		n_ri/d_r
 *	n_jc/d_j		n_ji/d_j
 *
 * where i is any other column and j is any other row,
 * is therefore transformed into
 *
 * s(n_rc)d_r/|n_rc|		-s(n_rc)n_ri/|n_rc|
 * s(n_rc)d_r n_jc/(|n_rc| d_j)	(n_ji |n_rc| - s(n_rc)n_jc n_ri)/(|n_rc| d_j)
 *
 * The transformation is performed along the following steps
 *
 *	d_r/n_rc		n_ri/n_rc
 *	n_jc/d_j		n_ji/d_j
 *
 *	s(n_rc)d_r/|n_rc|	-s(n_rc)n_ri/|n_rc|
 *	n_jc/d_j		n_ji/d_j
 *
 *	s(n_rc)d_r/|n_rc|	-s(n_rc)n_ri/|n_rc|
 *	n_jc/(|n_rc| d_j)	n_ji/(|n_rc| d_j)
 *
 *	s(n_rc)d_r/|n_rc|	-s(n_rc)n_ri/|n_rc|
 *	n_jc/(|n_rc| d_j)	(n_ji |n_rc|)/(|n_rc| d_j)
 *
 *	s(n_rc)d_r/|n_rc|	-s(n_rc)n_ri/|n_rc|
 *	n_jc/(|n_rc| d_j)	(n_ji |n_rc| - s(n_rc)n_jc n_ri)/(|n_rc| d_j)
 *
 * s(n_rc)d_r/|n_rc|		-s(n_rc)n_ri/|n_rc|
 * s(n_rc)d_r n_jc/(|n_rc| d_j)	(n_ji |n_rc| - s(n_rc)n_jc n_ri)/(|n_rc| d_j)
 *
 */
static void pivot(struct isl_tab *tab, int row, int col)
{
	int i, j;
	int sgn;
	int t;
	struct isl_mat *mat = tab->mat;
	struct isl_tab_var *var;

	isl_int_swap(mat->row[row][0], mat->row[row][2 + col]);
	sgn = isl_int_sgn(mat->row[row][0]);
	if (sgn < 0) {
		isl_int_neg(mat->row[row][0], mat->row[row][0]);
		isl_int_neg(mat->row[row][2 + col], mat->row[row][2 + col]);
	} else
		for (j = 0; j < 1 + tab->n_col; ++j) {
			if (j == 1 + col)
				continue;
			isl_int_neg(mat->row[row][1 + j], mat->row[row][1 + j]);
		}
	if (!isl_int_is_one(mat->row[row][0]))
		isl_seq_normalize(mat->row[row], 2 + tab->n_col);
	for (i = 0; i < tab->n_row; ++i) {
		if (i == row)
			continue;
		if (isl_int_is_zero(mat->row[i][2 + col]))
			continue;
		isl_int_mul(mat->row[i][0], mat->row[i][0], mat->row[row][0]);
		for (j = 0; j < 1 + tab->n_col; ++j) {
			if (j == 1 + col)
				continue;
			isl_int_mul(mat->row[i][1 + j],
				    mat->row[i][1 + j], mat->row[row][0]);
			isl_int_addmul(mat->row[i][1 + j],
				    mat->row[i][2 + col], mat->row[row][1 + j]);
		}
		isl_int_mul(mat->row[i][2 + col],
			    mat->row[i][2 + col], mat->row[row][2 + col]);
		if (!isl_int_is_one(mat->row[row][0]))
			isl_seq_normalize(mat->row[i], 2 + tab->n_col);
	}
	t = tab->row_var[row];
	tab->row_var[row] = tab->col_var[col];
	tab->col_var[col] = t;
	var = var_from_row(tab, row);
	var->is_row = 1;
	var->index = row;
	var = var_from_col(tab, col);
	var->is_row = 0;
	var->index = col;
	if (tab->in_undo)
		return;
	for (i = tab->n_redundant; i < tab->n_row; ++i) {
		if (isl_int_is_zero(mat->row[i][2 + col]))
			continue;
		if (!var_from_row(tab, i)->frozen &&
		    is_redundant(tab, i))
			if (mark_redundant(tab, i))
				--i;
	}
}

/* If "var" represents a column variable, then pivot is up (sgn > 0)
 * or down (sgn < 0) to a row.  The variable is assumed not to be
 * unbounded in the specified direction.
 */
static void to_row(struct isl_tab *tab, struct isl_tab_var *var, int sign)
{
	int r;

	if (var->is_row)
		return;

	r = pivot_row(tab, NULL, sign, var->index);
	isl_assert(tab->mat->ctx, r >= 0, return);
	pivot(tab, r, var->index);
}

static void check_table(struct isl_tab *tab)
{
	int i;

	if (tab->empty)
		return;
	for (i = 0; i < tab->n_row; ++i) {
		if (!var_from_row(tab, i)->is_nonneg)
			continue;
		assert(!isl_int_is_neg(tab->mat->row[i][1]));
	}
}

/* Return the sign of the maximal value of "var".
 * If the sign is not negative, then on return from this function,
 * the sample value will also be non-negative.
 *
 * If "var" is manifestly unbounded wrt positive values, we are done.
 * Otherwise, we pivot the variable up to a row if needed
 * Then we continue pivoting down until either
 *	- no more down pivots can be performed
 *	- the sample value is positive
 *	- the variable is pivoted into a manifestly unbounded column
 */
static int sign_of_max(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;

	if (max_is_manifestly_unbounded(tab, var))
		return 1;
	to_row(tab, var, 1);
	while (!isl_int_is_pos(tab->mat->row[var->index][1])) {
		find_pivot(tab, var, var, 1, &row, &col);
		if (row == -1)
			return isl_int_sgn(tab->mat->row[var->index][1]);
		pivot(tab, row, col);
		if (!var->is_row) /* manifestly unbounded */
			return 1;
	}
	return 1;
}

/* Perform pivots until the row variable "var" has a non-negative
 * sample value or until no more upward pivots can be performed.
 * Return the sign of the sample value after the pivots have been
 * performed.
 */
static int restore_row(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;

	while (isl_int_is_neg(tab->mat->row[var->index][1])) {
		find_pivot(tab, var, var, 1, &row, &col);
		if (row == -1)
			break;
		pivot(tab, row, col);
		if (!var->is_row) /* manifestly unbounded */
			return 1;
	}
	return isl_int_sgn(tab->mat->row[var->index][1]);
}

/* Perform pivots until we are sure that the row variable "var"
 * can attain non-negative values.  After return from this
 * function, "var" is still a row variable, but its sample
 * value may not be non-negative, even if the function returns 1.
 */
static int at_least_zero(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;

	while (isl_int_is_neg(tab->mat->row[var->index][1])) {
		find_pivot(tab, var, var, 1, &row, &col);
		if (row == -1)
			break;
		if (row == var->index) /* manifestly unbounded */
			return 1;
		pivot(tab, row, col);
	}
	return !isl_int_is_neg(tab->mat->row[var->index][1]);
}

/* Return a negative value if "var" can attain negative values.
 * Return a non-negative value otherwise.
 *
 * If "var" is manifestly unbounded wrt negative values, we are done.
 * Otherwise, if var is in a column, we can pivot it down to a row.
 * Then we continue pivoting down until either
 *	- the pivot would result in a manifestly unbounded column
 *	  => we don't perform the pivot, but simply return -1
 *	- no more down pivots can be performed
 *	- the sample value is negative
 * If the sample value becomes negative and the variable is supposed
 * to be nonnegative, then we undo the last pivot.
 * However, if the last pivot has made the pivoting variable
 * obviously redundant, then it may have moved to another row.
 * In that case we look for upward pivots until we reach a non-negative
 * value again.
 */
static int sign_of_min(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;
	struct isl_tab_var *pivot_var;

	if (min_is_manifestly_unbounded(tab, var))
		return -1;
	if (!var->is_row) {
		col = var->index;
		row = pivot_row(tab, NULL, -1, col);
		pivot_var = var_from_col(tab, col);
		pivot(tab, row, col);
		if (var->is_redundant)
			return 0;
		if (isl_int_is_neg(tab->mat->row[var->index][1])) {
			if (var->is_nonneg) {
				if (!pivot_var->is_redundant &&
				    pivot_var->index == row)
					pivot(tab, row, col);
				else
					restore_row(tab, var);
			}
			return -1;
		}
	}
	if (var->is_redundant)
		return 0;
	while (!isl_int_is_neg(tab->mat->row[var->index][1])) {
		find_pivot(tab, var, var, -1, &row, &col);
		if (row == var->index)
			return -1;
		if (row == -1)
			return isl_int_sgn(tab->mat->row[var->index][1]);
		pivot_var = var_from_col(tab, col);
		pivot(tab, row, col);
		if (var->is_redundant)
			return 0;
	}
	if (var->is_nonneg) {
		/* pivot back to non-negative value */
		if (!pivot_var->is_redundant && pivot_var->index == row)
			pivot(tab, row, col);
		else
			restore_row(tab, var);
	}
	return -1;
}

/* Return 1 if "var" can attain values <= -1.
 * Return 0 otherwise.
 *
 * The sample value of "var" is assumed to be non-negative when the
 * the function is called and will be made non-negative again before
 * the function returns.
 */
static int min_at_most_neg_one(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;
	struct isl_tab_var *pivot_var;

	if (min_is_manifestly_unbounded(tab, var))
		return 1;
	if (!var->is_row) {
		col = var->index;
		row = pivot_row(tab, NULL, -1, col);
		pivot_var = var_from_col(tab, col);
		pivot(tab, row, col);
		if (var->is_redundant)
			return 0;
		if (isl_int_is_neg(tab->mat->row[var->index][1]) &&
		    isl_int_abs_ge(tab->mat->row[var->index][1],
				   tab->mat->row[var->index][0])) {
			if (var->is_nonneg) {
				if (!pivot_var->is_redundant &&
				    pivot_var->index == row)
					pivot(tab, row, col);
				else
					restore_row(tab, var);
			}
			return 1;
		}
	}
	if (var->is_redundant)
		return 0;
	do {
		find_pivot(tab, var, var, -1, &row, &col);
		if (row == var->index)
			return 1;
		if (row == -1)
			return 0;
		pivot_var = var_from_col(tab, col);
		pivot(tab, row, col);
		if (var->is_redundant)
			return 0;
	} while (!isl_int_is_neg(tab->mat->row[var->index][1]) ||
		 isl_int_abs_lt(tab->mat->row[var->index][1],
				tab->mat->row[var->index][0]));
	if (var->is_nonneg) {
		/* pivot back to non-negative value */
		if (!pivot_var->is_redundant && pivot_var->index == row)
			pivot(tab, row, col);
		restore_row(tab, var);
	}
	return 1;
}

/* Return 1 if "var" can attain values >= 1.
 * Return 0 otherwise.
 */
static int at_least_one(struct isl_tab *tab, struct isl_tab_var *var)
{
	int row, col;
	isl_int *r;

	if (max_is_manifestly_unbounded(tab, var))
		return 1;
	to_row(tab, var, 1);
	r = tab->mat->row[var->index];
	while (isl_int_lt(r[1], r[0])) {
		find_pivot(tab, var, var, 1, &row, &col);
		if (row == -1)
			return isl_int_ge(r[1], r[0]);
		if (row == var->index) /* manifestly unbounded */
			return 1;
		pivot(tab, row, col);
	}
	return 1;
}

static void swap_cols(struct isl_tab *tab, int col1, int col2)
{
	int t;
	t = tab->col_var[col1];
	tab->col_var[col1] = tab->col_var[col2];
	tab->col_var[col2] = t;
	var_from_col(tab, col1)->index = col1;
	var_from_col(tab, col2)->index = col2;
	tab->mat = isl_mat_swap_cols(tab->mat, 2 + col1, 2 + col2);
}

/* Mark column with index "col" as representing a zero variable.
 * If we may need to undo the operation the column is kept,
 * but no longer considered.
 * Otherwise, the column is simply removed.
 *
 * The column may be interchanged with some other column.  If it
 * is interchanged with a later column, return 1.  Otherwise return 0.
 * If the columns are checked in order in the calling function,
 * then a return value of 1 means that the column with the given
 * column number may now contain a different column that
 * hasn't been checked yet.
 */
static int kill_col(struct isl_tab *tab, int col)
{
	var_from_col(tab, col)->is_zero = 1;
	if (tab->need_undo) {
		push(tab, isl_tab_undo_zero, var_from_col(tab, col));
		if (col != tab->n_dead)
			swap_cols(tab, col, tab->n_dead);
		tab->n_dead++;
		return 0;
	} else {
		if (col != tab->n_col - 1)
			swap_cols(tab, col, tab->n_col - 1);
		var_from_col(tab, tab->n_col - 1)->index = -1;
		tab->n_col--;
		return 1;
	}
}

/* Row variable "var" is non-negative and cannot attain any values
 * larger than zero.  This means that the coefficients of the unrestricted
 * column variables are zero and that the coefficients of the non-negative
 * column variables are zero or negative.
 * Each of the non-negative variables with a negative coefficient can
 * then also be written as the negative sum of non-negative variables
 * and must therefore also be zero.
 */
static void close_row(struct isl_tab *tab, struct isl_tab_var *var)
{
	int j;
	struct isl_mat *mat = tab->mat;

	isl_assert(tab->mat->ctx, var->is_nonneg, return);
	var->is_zero = 1;
	for (j = tab->n_dead; j < tab->n_col; ++j) {
		if (isl_int_is_zero(mat->row[var->index][2 + j]))
			continue;
		isl_assert(tab->mat->ctx,
			isl_int_is_neg(mat->row[var->index][2 + j]), return);
		if (kill_col(tab, j))
			--j;
	}
	mark_redundant(tab, var->index);
}

/* Add a row to the tableau.  The row is given as an affine combination
 * of the original variables and needs to be expressed in terms of the
 * column variables.
 *
 * We add each term in turn.
 * If r = n/d_r is the current sum and we need to add k x, then
 * 	if x is a column variable, we increase the numerator of
 *		this column by k d_r
 *	if x = f/d_x is a row variable, then the new representation of r is
 *
 *		 n    k f   d_x/g n + d_r/g k f   m/d_r n + m/d_g k f
 *		--- + --- = ------------------- = -------------------
 *		d_r   d_r        d_r d_x/g                m
 *
 *	with g the gcd of d_r and d_x and m the lcm of d_r and d_x.
 */
static int add_row(struct isl_tab *tab, isl_int *line)
{
	int i;
	unsigned r;
	isl_int *row;
	isl_int a, b;

	isl_assert(tab->mat->ctx, tab->n_row < tab->mat->n_row, return -1);

	isl_int_init(a);
	isl_int_init(b);
	r = tab->n_con;
	tab->con[r].index = tab->n_row;
	tab->con[r].is_row = 1;
	tab->con[r].is_nonneg = 0;
	tab->con[r].is_zero = 0;
	tab->con[r].is_redundant = 0;
	tab->con[r].frozen = 0;
	tab->row_var[tab->n_row] = ~r;
	row = tab->mat->row[tab->n_row];
	isl_int_set_si(row[0], 1);
	isl_int_set(row[1], line[0]);
	isl_seq_clr(row + 2, tab->n_col);
	for (i = 0; i < tab->n_var; ++i) {
		if (tab->var[i].is_zero)
			continue;
		if (tab->var[i].is_row) {
			isl_int_lcm(a,
				row[0], tab->mat->row[tab->var[i].index][0]);
			isl_int_swap(a, row[0]);
			isl_int_divexact(a, row[0], a);
			isl_int_divexact(b,
				row[0], tab->mat->row[tab->var[i].index][0]);
			isl_int_mul(b, b, line[1 + i]);
			isl_seq_combine(row + 1, a, row + 1,
			    b, tab->mat->row[tab->var[i].index] + 1,
			    1 + tab->n_col);
		} else
			isl_int_addmul(row[2 + tab->var[i].index],
							line[1 + i], row[0]);
	}
	isl_seq_normalize(row, 2 + tab->n_col);
	tab->n_row++;
	tab->n_con++;
	push(tab, isl_tab_undo_allocate, &tab->con[r]);
	isl_int_clear(a);
	isl_int_clear(b);

	return r;
}

static int drop_row(struct isl_tab *tab, int row)
{
	isl_assert(tab->mat->ctx, ~tab->row_var[row] == tab->n_con - 1, return -1);
	if (row != tab->n_row - 1)
		swap_rows(tab, row, tab->n_row - 1);
	tab->n_row--;
	tab->n_con--;
	return 0;
}

/* Add inequality "ineq" and check if it conflicts with the
 * previously added constraints or if it is obviously redundant.
 */
struct isl_tab *isl_tab_add_ineq(struct isl_tab *tab, isl_int *ineq)
{
	int r;
	int sgn;

	if (!tab)
		return NULL;
	r = add_row(tab, ineq);
	if (r < 0)
		goto error;
	tab->con[r].is_nonneg = 1;
	push(tab, isl_tab_undo_nonneg, &tab->con[r]);
	if (is_redundant(tab, tab->con[r].index)) {
		mark_redundant(tab, tab->con[r].index);
		return tab;
	}

	sgn = restore_row(tab, &tab->con[r]);
	if (sgn < 0)
		mark_empty(tab);
	else if (tab->con[r].is_row &&
		 is_redundant(tab, tab->con[r].index))
		mark_redundant(tab, tab->con[r].index);
	return tab;
error:
	isl_tab_free(tab);
	return NULL;
}

/* Pivot a non-negative variable down until it reaches the value zero
 * and then pivot the variable into a column position.
 */
static int to_col(struct isl_tab *tab, struct isl_tab_var *var)
{
	int i;
	int row, col;

	if (!var->is_row)
		return;

	while (isl_int_is_pos(tab->mat->row[var->index][1])) {
		find_pivot(tab, var, NULL, -1, &row, &col);
		isl_assert(tab->mat->ctx, row != -1, return -1);
		pivot(tab, row, col);
		if (!var->is_row)
			return;
	}

	for (i = tab->n_dead; i < tab->n_col; ++i)
		if (!isl_int_is_zero(tab->mat->row[var->index][2 + i]))
			break;

	isl_assert(tab->mat->ctx, i < tab->n_col, return -1);
	pivot(tab, var->index, i);

	return 0;
}

/* We assume Gaussian elimination has been performed on the equalities.
 * The equalities can therefore never conflict.
 * Adding the equalities is currently only really useful for a later call
 * to isl_tab_ineq_type.
 */
static struct isl_tab *add_eq(struct isl_tab *tab, isl_int *eq)
{
	int i;
	int r;

	if (!tab)
		return NULL;
	r = add_row(tab, eq);
	if (r < 0)
		goto error;

	r = tab->con[r].index;
	for (i = tab->n_dead; i < tab->n_col; ++i) {
		if (isl_int_is_zero(tab->mat->row[r][2 + i]))
			continue;
		pivot(tab, r, i);
		kill_col(tab, i);
		break;
	}
	tab->n_eq++;

	return tab;
error:
	isl_tab_free(tab);
	return NULL;
}

/* Add an equality that is known to be valid for the given tableau.
 */
struct isl_tab *isl_tab_add_valid_eq(struct isl_tab *tab, isl_int *eq)
{
	struct isl_tab_var *var;
	int i;
	int r;

	if (!tab)
		return NULL;
	r = add_row(tab, eq);
	if (r < 0)
		goto error;

	var = &tab->con[r];
	r = var->index;
	if (isl_int_is_neg(tab->mat->row[r][1]))
		isl_seq_neg(tab->mat->row[r] + 1, tab->mat->row[r] + 1,
			    1 + tab->n_col);
	var->is_nonneg = 1;
	if (to_col(tab, var) < 0)
		goto error;
	var->is_nonneg = 0;
	kill_col(tab, var->index);

	return tab;
error:
	isl_tab_free(tab);
	return NULL;
}

struct isl_tab *isl_tab_from_basic_map(struct isl_basic_map *bmap)
{
	int i;
	struct isl_tab *tab;

	if (!bmap)
		return NULL;
	tab = isl_tab_alloc(bmap->ctx,
			    isl_basic_map_total_dim(bmap) + bmap->n_ineq + 1,
			    isl_basic_map_total_dim(bmap));
	if (!tab)
		return NULL;
	tab->rational = ISL_F_ISSET(bmap, ISL_BASIC_MAP_RATIONAL);
	if (ISL_F_ISSET(bmap, ISL_BASIC_MAP_EMPTY)) {
		mark_empty(tab);
		return tab;
	}
	for (i = 0; i < bmap->n_eq; ++i) {
		tab = add_eq(tab, bmap->eq[i]);
		if (!tab)
			return tab;
	}
	for (i = 0; i < bmap->n_ineq; ++i) {
		tab = isl_tab_add_ineq(tab, bmap->ineq[i]);
		if (!tab || tab->empty)
			return tab;
	}
	return tab;
}

struct isl_tab *isl_tab_from_basic_set(struct isl_basic_set *bset)
{
	return isl_tab_from_basic_map((struct isl_basic_map *)bset);
}

/* Construct a tableau corresponding to the recession cone of "bmap".
 */
struct isl_tab *isl_tab_from_recession_cone(struct isl_basic_map *bmap)
{
	isl_int cst;
	int i;
	struct isl_tab *tab;

	if (!bmap)
		return NULL;
	tab = isl_tab_alloc(bmap->ctx, bmap->n_eq + bmap->n_ineq,
				isl_basic_map_total_dim(bmap));
	if (!tab)
		return NULL;
	tab->rational = ISL_F_ISSET(bmap, ISL_BASIC_MAP_RATIONAL);

	isl_int_init(cst);
	for (i = 0; i < bmap->n_eq; ++i) {
		isl_int_swap(bmap->eq[i][0], cst);
		tab = add_eq(tab, bmap->eq[i]);
		isl_int_swap(bmap->eq[i][0], cst);
		if (!tab)
			goto done;
	}
	for (i = 0; i < bmap->n_ineq; ++i) {
		int r;
		isl_int_swap(bmap->ineq[i][0], cst);
		r = add_row(tab, bmap->ineq[i]);
		isl_int_swap(bmap->ineq[i][0], cst);
		if (r < 0)
			goto error;
		tab->con[r].is_nonneg = 1;
		push(tab, isl_tab_undo_nonneg, &tab->con[r]);
	}
done:
	isl_int_clear(cst);
	return tab;
error:
	isl_int_clear(cst);
	isl_tab_free(tab);
	return NULL;
}

/* Assuming "tab" is the tableau of a cone, check if the cone is
 * bounded, i.e., if it is empty or only contains the origin.
 */
int isl_tab_cone_is_bounded(struct isl_tab *tab)
{
	int i;

	if (!tab)
		return -1;
	if (tab->empty)
		return 1;
	if (tab->n_dead == tab->n_col)
		return 1;

	for (;;) {
		for (i = tab->n_redundant; i < tab->n_row; ++i) {
			struct isl_tab_var *var;
			var = var_from_row(tab, i);
			if (!var->is_nonneg)
				continue;
			if (sign_of_max(tab, var) != 0)
				return 0;
			close_row(tab, var);
			break;
		}
		if (tab->n_dead == tab->n_col)
			return 1;
		if (i == tab->n_row)
			return 0;
	}
}

int isl_tab_sample_is_integer(struct isl_tab *tab)
{
	int i;

	if (!tab)
		return -1;

	for (i = 0; i < tab->n_var; ++i) {
		int row;
		if (!tab->var[i].is_row)
			continue;
		row = tab->var[i].index;
		if (!isl_int_is_divisible_by(tab->mat->row[row][1],
						tab->mat->row[row][0]))
			return 0;
	}
	return 1;
}

static struct isl_vec *extract_integer_sample(struct isl_tab *tab)
{
	int i;
	struct isl_vec *vec;

	vec = isl_vec_alloc(tab->mat->ctx, 1 + tab->n_var);
	if (!vec)
		return NULL;

	isl_int_set_si(vec->block.data[0], 1);
	for (i = 0; i < tab->n_var; ++i) {
		if (!tab->var[i].is_row)
			isl_int_set_si(vec->block.data[1 + i], 0);
		else {
			int row = tab->var[i].index;
			isl_int_divexact(vec->block.data[1 + i],
				tab->mat->row[row][1], tab->mat->row[row][0]);
		}
	}

	return vec;
}

struct isl_vec *isl_tab_get_sample_value(struct isl_tab *tab)
{
	int i;
	struct isl_vec *vec;
	isl_int m;

	if (!tab)
		return NULL;

	vec = isl_vec_alloc(tab->mat->ctx, 1 + tab->n_var);
	if (!vec)
		return NULL;

	isl_int_init(m);

	isl_int_set_si(vec->block.data[0], 1);
	for (i = 0; i < tab->n_var; ++i) {
		int row;
		if (!tab->var[i].is_row) {
			isl_int_set_si(vec->block.data[1 + i], 0);
			continue;
		}
		row = tab->var[i].index;
		isl_int_gcd(m, vec->block.data[0], tab->mat->row[row][0]);
		isl_int_divexact(m, tab->mat->row[row][0], m);
		isl_seq_scale(vec->block.data, vec->block.data, m, 1 + i);
		isl_int_divexact(m, vec->block.data[0], tab->mat->row[row][0]);
		isl_int_mul(vec->block.data[1 + i], m, tab->mat->row[row][1]);
	}
	isl_seq_normalize(vec->block.data, vec->size);

	isl_int_clear(m);
	return vec;
}

/* Update "bmap" based on the results of the tableau "tab".
 * In particular, implicit equalities are made explicit, redundant constraints
 * are removed and if the sample value happens to be integer, it is stored
 * in "bmap" (unless "bmap" already had an integer sample).
 *
 * The tableau is assumed to have been created from "bmap" using
 * isl_tab_from_basic_map.
 */
struct isl_basic_map *isl_basic_map_update_from_tab(struct isl_basic_map *bmap,
	struct isl_tab *tab)
{
	int i;
	unsigned n_eq;

	if (!bmap)
		return NULL;
	if (!tab)
		return bmap;

	n_eq = tab->n_eq;
	if (tab->empty)
		bmap = isl_basic_map_set_to_empty(bmap);
	else
		for (i = bmap->n_ineq - 1; i >= 0; --i) {
			if (isl_tab_is_equality(tab, n_eq + i))
				isl_basic_map_inequality_to_equality(bmap, i);
			else if (isl_tab_is_redundant(tab, n_eq + i))
				isl_basic_map_drop_inequality(bmap, i);
		}
	if (!tab->rational &&
	    !bmap->sample && isl_tab_sample_is_integer(tab))
		bmap->sample = extract_integer_sample(tab);
	return bmap;
}

struct isl_basic_set *isl_basic_set_update_from_tab(struct isl_basic_set *bset,
	struct isl_tab *tab)
{
	return (struct isl_basic_set *)isl_basic_map_update_from_tab(
		(struct isl_basic_map *)bset, tab);
}

/* Given a non-negative variable "var", add a new non-negative variable
 * that is the opposite of "var", ensuring that var can only attain the
 * value zero.
 * If var = n/d is a row variable, then the new variable = -n/d.
 * If var is a column variables, then the new variable = -var.
 * If the new variable cannot attain non-negative values, then
 * the resulting tableau is empty.
 * Otherwise, we know the value will be zero and we close the row.
 */
static struct isl_tab *cut_to_hyperplane(struct isl_tab *tab,
	struct isl_tab_var *var)
{
	unsigned r;
	isl_int *row;
	int sgn;

	if (extend_cons(tab, 1) < 0)
		goto error;

	r = tab->n_con;
	tab->con[r].index = tab->n_row;
	tab->con[r].is_row = 1;
	tab->con[r].is_nonneg = 0;
	tab->con[r].is_zero = 0;
	tab->con[r].is_redundant = 0;
	tab->con[r].frozen = 0;
	tab->row_var[tab->n_row] = ~r;
	row = tab->mat->row[tab->n_row];

	if (var->is_row) {
		isl_int_set(row[0], tab->mat->row[var->index][0]);
		isl_seq_neg(row + 1,
			    tab->mat->row[var->index] + 1, 1 + tab->n_col);
	} else {
		isl_int_set_si(row[0], 1);
		isl_seq_clr(row + 1, 1 + tab->n_col);
		isl_int_set_si(row[2 + var->index], -1);
	}

	tab->n_row++;
	tab->n_con++;
	push(tab, isl_tab_undo_allocate, &tab->con[r]);

	sgn = sign_of_max(tab, &tab->con[r]);
	if (sgn < 0)
		mark_empty(tab);
	else {
		tab->con[r].is_nonneg = 1;
		push(tab, isl_tab_undo_nonneg, &tab->con[r]);
		/* sgn == 0 */
		close_row(tab, &tab->con[r]);
	}

	return tab;
error:
	isl_tab_free(tab);
	return NULL;
}

/* Given a tableau "tab" and an inequality constraint "con" of the tableau,
 * relax the inequality by one.  That is, the inequality r >= 0 is replaced
 * by r' = r + 1 >= 0.
 * If r is a row variable, we simply increase the constant term by one
 * (taking into account the denominator).
 * If r is a column variable, then we need to modify each row that
 * refers to r = r' - 1 by substituting this equality, effectively
 * subtracting the coefficient of the column from the constant.
 */
struct isl_tab *isl_tab_relax(struct isl_tab *tab, int con)
{
	struct isl_tab_var *var;
	if (!tab)
		return NULL;

	var = &tab->con[con];

	if (!var->is_row && !max_is_manifestly_unbounded(tab, var))
		to_row(tab, var, 1);

	if (var->is_row)
		isl_int_add(tab->mat->row[var->index][1],
		    tab->mat->row[var->index][1], tab->mat->row[var->index][0]);
	else {
		int i;

		for (i = 0; i < tab->n_row; ++i) {
			if (isl_int_is_zero(tab->mat->row[i][2 + var->index]))
				continue;
			isl_int_sub(tab->mat->row[i][1], tab->mat->row[i][1],
			    tab->mat->row[i][2 + var->index]);
		}

	}

	push(tab, isl_tab_undo_relax, var);

	return tab;
}

struct isl_tab *isl_tab_select_facet(struct isl_tab *tab, int con)
{
	if (!tab)
		return NULL;

	return cut_to_hyperplane(tab, &tab->con[con]);
}

static int may_be_equality(struct isl_tab *tab, int row)
{
	return (tab->rational ? isl_int_is_zero(tab->mat->row[row][1])
			      : isl_int_lt(tab->mat->row[row][1],
					    tab->mat->row[row][0])) &&
		isl_seq_first_non_zero(tab->mat->row[row] + 2 + tab->n_dead,
					tab->n_col - tab->n_dead) != -1;
}

/* Check for (near) equalities among the constraints.
 * A constraint is an equality if it is non-negative and if
 * its maximal value is either
 *	- zero (in case of rational tableaus), or
 *	- strictly less than 1 (in case of integer tableaus)
 *
 * We first mark all non-redundant and non-dead variables that
 * are not frozen and not obviously not an equality.
 * Then we iterate over all marked variables if they can attain
 * any values larger than zero or at least one.
 * If the maximal value is zero, we mark any column variables
 * that appear in the row as being zero and mark the row as being redundant.
 * Otherwise, if the maximal value is strictly less than one (and the
 * tableau is integer), then we restrict the value to being zero
 * by adding an opposite non-negative variable.
 */
struct isl_tab *isl_tab_detect_equalities(struct isl_tab *tab)
{
	int i;
	unsigned n_marked;

	if (!tab)
		return NULL;
	if (tab->empty)
		return tab;
	if (tab->n_dead == tab->n_col)
		return tab;

	n_marked = 0;
	for (i = tab->n_redundant; i < tab->n_row; ++i) {
		struct isl_tab_var *var = var_from_row(tab, i);
		var->marked = !var->frozen && var->is_nonneg &&
			may_be_equality(tab, i);
		if (var->marked)
			n_marked++;
	}
	for (i = tab->n_dead; i < tab->n_col; ++i) {
		struct isl_tab_var *var = var_from_col(tab, i);
		var->marked = !var->frozen && var->is_nonneg;
		if (var->marked)
			n_marked++;
	}
	while (n_marked) {
		struct isl_tab_var *var;
		for (i = tab->n_redundant; i < tab->n_row; ++i) {
			var = var_from_row(tab, i);
			if (var->marked)
				break;
		}
		if (i == tab->n_row) {
			for (i = tab->n_dead; i < tab->n_col; ++i) {
				var = var_from_col(tab, i);
				if (var->marked)
					break;
			}
			if (i == tab->n_col)
				break;
		}
		var->marked = 0;
		n_marked--;
		if (sign_of_max(tab, var) == 0)
			close_row(tab, var);
		else if (!tab->rational && !at_least_one(tab, var)) {
			tab = cut_to_hyperplane(tab, var);
			return isl_tab_detect_equalities(tab);
		}
		for (i = tab->n_redundant; i < tab->n_row; ++i) {
			var = var_from_row(tab, i);
			if (!var->marked)
				continue;
			if (may_be_equality(tab, i))
				continue;
			var->marked = 0;
			n_marked--;
		}
	}

	return tab;
}

/* Check for (near) redundant constraints.
 * A constraint is redundant if it is non-negative and if
 * its minimal value (temporarily ignoring the non-negativity) is either
 *	- zero (in case of rational tableaus), or
 *	- strictly larger than -1 (in case of integer tableaus)
 *
 * We first mark all non-redundant and non-dead variables that
 * are not frozen and not obviously negatively unbounded.
 * Then we iterate over all marked variables if they can attain
 * any values smaller than zero or at most negative one.
 * If not, we mark the row as being redundant (assuming it hasn't
 * been detected as being obviously redundant in the mean time).
 */
struct isl_tab *isl_tab_detect_redundant(struct isl_tab *tab)
{
	int i;
	unsigned n_marked;

	if (!tab)
		return NULL;
	if (tab->empty)
		return tab;
	if (tab->n_redundant == tab->n_row)
		return tab;

	n_marked = 0;
	for (i = tab->n_redundant; i < tab->n_row; ++i) {
		struct isl_tab_var *var = var_from_row(tab, i);
		var->marked = !var->frozen && var->is_nonneg;
		if (var->marked)
			n_marked++;
	}
	for (i = tab->n_dead; i < tab->n_col; ++i) {
		struct isl_tab_var *var = var_from_col(tab, i);
		var->marked = !var->frozen && var->is_nonneg &&
			!min_is_manifestly_unbounded(tab, var);
		if (var->marked)
			n_marked++;
	}
	while (n_marked) {
		struct isl_tab_var *var;
		for (i = tab->n_redundant; i < tab->n_row; ++i) {
			var = var_from_row(tab, i);
			if (var->marked)
				break;
		}
		if (i == tab->n_row) {
			for (i = tab->n_dead; i < tab->n_col; ++i) {
				var = var_from_col(tab, i);
				if (var->marked)
					break;
			}
			if (i == tab->n_col)
				break;
		}
		var->marked = 0;
		n_marked--;
		if ((tab->rational ? (sign_of_min(tab, var) >= 0)
				   : !min_at_most_neg_one(tab, var)) &&
		    !var->is_redundant)
			mark_redundant(tab, var->index);
		for (i = tab->n_dead; i < tab->n_col; ++i) {
			var = var_from_col(tab, i);
			if (!var->marked)
				continue;
			if (!min_is_manifestly_unbounded(tab, var))
				continue;
			var->marked = 0;
			n_marked--;
		}
	}

	return tab;
}

int isl_tab_is_equality(struct isl_tab *tab, int con)
{
	int row;

	if (!tab)
		return -1;
	if (tab->con[con].is_zero)
		return 1;
	if (tab->con[con].is_redundant)
		return 0;
	if (!tab->con[con].is_row)
		return tab->con[con].index < tab->n_dead;

	row = tab->con[con].index;

	return isl_int_is_zero(tab->mat->row[row][1]) &&
		isl_seq_first_non_zero(tab->mat->row[row] + 2 + tab->n_dead,
					tab->n_col - tab->n_dead) == -1;
}

/* Return the minimial value of the affine expression "f" with denominator
 * "denom" in *opt, *opt_denom, assuming the tableau is not empty and
 * the expression cannot attain arbitrarily small values.
 * If opt_denom is NULL, then *opt is rounded up to the nearest integer.
 * The return value reflects the nature of the result (empty, unbounded,
 * minmimal value returned in *opt).
 */
enum isl_lp_result isl_tab_min(struct isl_tab *tab,
	isl_int *f, isl_int denom, isl_int *opt, isl_int *opt_denom,
	unsigned flags)
{
	int r;
	enum isl_lp_result res = isl_lp_ok;
	struct isl_tab_var *var;
	struct isl_tab_undo *snap;

	if (tab->empty)
		return isl_lp_empty;

	snap = isl_tab_snap(tab);
	r = add_row(tab, f);
	if (r < 0)
		return isl_lp_error;
	var = &tab->con[r];
	isl_int_mul(tab->mat->row[var->index][0],
		    tab->mat->row[var->index][0], denom);
	for (;;) {
		int row, col;
		find_pivot(tab, var, var, -1, &row, &col);
		if (row == var->index) {
			res = isl_lp_unbounded;
			break;
		}
		if (row == -1)
			break;
		pivot(tab, row, col);
	}
	if (isl_tab_rollback(tab, snap) < 0)
		return isl_lp_error;
	if (ISL_FL_ISSET(flags, ISL_TAB_SAVE_DUAL)) {
		int i;

		isl_vec_free(tab->dual);
		tab->dual = isl_vec_alloc(tab->mat->ctx, 1 + tab->n_con);
		if (!tab->dual)
			return isl_lp_error;
		isl_int_set(tab->dual->el[0], tab->mat->row[var->index][0]);
		for (i = 0; i < tab->n_con; ++i) {
			if (tab->con[i].is_row)
				isl_int_set_si(tab->dual->el[1 + i], 0);
			else {
				int pos = 2 + tab->con[i].index;
				isl_int_set(tab->dual->el[1 + i],
					    tab->mat->row[var->index][pos]);
			}
		}
	}
	if (res == isl_lp_ok) {
		if (opt_denom) {
			isl_int_set(*opt, tab->mat->row[var->index][1]);
			isl_int_set(*opt_denom, tab->mat->row[var->index][0]);
		} else
			isl_int_cdiv_q(*opt, tab->mat->row[var->index][1],
					     tab->mat->row[var->index][0]);
	}
	return res;
}

int isl_tab_is_redundant(struct isl_tab *tab, int con)
{
	int row;
	unsigned n_col;

	if (!tab)
		return -1;
	if (tab->con[con].is_zero)
		return 0;
	if (tab->con[con].is_redundant)
		return 1;
	return tab->con[con].is_row && tab->con[con].index < tab->n_redundant;
}

/* Take a snapshot of the tableau that can be restored by s call to
 * isl_tab_rollback.
 */
struct isl_tab_undo *isl_tab_snap(struct isl_tab *tab)
{
	if (!tab)
		return NULL;
	tab->need_undo = 1;
	return tab->top;
}

/* Undo the operation performed by isl_tab_relax.
 */
static void unrelax(struct isl_tab *tab, struct isl_tab_var *var)
{
	if (!var->is_row && !max_is_manifestly_unbounded(tab, var))
		to_row(tab, var, 1);

	if (var->is_row)
		isl_int_sub(tab->mat->row[var->index][1],
		    tab->mat->row[var->index][1], tab->mat->row[var->index][0]);
	else {
		int i;

		for (i = 0; i < tab->n_row; ++i) {
			if (isl_int_is_zero(tab->mat->row[i][2 + var->index]))
				continue;
			isl_int_add(tab->mat->row[i][1], tab->mat->row[i][1],
			    tab->mat->row[i][2 + var->index]);
		}

	}
}

static void perform_undo(struct isl_tab *tab, struct isl_tab_undo *undo)
{
	switch(undo->type) {
	case isl_tab_undo_empty:
		tab->empty = 0;
		break;
	case isl_tab_undo_nonneg:
		undo->var->is_nonneg = 0;
		break;
	case isl_tab_undo_redundant:
		undo->var->is_redundant = 0;
		tab->n_redundant--;
		break;
	case isl_tab_undo_zero:
		undo->var->is_zero = 0;
		tab->n_dead--;
		break;
	case isl_tab_undo_allocate:
		if (!undo->var->is_row) {
			if (max_is_manifestly_unbounded(tab, undo->var))
				to_row(tab, undo->var, -1);
			else
				to_row(tab, undo->var, 1);
		}
		drop_row(tab, undo->var->index);
		break;
	case isl_tab_undo_relax:
		unrelax(tab, undo->var);
		break;
	}
}

/* Return the tableau to the state it was in when the snapshot "snap"
 * was taken.
 */
int isl_tab_rollback(struct isl_tab *tab, struct isl_tab_undo *snap)
{
	struct isl_tab_undo *undo, *next;

	if (!tab)
		return -1;

	tab->in_undo = 1;
	for (undo = tab->top; undo && undo != &tab->bottom; undo = next) {
		next = undo->next;
		if (undo == snap)
			break;
		perform_undo(tab, undo);
		free(undo);
	}
	tab->in_undo = 0;
	tab->top = undo;
	if (!undo)
		return -1;
	return 0;
}

/* The given row "row" represents an inequality violated by all
 * points in the tableau.  Check for some special cases of such
 * separating constraints.
 * In particular, if the row has been reduced to the constant -1,
 * then we know the inequality is adjacent (but opposite) to
 * an equality in the tableau.
 * If the row has been reduced to r = -1 -r', with r' an inequality
 * of the tableau, then the inequality is adjacent (but opposite)
 * to the inequality r'.
 */
static enum isl_ineq_type separation_type(struct isl_tab *tab, unsigned row)
{
	int pos;

	if (tab->rational)
		return isl_ineq_separate;

	if (!isl_int_is_one(tab->mat->row[row][0]))
		return isl_ineq_separate;
	if (!isl_int_is_negone(tab->mat->row[row][1]))
		return isl_ineq_separate;

	pos = isl_seq_first_non_zero(tab->mat->row[row] + 2 + tab->n_dead,
					tab->n_col - tab->n_dead);
	if (pos == -1)
		return isl_ineq_adj_eq;

	if (!isl_int_is_negone(tab->mat->row[row][2 + tab->n_dead + pos]))
		return isl_ineq_separate;

	pos = isl_seq_first_non_zero(
			tab->mat->row[row] + 2 + tab->n_dead + pos + 1,
			tab->n_col - tab->n_dead - pos - 1);

	return pos == -1 ? isl_ineq_adj_ineq : isl_ineq_separate;
}

/* Check the effect of inequality "ineq" on the tableau "tab".
 * The result may be
 *	isl_ineq_redundant:	satisfied by all points in the tableau
 *	isl_ineq_separate:	satisfied by no point in the tableau
 *	isl_ineq_cut:		satisfied by some by not all points
 *	isl_ineq_adj_eq:	adjacent to an equality
 *	isl_ineq_adj_ineq:	adjacent to an inequality.
 */
enum isl_ineq_type isl_tab_ineq_type(struct isl_tab *tab, isl_int *ineq)
{
	enum isl_ineq_type type = isl_ineq_error;
	struct isl_tab_undo *snap = NULL;
	int con;
	int row;

	if (!tab)
		return isl_ineq_error;

	if (extend_cons(tab, 1) < 0)
		return isl_ineq_error;

	snap = isl_tab_snap(tab);

	con = add_row(tab, ineq);
	if (con < 0)
		goto error;

	row = tab->con[con].index;
	if (is_redundant(tab, row))
		type = isl_ineq_redundant;
	else if (isl_int_is_neg(tab->mat->row[row][1]) &&
		 (tab->rational ||
		    isl_int_abs_ge(tab->mat->row[row][1],
				   tab->mat->row[row][0]))) {
		if (at_least_zero(tab, &tab->con[con]))
			type = isl_ineq_cut;
		else
			type = separation_type(tab, row);
	} else if (tab->rational ? (sign_of_min(tab, &tab->con[con]) < 0)
			     : min_at_most_neg_one(tab, &tab->con[con]))
		type = isl_ineq_cut;
	else
		type = isl_ineq_redundant;

	if (isl_tab_rollback(tab, snap))
		return isl_ineq_error;
	return type;
error:
	isl_tab_rollback(tab, snap);
	return isl_ineq_error;
}

void isl_tab_dump(struct isl_tab *tab, FILE *out, int indent)
{
	unsigned r, c;
	int i;

	if (!tab) {
		fprintf(out, "%*snull tab\n", indent, "");
		return;
	}
	fprintf(out, "%*sn_redundant: %d, n_dead: %d", indent, "",
		tab->n_redundant, tab->n_dead);
	if (tab->rational)
		fprintf(out, ", rational");
	if (tab->empty)
		fprintf(out, ", empty");
	fprintf(out, "\n");
	fprintf(out, "%*s[", indent, "");
	for (i = 0; i < tab->n_var; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "%c%d%s", tab->var[i].is_row ? 'r' : 'c',
					tab->var[i].index,
					tab->var[i].is_zero ? " [=0]" :
					tab->var[i].is_redundant ? " [R]" : "");
	}
	fprintf(out, "]\n");
	fprintf(out, "%*s[", indent, "");
	for (i = 0; i < tab->n_con; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "%c%d%s", tab->con[i].is_row ? 'r' : 'c',
					tab->con[i].index,
					tab->con[i].is_zero ? " [=0]" :
					tab->con[i].is_redundant ? " [R]" : "");
	}
	fprintf(out, "]\n");
	fprintf(out, "%*s[", indent, "");
	for (i = 0; i < tab->n_row; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "r%d: %d%s", i, tab->row_var[i],
		    var_from_row(tab, i)->is_nonneg ? " [>=0]" : "");
	}
	fprintf(out, "]\n");
	fprintf(out, "%*s[", indent, "");
	for (i = 0; i < tab->n_col; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "c%d: %d%s", i, tab->col_var[i],
		    var_from_col(tab, i)->is_nonneg ? " [>=0]" : "");
	}
	fprintf(out, "]\n");
	r = tab->mat->n_row;
	tab->mat->n_row = tab->n_row;
	c = tab->mat->n_col;
	tab->mat->n_col = 2 + tab->n_col;
	isl_mat_dump(tab->mat, out, indent);
	tab->mat->n_row = r;
	tab->mat->n_col = c;
}
