#ifndef ISL_TAB_H
#define ISL_TAB_H

#include "isl_lp.h"
#include "isl_map.h"
#include "isl_mat.h"

struct isl_tab_var {
	int index;
	unsigned is_row : 1;
	unsigned is_nonneg : 1;
	unsigned is_zero : 1;
	unsigned is_redundant : 1;
	unsigned marked : 1;
	unsigned frozen : 1;
};

enum isl_tab_undo_type {
	isl_tab_undo_bottom,
	isl_tab_undo_empty,
	isl_tab_undo_nonneg,
	isl_tab_undo_redundant,
	isl_tab_undo_zero,
	isl_tab_undo_allocate,
	isl_tab_undo_relax,
};

struct isl_tab_undo {
	enum isl_tab_undo_type	type;
	struct isl_tab_var	*var;
	struct isl_tab_undo	*next;
};

/* The tableau maintains equality relations.
 * Each column and each row is associated to a variable or a constraint.
 * The "value" of an inequality constraint is the value of the corresponding
 * slack variable.
 * The "row_var" and "col_var" arrays map column and row indices
 * to indices in the "var" and "con" arrays.  The elements of these
 * arrays maintain extra information about the variables and the constraints.
 * Each row expresses the corresponding row variable as an affine expression
 * of the column variables.
 * The first two columns in the matrix contain the common denominator of
 * the row and the numerator of the constant term.  The third column
 * in the matrix is called column 0 with respect to the col_var array.
 * The sample value of the tableau is the value that assigns zero
 * to all the column variables and the constant term of each affine
 * expression to the corresponding row variable.
 * The operations on the tableau maintain the property that the sample
 * value satisfies the non-negativity constraints (usually on the slack
 * variables).
 *
 * The first n_dead column variables have their values fixed to zero.
 * The corresponding tab_vars are flagged "is_zero".
 * Some of the rows that have have zero coefficients in all but
 * the dead columns are also flagged "is_zero".
 *
 * The first n_redundant rows correspond to inequality constraints
 * that are always satisfied for any value satisfying the non-redundant
 * rows.  The corresponding tab_vars are flagged "is_redundant".
 * A row variable that is flagged "is_zero" is also flagged "is_redundant"
 * since the constraint has been reduced to 0 = 0 and is therefore always
 * satisfied.
 *
 * Dead columns and redundant rows are detected on the fly.
 * However, the basic operations do not ensure that all dead columns
 * or all redundant rows are detected.
 * isl_tab_detect_equalities and isl_tab_detect_redundant can be used
 * to peform and exhaustive search for dead columns and redundant rows.
 */
struct isl_tab {
	struct isl_mat *mat;

	unsigned n_row;
	unsigned n_col;
	unsigned n_dead;
	unsigned n_redundant;

	unsigned n_var;
	unsigned n_con;
	unsigned n_eq;
	unsigned max_con;
	struct isl_tab_var *var;
	struct isl_tab_var *con;
	int *row_var;	/* v >= 0 -> var v;	v < 0 -> con ~v */
	int *col_var;	/* v >= 0 -> var v;	v < 0 -> con ~v */

	struct isl_tab_undo bottom;
	struct isl_tab_undo *top;

	struct isl_vec *dual;

	unsigned need_undo : 1;
	unsigned rational : 1;
	unsigned empty : 1;
	unsigned in_undo : 1;
};

struct isl_tab *isl_tab_alloc(struct isl_ctx *ctx,
	unsigned n_row, unsigned n_var);
void isl_tab_free(struct isl_tab *tab);

struct isl_tab *isl_tab_from_basic_map(struct isl_basic_map *bmap);
struct isl_tab *isl_tab_from_basic_set(struct isl_basic_set *bset);
struct isl_tab *isl_tab_from_recession_cone(struct isl_basic_map *bmap);
int isl_tab_cone_is_bounded(struct isl_tab *tab);
struct isl_basic_map *isl_basic_map_update_from_tab(struct isl_basic_map *bmap,
	struct isl_tab *tab);
struct isl_basic_set *isl_basic_set_update_from_tab(struct isl_basic_set *bset,
	struct isl_tab *tab);
struct isl_tab *isl_tab_detect_equalities(struct isl_tab *tab);
struct isl_tab *isl_tab_detect_redundant(struct isl_tab *tab);
#define ISL_TAB_SAVE_DUAL	(1 << 0)
enum isl_lp_result isl_tab_min(struct isl_tab *tab,
	isl_int *f, isl_int denom, isl_int *opt, isl_int *opt_denom,
	unsigned flags);

struct isl_tab *isl_tab_extend(struct isl_tab *tab, unsigned n_new);
struct isl_tab *isl_tab_add_ineq(struct isl_tab *tab, isl_int *ineq);
struct isl_tab *isl_tab_add_valid_eq(struct isl_tab *tab, isl_int *eq);

int isl_tab_is_equality(struct isl_tab *tab, int con);
int isl_tab_is_redundant(struct isl_tab *tab, int con);

int isl_tab_sample_is_integer(struct isl_tab *tab);
struct isl_vec *isl_tab_get_sample_value(struct isl_tab *tab);

enum isl_ineq_type {
	isl_ineq_error = -1,
	isl_ineq_redundant,
	isl_ineq_separate,
	isl_ineq_cut,
	isl_ineq_adj_eq,
	isl_ineq_adj_ineq,
};

enum isl_ineq_type isl_tab_ineq_type(struct isl_tab *tab, isl_int *ineq);

struct isl_tab_undo *isl_tab_snap(struct isl_tab *tab);
int isl_tab_rollback(struct isl_tab *tab, struct isl_tab_undo *snap);

struct isl_tab *isl_tab_relax(struct isl_tab *tab, int con);
struct isl_tab *isl_tab_select_facet(struct isl_tab *tab, int con);

void isl_tab_dump(struct isl_tab *tab, FILE *out, int indent);

#endif
