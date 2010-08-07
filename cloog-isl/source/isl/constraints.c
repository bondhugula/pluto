#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <cloog/isl/cloog.h>
#include <cloog/isl/backend.h>
#include <isl_set.h>


#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n)*sizeof(type))


/******************************************************************************
 *                             Memory leaks hunting                           *
 ******************************************************************************/



void cloog_constraint_set_free(CloogConstraintSet *constraints)
{
	isl_basic_set_free(constraints);
}


int cloog_constraint_set_contains_level(CloogConstraintSet *constraints,
			int level, int nb_parameters)
{
	return isl_basic_set_n_dim(constraints) >= level;
}

struct cloog_isl_dim {
	enum isl_dim_type type;
	int		  pos;
};

static struct cloog_isl_dim set_cloog_dim_to_isl_dim(
					CloogConstraintSet *bset, int pos)
{
	enum isl_dim_type types[] = { isl_dim_set, isl_dim_div, isl_dim_param };
	int i;
	struct cloog_isl_dim ci_dim;

	for (i = 0; i < 3; ++i) {
		unsigned dim = isl_basic_set_dim(bset, types[i]);
		if (pos < dim) {
			ci_dim.type = types[i];
			ci_dim.pos = pos;
			return ci_dim;
		}
		pos -= dim;
	}
	assert(0);
}

/* Check if the variable at position level is defined by an
 * equality.  If so, return the row number.  Otherwise, return -1.
 */
CloogConstraint cloog_constraint_set_defining_equality(
	CloogConstraintSet *bset, int level)
{
	struct isl_constraint *c;
	struct cloog_isl_dim dim;

	dim = set_cloog_dim_to_isl_dim(bset, level - 1);
	if (isl_basic_set_has_defining_equality(bset, dim.type, dim.pos, &c))
		return c;
	else
		return NULL;
}


/* Check if the variable (e) at position level is defined by a
 * pair of inequalities
 *		 <a, i> + -m e +  <b, p> + k1 >= 0
 *		<-a, i> +  m e + <-b, p> + k2 >= 0
 * with 0 <= k1 + k2 < m
 * If so return the row number of the upper bound and set *lower
 * to the row number of the lower bound.  If not, return -1.
 *
 * If the variable at position level occurs in any other constraint,
 * then we currently return -1.  The modulo guard that we would generate
 * would still be correct, but we would also need to generate
 * guards corresponding to the other constraints, and this has not
 * been implemented yet.
 */
CloogConstraint cloog_constraint_set_defining_inequalities(
	CloogConstraintSet *bset,
	int level, CloogConstraint *lower, int nb_par)
{
	struct isl_constraint *upper;
	struct isl_constraint *c;
	struct cloog_isl_dim dim;

	dim = set_cloog_dim_to_isl_dim(bset, level - 1);
	if (!isl_basic_set_has_defining_inequalities(bset, dim.type, dim.pos,
								lower, &upper))
		return cloog_constraint_invalid();
	for (c = isl_basic_set_first_constraint(isl_basic_set_copy(bset)); c;
	     c = isl_constraint_next(c)) {
		if (isl_constraint_is_equal(c, *lower))
			continue;
		if (isl_constraint_is_equal(c, upper))
			continue;
		if (cloog_constraint_involves(c, level-1)) {
			isl_constraint_free(*lower);
			isl_constraint_free(upper);
			*lower = NULL;
			isl_constraint_free(c);
			return NULL;
		}
	}
	return upper;
}

int cloog_constraint_set_total_dimension(CloogConstraintSet *constraints)
{
	return isl_basic_set_total_dim(constraints);
}

int cloog_constraint_set_n_iterators(CloogConstraintSet *constraints, int n_par)
{
	return isl_basic_set_n_dim(constraints);
}


/******************************************************************************
 *                        Equalities spreading functions                      *
 ******************************************************************************/


/* Equalities are stored inside a CloogMatrix data structure called "equal".
 * This matrix has (nb_scattering + nb_iterators + 1) rows (i.e. total
 * dimensions + 1, the "+ 1" is because a statement can be included inside an
 * external loop without iteration domain), and (nb_scattering + nb_iterators +
 * nb_parameters + 2) columns (all unknowns plus the scalar plus the equality
 * type). The ith row corresponds to the equality "= 0" for the ith dimension
 * iterator. The first column gives the equality type (0: no equality, then
 * EQTYPE_* -see pprint.h-). At each recursion of pprint, if an equality for
 * the current level is found, the corresponding row is updated. Then the
 * equality if it exists is used to simplify expressions (e.g. if we have 
 * "i+1" while we know that "i=2", we simplify it in "3"). At the end of
 * the pprint call, the corresponding row is reset to zero.
 */

CloogEqualities *cloog_equal_alloc(int n, int nb_levels, int nb_parameters)
{
	int i;
	CloogEqualities *equal = ALLOC(CloogEqualities);

	equal->total_dim = nb_levels - 1 + nb_parameters;
	equal->n = n;
	equal->constraints = ALLOCN(struct isl_basic_set *, n);
	equal->types = ALLOCN(int, n);
	for (i = 0; i < n; ++i) {
		equal->constraints[i] = NULL;
		equal->types[i] = EQTYPE_NONE;
	}
	return equal;
}

int cloog_equal_total_dimension(CloogEqualities *equal)
{
	return equal->total_dim;
}

void cloog_equal_free(CloogEqualities *equal)
{
	int i;

	for (i = 0; i < equal->n; ++i)
		isl_basic_set_free(equal->constraints[i]);
	free(equal->constraints);
	free(equal->types);
	free(equal);
}

int cloog_equal_count(CloogEqualities *equal)
{
	return equal->n;
}


/**
 * cloog_constraint_equal_type function :
 * This function returns the type of the equality in the constraint (line) of
 * (constraints) for the element (level). An equality is 'constant' iff all
 * other factors are null except the constant one. It is a 'pure item' iff
 * it is equal or opposite to a single variable or parameter.
 * Otherwise it is an 'affine expression'.
 * For instance:
 *   i = -13 is constant, i = j, j = -M are pure items,
 *   j = 2*M, i = j+1, 2*j = M are affine expressions.
 *
 * - constraints is the matrix of constraints,
 * - level is the column number in equal of the element which is 'equal to',
 */
static int cloog_constraint_equal_type(CloogConstraint constraint, int level)
{ 
	int i;
	isl_int c;
	int type = EQTYPE_NONE;
    
	isl_int_init(c);
	isl_constraint_get_constant(constraint, &c);
	if (!isl_int_is_zero(c))
		type = EQTYPE_CONSTANT;
	isl_constraint_get_coefficient(constraint, isl_dim_set, level - 1, &c);
	if (!isl_int_is_one(c) && !isl_int_is_negone(c))
		type = EQTYPE_EXAFFINE;
	for (i = 0; i < isl_constraint_dim(constraint, isl_dim_param); ++i) {
		isl_constraint_get_coefficient(constraint, isl_dim_param, i, &c);
		if (isl_int_is_zero(c))
			continue;
		if ((!isl_int_is_one(c) && !isl_int_is_negone(c)) ||
		    type != EQTYPE_NONE) {
			type = EQTYPE_EXAFFINE;
			break;
		}
		type = EQTYPE_PUREITEM;
	}
	for (i = 0; i < isl_constraint_dim(constraint, isl_dim_set); ++i) {
		if (i == level - 1)
			continue;
		isl_constraint_get_coefficient(constraint, isl_dim_set, i, &c);
		if (isl_int_is_zero(c))
			continue;
		if ((!isl_int_is_one(c) && !isl_int_is_negone(c)) ||
		    type != EQTYPE_NONE) {
			type = EQTYPE_EXAFFINE;
			break;
		}
		type = EQTYPE_PUREITEM;
	}
	for (i = 0; i < isl_constraint_dim(constraint, isl_dim_div); ++i) {
		isl_constraint_get_coefficient(constraint, isl_dim_div, i, &c);
		if (isl_int_is_zero(c))
			continue;
		if ((!isl_int_is_one(c) && !isl_int_is_negone(c)) ||
		    type != EQTYPE_NONE) {
			type = EQTYPE_EXAFFINE;
			break;
		}
		type = EQTYPE_PUREITEM;
	}
	isl_int_clear(c);

	if (type == EQTYPE_NONE)
		type = EQTYPE_CONSTANT;

	return type;
}


int cloog_equal_type(CloogEqualities *equal, int level)
{
	return equal->types[level-1];
}


/**
 * cloog_equal_add function:
 * This function updates the row (level-1) of the equality matrix (equal) with
 * the row that corresponds to the row (line) of the matrix (matrix).
 * - equal is the matrix of equalities,
 * - matrix is the matrix of constraints,
 * - level is the column number in matrix of the element which is 'equal to',
 * - line is the line number in matrix of the constraint we want to study,
 * - the infos structure gives the user all options on code printing and more.
 **
 * line is set to and invalid constraint for equalities that CLooG itself has
 * discovered because the lower and upper bound of a loop happened to be equal.
 * This situation shouldn't happen in the isl port since isl should
 * have found the equality itself.
 */
void cloog_equal_add(CloogEqualities *equal, CloogConstraintSet *matrix,
			int level, CloogConstraint line, int nb_par)
{ 
	struct isl_basic_set *bset;
	unsigned nparam;
	assert(cloog_constraint_is_valid(line));
  
	equal->types[level-1] = cloog_constraint_equal_type(line, level);
	bset = isl_basic_set_from_constraint(isl_constraint_copy(line));
	nparam = isl_basic_set_n_param(bset);
	bset = isl_basic_set_extend(bset, nparam,
				    equal->total_dim - nparam, 0, 0, 0);
	bset = isl_basic_set_finalize(bset);
	equal->constraints[level-1] = bset;
}


/**
 * cloog_equal_del function :
 * This function reset the equality corresponding to the iterator (level)
 * in the equality matrix (equal).
 * - July 2nd 2002: first version. 
 */
void cloog_equal_del(CloogEqualities *equal, int level)
{ 
	equal->types[level-1] = EQTYPE_NONE;
	isl_basic_set_free(equal->constraints[level-1]);
	equal->constraints[level-1] = NULL;
}



/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

/**
 * Function cloog_constraint_set_normalize:
 * This function will modify the constraint system in such a way that when
 * there is an equality depending on the element at level 'level', there are
 * no more (in)equalities depending on this element.
 *
 * The simplified form of isl automatically satisfies this condition.
 */
void cloog_constraint_set_normalize(CloogConstraintSet *matrix, int level)
{
}



/**
 * cloog_constraint_set_copy function:
 * this functions builds and returns a "hard copy" (not a pointer copy) of a
 * CloogConstraintSet data structure.
 */
CloogConstraintSet *cloog_constraint_set_copy(CloogConstraintSet *bset)
{
	return isl_basic_set_dup(bset);
}


/**
 * cloog_constraint_set_simplify function:
 * this function simplify all constraints inside the matrix "matrix" thanks to
 * an equality matrix "equal" that gives for some elements of the affine
 * constraint an equality with other elements, preferably constants.
 * For instance, if a row of the matrix contains i+j+3>=0 and the equality
 * matrix gives i=n and j=2, the constraint is simplified to n+3>=0. The
 * simplified constraints are returned back inside a new simplified matrix.
 * - matrix is the set of constraints to simplify,
 * - equal is the matrix of equalities,
 * - level is a level we don't want to simplify (-1 if none),
 * - nb_par is the number of parameters of the program.
 **
 * isl should have performed these simplifications already in isl_set_gist.
 */
CloogConstraintSet *cloog_constraint_set_simplify(CloogConstraintSet *matrix,
	CloogEqualities *equal, int level, int nb_par)
{
	return cloog_constraint_set_copy(matrix);
}


static struct cloog_isl_dim constraint_cloog_dim_to_isl_dim(
					CloogConstraint constraint, int pos)
{
	enum isl_dim_type types[] = { isl_dim_set, isl_dim_div, isl_dim_param };
	int i;
	struct cloog_isl_dim ci_dim;

	for (i = 0; i < 3; ++i) {
		unsigned dim = isl_constraint_dim(constraint, types[i]);
		if (pos < dim) {
			ci_dim.type = types[i];
			ci_dim.pos = pos;
			return ci_dim;
		}
		pos -= dim;
	}
	assert(0);
}

static struct clast_expr *div_expr(CloogConstraint constraint, int pos,
					CloogNames *names)
{
	int i, nb_elts;
	unsigned dim = cloog_constraint_total_dimension(constraint);
	cloog_int_t c;
	struct clast_reduction *r;
	struct clast_expr *e = NULL;
	struct isl_div *div;

	div = isl_constraint_div(constraint, pos);

	cloog_int_init(c);
	for (i = 0, nb_elts = 0; i < dim; ++i) {
		struct cloog_isl_dim dim;

		dim = constraint_cloog_dim_to_isl_dim(constraint, i);
		isl_div_get_coefficient(div, dim.type, dim.pos, &c);
		if (!cloog_int_is_zero(c))
			++nb_elts;
	}
	isl_div_get_constant(div, &c);
	if (!cloog_int_is_zero(c))
		++nb_elts;

	r = new_clast_reduction(clast_red_sum, nb_elts);
	for (i = 0, nb_elts = 0; i < dim; ++i) {
		struct clast_expr *v;
		struct cloog_isl_dim dim;

		dim = constraint_cloog_dim_to_isl_dim(constraint, i);
		isl_div_get_coefficient(div, dim.type, dim.pos, &c);
		if (cloog_int_is_zero(c))
			continue;

		v = cloog_constraint_variable_expr(constraint, 1 + i, names);

		r->elts[nb_elts++] = &new_clast_term(c, v)->expr;
	}
	isl_div_get_constant(div, &c);
	if (!cloog_int_is_zero(c))
		r->elts[nb_elts++] = &new_clast_term(c, NULL)->expr;

	isl_div_get_denominator(div, &c);
	e = &new_clast_binary(clast_bin_fdiv, &r->expr, c)->expr;

	cloog_int_clear(c);

	isl_div_free(div);

	return e;
}

/**
 * Return clast_expr corresponding to the variable "level" (1 based) in
 * the given constraint.
 */
struct clast_expr *cloog_constraint_variable_expr(CloogConstraint constraint,
	int level, CloogNames *names)
{
	struct cloog_isl_dim dim;
	const char *name;

	assert(constraint);

	dim = constraint_cloog_dim_to_isl_dim(constraint, level - 1);
	if (dim.type == isl_dim_div)
		return div_expr(constraint, dim.pos, names);

	if (dim.type == isl_dim_set)
		name = cloog_names_name_at_level(names, level);
	else
		name = names->parameters[dim.pos];

	return &new_clast_name(name)->expr;
}


/**
 * Return true if constraint c involves variable v (zero-based).
 */
int cloog_constraint_involves(CloogConstraint constraint, int v)
{
	isl_int c;
	int res;

	isl_int_init(c);
	cloog_constraint_coefficient_get(constraint, v, &c);
	res = !isl_int_is_zero(c);
	isl_int_clear(c);
	return res;
}

int cloog_constraint_is_lower_bound(CloogConstraint constraint, int v)
{
	isl_int c;
	int res;

	isl_int_init(c);
	cloog_constraint_coefficient_get(constraint, v, &c);
	res = isl_int_is_pos(c);
	isl_int_clear(c);
	return res;
}

int cloog_constraint_is_upper_bound(CloogConstraint constraint, int v)
{
	isl_int c;
	int res;

	isl_int_init(c);
	cloog_constraint_coefficient_get(constraint, v, &c);
	res = isl_int_is_neg(c);
	isl_int_clear(c);
	return res;
}

int cloog_constraint_is_equality(CloogConstraint constraint)
{
	return isl_constraint_is_equality(constraint);
}

void cloog_constraint_clear(CloogConstraint constraint)
{
	isl_constraint_clear(constraint);
}

void cloog_constraint_coefficient_get(CloogConstraint constraint,
			int var, cloog_int_t *val)
{
	struct cloog_isl_dim dim;

	if (!constraint)
		return;

	dim = constraint_cloog_dim_to_isl_dim(constraint, var);
	isl_constraint_get_coefficient(constraint, dim.type, dim.pos, val);
}

void cloog_constraint_coefficient_set(CloogConstraint constraint,
			int var, cloog_int_t val)
{
	struct cloog_isl_dim dim;

	assert(constraint);

	dim = constraint_cloog_dim_to_isl_dim(constraint, var);
	isl_constraint_set_coefficient(constraint, dim.type, dim.pos, val);
}

void cloog_constraint_constant_get(CloogConstraint constraint, cloog_int_t *val)
{
	isl_constraint_get_constant(constraint, val);
}

/**
 * Copy the coefficient of constraint c into dst in PolyLib order,
 * i.e., first the coefficients of the variables, then the coefficients
 * of the parameters and finally the constant.
 */
void cloog_constraint_copy_coefficients(CloogConstraint constraint,
					cloog_int_t *dst)
{
	int i;
	unsigned dim;

	dim = isl_constraint_dim(constraint, isl_dim_all);

	for (i = 0; i < dim; ++i)
		cloog_constraint_coefficient_get(constraint, i, dst+i);
	cloog_constraint_constant_get(constraint, dst+dim);
}

CloogConstraint cloog_constraint_invalid(void)
{
	return NULL;
}

int cloog_constraint_is_valid(CloogConstraint constraint)
{
	return constraint != NULL;
}

int cloog_constraint_total_dimension(CloogConstraint constraint)
{
	return isl_constraint_dim(constraint, isl_dim_all);
}

/**
 * Create a CloogConstraintSet containing enough information to perform
 * a reduction on the upper equality (in this case lower is an invalid
 * CloogConstraint) or the pair of inequalities upper and lower
 * from within insert_modulo_guard.
 * In the isl backend, we return a CloogConstraintSet containting both
 * bounds, as the stride may change during the reduction and we may
 * need to recompute the bound on the modulo expression.
 */
CloogConstraintSet *cloog_constraint_set_for_reduction(CloogConstraint upper,
	CloogConstraint lower)
{
	CloogConstraintSet *set;

	set = isl_basic_set_from_constraint(isl_constraint_copy(upper));
	if (cloog_constraint_is_valid(lower))
		set = isl_basic_set_add_constraint(set,
						isl_constraint_copy(lower));
	return set;
}

/**
 * Reduce the modulo guard expressed by "contraints" using equalities
 * found in outer nesting levels (stored in "equal").
 * The modulo guard may be an equality or a pair of inequalities.
 * In case of a pair of inequalities, *bound contains the bound on the
 * corresponding modulo expression.  If any reduction is performed
 * then this bound is recomputed.
 *
 * We first check if "level" corresponds to an existentially quantified
 * variable.  If so, there is no need to reduce it as it would have
 * been removed already if it had been redundant.
 * Then we check if there are any equalities we can use.  If not,
 * there is again nothing to reduce.
 * For the actual reduction, we use isl_basic_set_gist, but this
 * function will only perform the reduction we want hear if the
 * the variable that imposes the modulo constraint has been projected
 * out (i.e., turned into an existentially quantified variable).
 * After the call to isl_basic_set_gist, we need to move the
 * existential variable back into the position where the calling
 * function expects it (assuming there are any constraints left).
 * We do this by adding equality between the given dimension and
 * the existentially quantified variable.
 */
CloogConstraintSet *cloog_constraint_set_reduce(CloogConstraintSet *constraints,
	int level, CloogEqualities *equal, int nb_par, cloog_int_t *bound)
{
	int j;
	struct isl_basic_set *eq;
	struct isl_basic_map *id;
	struct cloog_isl_dim dim;
	struct isl_constraint *c;
	struct isl_div *div;
	unsigned constraints_dim;
	int pos;
	isl_int v;

	dim = set_cloog_dim_to_isl_dim(constraints, level - 1);
	if (dim.type != isl_dim_set)
		return constraints;

	eq = NULL;
	for (j = 0; j < level - 1; ++j) {
		if (equal->types[j] != EQTYPE_EXAFFINE)
			continue;
		if (!eq)
			eq = isl_basic_set_copy(equal->constraints[j]);
		else
			eq = isl_basic_set_intersect(eq,
				isl_basic_set_copy(equal->constraints[j]));
	}
	if (!eq)
		return constraints;

	id = isl_basic_map_identity(isl_basic_set_get_dim(constraints));
	id = isl_basic_map_remove(id, isl_dim_out, dim.pos, 1);
	constraints = isl_basic_set_apply(constraints, isl_basic_map_copy(id));
	constraints = isl_basic_set_apply(constraints,
						isl_basic_map_reverse(id));

	constraints_dim = isl_basic_set_dim(constraints, isl_dim_set);
	eq = isl_basic_set_remove_dims(eq, constraints_dim,
			isl_basic_set_dim(eq, isl_dim_set) - constraints_dim);
	constraints = isl_basic_set_gist(constraints, eq);
	if (isl_basic_set_dim(constraints, isl_dim_div) != 1)
		return constraints;

	div = isl_basic_set_div(isl_basic_set_copy(constraints), 0);
	c = isl_equality_alloc(isl_basic_set_get_dim(constraints));
	c = isl_constraint_add_div(c, div, &pos);
	isl_constraint_set_coefficient(c, isl_dim_set, dim.pos,
					constraints->ctx->one);
	isl_constraint_set_coefficient(c, isl_dim_div, pos,
					constraints->ctx->negone);
	constraints = isl_basic_set_add_constraint(constraints, c);

	isl_int_init(v);
	isl_int_set_si(*bound, 0);
	for (c = cloog_constraint_first(constraints);
	     cloog_constraint_is_valid(c); c = cloog_constraint_next(c)) {
		cloog_constraint_constant_get(c, &v);
		isl_int_add(*bound, *bound, v);
	}
	isl_int_clear(v);

	return constraints;
}

CloogConstraint cloog_constraint_first(CloogConstraintSet *constraints)
{
	return isl_basic_set_first_constraint(isl_basic_set_copy(constraints));
}

CloogConstraint cloog_constraint_next(CloogConstraint constraint)
{
	return isl_constraint_next(constraint);
}

CloogConstraint cloog_constraint_copy(CloogConstraint constraint)
{
	return isl_constraint_copy(constraint);
}

void cloog_constraint_release(CloogConstraint constraint)
{
	isl_constraint_free(constraint);
}

CloogConstraint cloog_equal_constraint(CloogEqualities *equal, int j)
{
	return isl_basic_set_first_constraint(
			isl_basic_set_copy(equal->constraints[j]));
}
