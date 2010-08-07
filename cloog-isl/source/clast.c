#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../include/cloog/cloog.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n)*sizeof(type))

/**
 * CloogInfos structure:
 * this structure contains all the informations necessary for pretty printing,
 * they come from the original CloogProgram structure (language, names), from
 * genereral options (options) or are built only for pretty printing (stride).
 * This structure is mainly there to reduce the number of function parameters,
 * since most pprint.c functions need most of its field.
 */
struct clooginfos {
  CloogState *state;         /**< State. */
  cloog_int_t *stride;       /**< The stride for each iterator. */
  int  nb_scattdims ;        /**< Scattering dimension number. */
  int * scaldims ;           /**< Boolean array saying whether a given
                              *   scattering dimension is scalar or not.
                              */
  CloogNames * names ;       /**< Names of iterators and parameters. */
  CloogOptions * options ;   /**< Options on CLooG's behaviour. */
  CloogEqualities *equal;    /**< Matrix of equalities. */
} ;

typedef struct clooginfos CloogInfos ;

static int clast_expr_cmp(struct clast_expr *e1, struct clast_expr *e2);
static int clast_term_cmp(struct clast_term *t1, struct clast_term *t2);
static int clast_binary_cmp(struct clast_binary *b1, struct clast_binary *b2);
static int clast_reduction_cmp(struct clast_reduction *r1, 
				 struct clast_reduction *r2);

static int clast_equal_add(CloogEqualities *equal,
				CloogConstraintSet *constraints,
				int level, CloogConstraint constraint,
				CloogInfos *infos);

static struct clast_stmt *clast_equal(int level, CloogInfos *infos);
static struct clast_expr *clast_minmax(CloogConstraintSet *constraints,
					int level, int max, int guard, 
					CloogInfos *infos);
static void insert_guard(CloogConstraintSet *constraints, int level,
			  struct clast_stmt ***next, CloogInfos *infos);
static void insert_modulo_guard(CloogConstraint upper,
				CloogConstraint lower, int level,
			        struct clast_stmt ***next, CloogInfos *infos);
static void insert_equation(CloogConstraint upper, CloogConstraint lower,
			 int level, struct clast_stmt ***next, CloogInfos *infos);
static int insert_for(CloogConstraintSet *constraints, int level,
			struct clast_stmt ***next, CloogInfos *infos);
static void insert_block(CloogBlock *block, int level,
			  struct clast_stmt ***next, CloogInfos *infos);
static void insert_loop(CloogLoop * loop, int level, int scalar,
			struct clast_stmt ***next, CloogInfos *infos);


struct clast_name *new_clast_name(const char *name)
{
    struct clast_name *n = malloc(sizeof(struct clast_name));
    n->expr.type = clast_expr_name;
    n->name = name;
    return n;
}

struct clast_term *new_clast_term(cloog_int_t c, struct clast_expr *v)
{
    struct clast_term *t = malloc(sizeof(struct clast_term));
    t->expr.type = clast_expr_term;
    cloog_int_init(t->val);
    cloog_int_set(t->val, c);
    t->var = v;
    return t;
}

struct clast_binary *new_clast_binary(enum clast_bin_type t, 
				      struct clast_expr *lhs, cloog_int_t rhs)
{
    struct clast_binary *b = malloc(sizeof(struct clast_binary));
    b->expr.type = clast_expr_bin;
    b->type = t;
    b->LHS = lhs;
    cloog_int_init(b->RHS);
    cloog_int_set(b->RHS, rhs);
    return b;
}

struct clast_reduction *new_clast_reduction(enum clast_red_type t, int n)
{
    int i;
    struct clast_reduction *r;
    r = malloc(sizeof(struct clast_reduction)+(n-1)*sizeof(struct clast_expr *));
    r->expr.type = clast_expr_red;
    r->type = t;
    r->n = n;
    for (i = 0; i < n; ++i)
	r->elts[i] = NULL;
    return r;
}

static void free_clast_root(struct clast_stmt *s);

const struct clast_stmt_op stmt_root = { free_clast_root };

static void free_clast_root(struct clast_stmt *s)
{
    struct clast_root *r = (struct clast_root *)s;
    assert(CLAST_STMT_IS_A(s, stmt_root));
    cloog_names_free(r->names);
    free(r);
}

struct clast_root *new_clast_root(CloogNames *names)
{
    struct clast_root *r = malloc(sizeof(struct clast_root));
    r->stmt.op = &stmt_root;
    r->stmt.next = NULL;
    r->names = cloog_names_copy(names);
    return r;
}

static void free_clast_assignment(struct clast_stmt *s);

const struct clast_stmt_op stmt_ass = { free_clast_assignment };

static void free_clast_assignment(struct clast_stmt *s)
{
    struct clast_assignment *a = (struct clast_assignment *)s;
    assert(CLAST_STMT_IS_A(s, stmt_ass));
    free_clast_expr(a->RHS);
    free(a);
}

struct clast_assignment *new_clast_assignment(const char *lhs,
					      struct clast_expr *rhs)
{
    struct clast_assignment *a = malloc(sizeof(struct clast_assignment));
    a->stmt.op = &stmt_ass;
    a->stmt.next = NULL;
    a->LHS = lhs;
    a->RHS = rhs;
    return a;
}

static void free_clast_user_stmt(struct clast_stmt *s);

const struct clast_stmt_op stmt_user = { free_clast_user_stmt };

static void free_clast_user_stmt(struct clast_stmt *s)
{
    struct clast_user_stmt *u = (struct clast_user_stmt *)s;
    assert(CLAST_STMT_IS_A(s, stmt_user));
    cloog_statement_free(u->statement);
    cloog_clast_free(u->substitutions);
    free(u);
}

struct clast_user_stmt *new_clast_user_stmt(CloogStatement *stmt, 
					    struct clast_stmt *subs)
{
    struct clast_user_stmt *u = malloc(sizeof(struct clast_user_stmt));
    u->stmt.op = &stmt_user;
    u->stmt.next = NULL;
    u->statement = cloog_statement_copy(stmt);
    u->substitutions = subs;
    return u;
}

static void free_clast_block(struct clast_stmt *b);

const struct clast_stmt_op stmt_block = { free_clast_block };

static void free_clast_block(struct clast_stmt *s)
{
    struct clast_block *b = (struct clast_block *)s;
    assert(CLAST_STMT_IS_A(s, stmt_block));
    cloog_clast_free(b->body);
    free(b);
}

struct clast_block *new_clast_block()
{
    struct clast_block *b = malloc(sizeof(struct clast_block));
    b->stmt.op = &stmt_block;
    b->stmt.next = NULL;
    b->body = NULL;
    return b;
}

static void free_clast_for(struct clast_stmt *s);

const struct clast_stmt_op stmt_for = { free_clast_for };

static void free_clast_for(struct clast_stmt *s)
{
    struct clast_for *f = (struct clast_for *)s;
    assert(CLAST_STMT_IS_A(s, stmt_for));
    free_clast_expr(f->LB);
    free_clast_expr(f->UB);
    cloog_int_clear(f->stride);
    cloog_clast_free(f->body);
    free(f);
}

struct clast_for *new_clast_for(const char *it, struct clast_expr *LB, 
				struct clast_expr *UB, cloog_int_t stride)
{
    struct clast_for *f = malloc(sizeof(struct clast_for));
    f->stmt.op = &stmt_for;
    f->stmt.next = NULL;
    f->iterator = it;
    f->LB = LB;
    f->UB = UB;
    f->body = NULL;
    cloog_int_init(f->stride);
    cloog_int_set(f->stride, stride);
    return f;
}

static void free_clast_guard(struct clast_stmt *s);

const struct clast_stmt_op stmt_guard = { free_clast_guard };

static void free_clast_guard(struct clast_stmt *s)
{
    int i;
    struct clast_guard *g = (struct clast_guard *)s;
    assert(CLAST_STMT_IS_A(s, stmt_guard));
    cloog_clast_free(g->then);
    for (i = 0; i < g->n; ++i) {
	free_clast_expr(g->eq[i].LHS);
	free_clast_expr(g->eq[i].RHS);
    }
    free(g);
}

struct clast_guard *new_clast_guard(int n)
{
    int i;
    struct clast_guard *g = malloc(sizeof(struct clast_guard) + 
				   (n-1) * sizeof(struct clast_equation));
    g->stmt.op = &stmt_guard;
    g->stmt.next = NULL;
    g->then = NULL;
    g->n = n;
    for (i = 0; i < n; ++i) {
	g->eq[i].LHS = NULL;
	g->eq[i].RHS = NULL;
    }
    return g;
}

void free_clast_name(struct clast_name *n)
{
    free(n);
}

void free_clast_term(struct clast_term *t)
{
    cloog_int_clear(t->val);
    free_clast_expr(t->var);
    free(t);
}

void free_clast_binary(struct clast_binary *b)
{
    cloog_int_clear(b->RHS);
    free_clast_expr(b->LHS);
    free(b);
}

void free_clast_reduction(struct clast_reduction *r)
{
    int i;
    for (i = 0; i < r->n; ++i)
	free_clast_expr(r->elts[i]);
    free(r);
}

void free_clast_expr(struct clast_expr *e)
{
    if (!e)
	return;
    switch (e->type) {
    case clast_expr_name:
	free_clast_name((struct clast_name*) e);
	break;
    case clast_expr_term:
	free_clast_term((struct clast_term*) e);
	break;
    case clast_expr_red:
	free_clast_reduction((struct clast_reduction*) e);
	break;
    case clast_expr_bin:
	free_clast_binary((struct clast_binary*) e);
	break;
    default:
	assert(0);
    }
}

void free_clast_stmt(struct clast_stmt *s)
{
    assert(s->op);
    assert(s->op->free);
    s->op->free(s);
}

void cloog_clast_free(struct clast_stmt *s)
{
    struct clast_stmt *next;
    while (s) {
	next = s->next;
	free_clast_stmt(s);
	s = next;
    }
}

static int clast_name_cmp(struct clast_name *n1, struct clast_name *n2)
{
    return n1->name == n2->name ? 0 : strcmp(n1->name, n2->name);
}

static int clast_term_cmp(struct clast_term *t1, struct clast_term *t2)
{
    int c;
    if (!t1->var && t2->var)
	return -1;
    if (t1->var && !t2->var)
	return 1;
    c = clast_expr_cmp(t1->var, t2->var);
    if (c)
	return c;
    return cloog_int_cmp(t1->val, t2->val);
}

static int clast_binary_cmp(struct clast_binary *b1, struct clast_binary *b2)
{
    int c;

    if (b1->type != b2->type)
	return b1->type - b2->type;
    if ((c = cloog_int_cmp(b1->RHS, b2->RHS)))
	return c;
    return clast_expr_cmp(b1->LHS, b2->LHS);
}

static int clast_reduction_cmp(struct clast_reduction *r1, struct clast_reduction *r2)
{
    int i;
    int c;

    if (r1->n == 1 && r2->n == 1)
	return clast_expr_cmp(r1->elts[0], r2->elts[0]);
    if (r1->type != r2->type)
	return r1->type - r2->type;
    if (r1->n != r2->n)
	return r1->n - r2->n;
    for (i = 0; i < r1->n; ++i)
	if ((c = clast_expr_cmp(r1->elts[i], r2->elts[i])))
	    return c;
    return 0;
}

static int clast_expr_cmp(struct clast_expr *e1, struct clast_expr *e2)
{
    if (!e1 && !e2)
	return 0;
    if (!e1)
	return -1;
    if (!e2)
	return 1;
    if (e1->type != e2->type)
	return e1->type - e2->type;
    switch (e1->type) {
    case clast_expr_name:
	return clast_name_cmp((struct clast_name*) e1, 
				(struct clast_name*) e2);
    case clast_expr_term:
	return clast_term_cmp((struct clast_term*) e1, 
				(struct clast_term*) e2);
    case clast_expr_bin:
	return clast_binary_cmp((struct clast_binary*) e1, 
				(struct clast_binary*) e2);
    case clast_expr_red:
	return clast_reduction_cmp((struct clast_reduction*) e1, 
				   (struct clast_reduction*) e2);
    default:
	assert(0);
    }
}

int clast_expr_equal(struct clast_expr *e1, struct clast_expr *e2)
{
    return clast_expr_cmp(e1, e2) == 0;
}

/**
 * Return 1 is both expressions are constant terms and e1 is bigger than e2.
 */
int clast_expr_is_bigger_constant(struct clast_expr *e1, struct clast_expr *e2)
{
    struct clast_term *t1, *t2;
    struct clast_reduction *r;

    if (!e1 || !e2)
	return 0;
    if (e1->type == clast_expr_red) {
	r = (struct clast_reduction *)e1;
	return r->n == 1 && clast_expr_is_bigger_constant(r->elts[0], e2);
    }
    if (e2->type == clast_expr_red) {
	r = (struct clast_reduction *)e2;
	return r->n == 1 && clast_expr_is_bigger_constant(e1, r->elts[0]);
    }
    if (e1->type != clast_expr_term || e2->type != clast_expr_term)
	return 0;
    t1 = (struct clast_term *)e1;
    t2 = (struct clast_term *)e2;
    if (t1->var || t2->var)
	return 0;
    return cloog_int_gt(t1->val, t2->val);
}

static int qsort_expr_cmp(const void *p1, const void *p2)
{
    return clast_expr_cmp(*(struct clast_expr **)p1, *(struct clast_expr **)p2);
}

static void clast_reduction_sort(struct clast_reduction *r)
{
    qsort(&r->elts[0], r->n, sizeof(struct clast_expr *), qsort_expr_cmp);
}

static int qsort_eq_cmp(const void *p1, const void *p2)
{
    struct clast_equation *eq1 = (struct clast_equation *)p1;
    struct clast_equation *eq2 = (struct clast_equation *)p2;
    int cmp;

    cmp = clast_expr_cmp(eq1->LHS, eq2->LHS);
    if (cmp)
	return cmp;

    cmp = clast_expr_cmp(eq1->RHS, eq2->RHS);
    if (cmp)
	return cmp;

    return eq1->sign - eq2->sign;
}

/**
 * Sort equations in a clast_guard.
 */
static void clast_guard_sort(struct clast_guard *g)
{
    qsort(&g->eq[0], g->n, sizeof(struct clast_equation), qsort_eq_cmp);
}


/******************************************************************************
 *                        Equalities spreading functions                      *
 ******************************************************************************/


/**
 * clast_equal_allow function:
 * This function checks whether the options allow us to spread the equality or
 * not. It returns 1 if so, 0 otherwise.
 * - equal is the matrix of equalities,
 * - level is the column number in equal of the element which is 'equal to',
 * - line is the line number in equal of the constraint we want to study,
 * - the infos structure gives the user all options on code printing and more.
 **
 * - October 27th 2005: first version (extracted from old pprint_equal_add).
 */
static int clast_equal_allow(CloogEqualities *equal, int level, int line,
				CloogInfos *infos)
{ 
  if (level < infos->options->fsp)
  return 0 ;
  
  if ((cloog_equal_type(equal, level) == EQTYPE_EXAFFINE) &&
      !infos->options->esp)
  return 0 ;

  return 1 ;
}


/**
 * clast_equal_add function:
 * This function updates the row (level-1) of the equality matrix (equal) with
 * the row that corresponds to the row (line) of the matrix (matrix). It returns
 * 1 if the row can be updated, 0 otherwise.
 * - equal is the matrix of equalities,
 * - matrix is the matrix of constraints,
 * - level is the column number in matrix of the element which is 'equal to',
 * - line is the line number in matrix of the constraint we want to study,
 * - the infos structure gives the user all options on code printing and more.
 */
static int clast_equal_add(CloogEqualities *equal,
				CloogConstraintSet *constraints,
				int level, CloogConstraint constraint,
				CloogInfos *infos)
{
    cloog_equal_add(equal, constraints, level, constraint,
		    infos->names->nb_parameters);
  
    return clast_equal_allow(equal, level, level-1, infos);
}



/**
 * clast_equal function:
 * This function prints the substitution data of a statement into a clast_stmt.
 * Using this function instead of pprint_equal is useful for generating
 * a compilable pseudo-code by using preprocessor macro for each statement.
 * By opposition to pprint_equal, the result is less human-readable. For
 * instance this function will print (i,i+3,k,3) where pprint_equal would
 * return (j=i+3,l=3).
 * - level is the number of loops enclosing the statement,
 * - the infos structure gives the user all options on code printing and more.
 **
 * - March    12th 2004: first version. 
 * - November 21th 2005: (debug) now works well with GMP version. 
 */
static struct clast_stmt *clast_equal(int level, CloogInfos *infos)
{ 
  int i ;
  struct clast_expr *e;
  struct clast_stmt *a = NULL;
  struct clast_stmt **next = &a;
  CloogEqualities *equal = infos->equal;
  CloogConstraint equal_constraint;

  for (i=infos->names->nb_scattering;i<level-1;i++)
  { if (cloog_equal_type(equal, i+1)) {
      equal_constraint = cloog_equal_constraint(equal, i);
      e = clast_bound_from_constraint(equal_constraint, i+1, infos->names);
      cloog_constraint_release(equal_constraint);
    } else {
      e = &new_clast_term(infos->state->one, &new_clast_name(
		 cloog_names_name_at_level(infos->names, i+1))->expr)->expr;
    }
    *next = &new_clast_assignment(NULL, e)->stmt;
    next = &(*next)->next;
  }

  return a;
}

 
/**
 * clast_bound_from_constraint function:
 * This function returns a clast_expr containing the printing of the
 * 'right part' of a constraint according to an element.
 * For instance, for the constraint -3*i + 2*j - M >=0 and the element j,
 * we have j >= (3*i + M)/2. As we are looking for integral solutions, this
 * function should return 'ceild(3*i+M,2)'.
 * - matrix is the polyhedron containing all the constraints,
 * - line_num is the line number in domain of the constraint we want to print,
 * - level is the column number in domain of the element we want to use,
 * - names structure gives the user some options about code printing,
 *   the number of parameters in domain (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - November 2nd 2001: first version. 
 * - June    27th 2003: 64 bits version ready.
 */
struct clast_expr *clast_bound_from_constraint(CloogConstraint constraint,
					       int level, CloogNames *names)
{ 
  int i, sign, nb_elts=0, len;
  cloog_int_t *line, numerator, denominator, temp, division;
  struct clast_expr *e = NULL;
  struct cloog_vec *line_vector;

  len = cloog_constraint_total_dimension(constraint) + 2;
  line_vector = cloog_vec_alloc(len);
  line = line_vector->p;
  cloog_constraint_copy_coefficients(constraint, line+1);
  cloog_int_init(temp);
  cloog_int_init(numerator);
  cloog_int_init(denominator);

  if (!cloog_int_is_zero(line[level])) {
    struct clast_reduction *r;
    /* Maybe we need to invert signs in such a way that the element sign is>0.*/
    sign = -cloog_int_sgn(line[level]);

    for (i = 1, nb_elts = 0; i <= len - 1; ++i)
	if (i != level && !cloog_int_is_zero(line[i]))
	    nb_elts++;
    r = new_clast_reduction(clast_red_sum, nb_elts);
    nb_elts = 0;

    /* First, we have to print the iterators and the parameters. */
    for (i = 1; i <= len - 2; i++) {
      struct clast_expr *v;

      if (i == level || cloog_int_is_zero(line[i]))
	continue;

      v = cloog_constraint_variable_expr(constraint, i, names);
      
      if (sign == -1)
	cloog_int_neg(temp,line[i]);
      else
	cloog_int_set(temp,line[i]);
      
      r->elts[nb_elts++] = &new_clast_term(temp, v)->expr;
    }    

    if (sign == -1) {
      cloog_int_neg(numerator, line[len - 1]);
      cloog_int_set(denominator, line[level]);
    }
    else {
      cloog_int_set(numerator, line[len - 1]);
      cloog_int_neg(denominator, line[level]);
    }
        
    /* Finally, the constant, and the final printing. */
    if (nb_elts) {
      if (!cloog_int_is_zero(numerator))
	  r->elts[nb_elts++] = &new_clast_term(numerator, NULL)->expr;
    
      if (!cloog_int_is_one(line[level]) && !cloog_int_is_neg_one(line[level]))
      { if (!cloog_constraint_is_equality(constraint))
        { if (cloog_int_is_pos(line[level]))
	    e = &new_clast_binary(clast_bin_cdiv, &r->expr, denominator)->expr;
          else
	    e = &new_clast_binary(clast_bin_fdiv, &r->expr, denominator)->expr;
        } else
	    e = &new_clast_binary(clast_bin_div, &r->expr, denominator)->expr;
      }
      else
	e = &r->expr;
    } else { 
      free_clast_reduction(r);
      if (cloog_int_is_zero(numerator))
	e = &new_clast_term(numerator, NULL)->expr;
      else
      { if (!cloog_int_is_one(denominator))
        { if (!cloog_constraint_is_equality(constraint)) { /* useful? */
            if (cloog_int_is_divisible_by(numerator, denominator)) {
              cloog_int_divexact(temp, numerator, denominator);
	      e = &new_clast_term(temp, NULL)->expr;
            }
            else {
              cloog_int_init(division);
	      cloog_int_tdiv_q(division, numerator, denominator);
	      if (cloog_int_is_neg(numerator)) {
                if (cloog_int_is_pos(line[level])) {
		    /* nb<0 need max */
		    e = &new_clast_term(division, NULL)->expr;
		} else {
                  /* nb<0 need min */
                  cloog_int_sub_ui(temp, division, 1);
		  e = &new_clast_term(temp, NULL)->expr;
                }
	      }
              else
              { if (cloog_int_is_pos(line[level]))
	        { /* nb>0 need max */
                  cloog_int_add_ui(temp, division, 1);
		  e = &new_clast_term(temp, NULL)->expr;
                }
		else
		    /* nb>0 need min */
		    e = &new_clast_term(division, NULL)->expr;
              }
	      cloog_int_clear(division);
            }
          }
          else
	    e = &new_clast_binary(clast_bin_div, 
				  &new_clast_term(numerator, NULL)->expr,
				  denominator)->expr;
        }
        else
	    e = &new_clast_term(numerator, NULL)->expr;
      }
    }
  }

  cloog_vec_free(line_vector);

  cloog_int_clear(temp);
  cloog_int_clear(numerator);
  cloog_int_clear(denominator);
    
  return e;
}


/**
 * clast_minmax function:
 * This function returns a clast_expr containing the printing of a minimum or a
 * maximum of the 'right parts' of all constraints according to an element.
 * For instance consider the constraints:
 * -3*i  +2*j   -M >= 0
 *  2*i    +j      >= 0
 *   -i    -j +2*M >= 0
 * if we are looking for the minimum for the element j, the function should
 * return 'max(ceild(3*i+M,2),-2*i)'.
 * - constraints is the constraints,
 * - level is the column number in domain of the element we want to use,
 * - max is a boolean set to 1 if we are looking for a maximum, 0 for a minimum,
 * - guard is set to 0 if there is no guard, and set to the level of the element
 *   with a guard otherwise (then the function gives the max or the min only
 *   for the constraint where the guarded coefficient is 0), 
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in domain (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - November 2nd 2001: first version. 
 */
static struct clast_expr *clast_minmax(CloogConstraintSet *constraints,
				       int level, int max, int guard,
				       CloogInfos *infos)
{ int n;
  struct clast_reduction *r;
  CloogConstraint constraint;
  
  n = 0;
  for (constraint = cloog_constraint_first(constraints);
       cloog_constraint_is_valid(constraint);
       constraint = cloog_constraint_next(constraint))
      if (((max && cloog_constraint_is_lower_bound(constraint, level-1)) ||
	   (!max && cloog_constraint_is_upper_bound(constraint, level-1))) &&
	  (!guard || !cloog_constraint_involves(constraint, guard-1)) &&
	  (!cloog_constraint_is_equality(constraint)))
	n++;
  if (!n)
    return NULL;
  r = new_clast_reduction(max ? clast_red_max : clast_red_min, n);

  n = 0;
  for (constraint = cloog_constraint_first(constraints);
       cloog_constraint_is_valid(constraint);
       constraint = cloog_constraint_next(constraint))
      if (((max && cloog_constraint_is_lower_bound(constraint, level-1)) ||
	   (!max && cloog_constraint_is_upper_bound(constraint, level-1))) &&
	  (!guard || !cloog_constraint_involves(constraint, guard-1)) &&
	  (!cloog_constraint_is_equality(constraint)))
	r->elts[n++] = clast_bound_from_constraint(constraint, level,
								infos->names);

  clast_reduction_sort(r);
  return &r->expr;
}


/**
 * Insert modulo guards defined by existentially quantified dimensions,
 * not involving the given level.
 *
 * This function is called from within insert_guard or insert_for.
 * Any constraint used in constructing * a modulo guard is removed
 * from the constraint set to avoid insert_guard or insert_for
 * adding a duplicate (pair of) constraint(s).
 */
static void insert_extra_modulo_guards(CloogConstraintSet *constraints,
		int level, struct clast_stmt ***next, CloogInfos *infos)
{
    int i;
    int nb_iter;
    int total_dim;
    CloogConstraint upper, lower;

    total_dim = cloog_constraint_set_total_dimension(constraints);
    nb_iter = cloog_constraint_set_n_iterators(constraints,
						infos->names->nb_parameters);

    for (i = total_dim - infos->names->nb_parameters; i >= nb_iter + 1; i--) {
	if (cloog_constraint_is_valid(upper =
		cloog_constraint_set_defining_equality(constraints, i))) {
	    if (!level || (nb_iter < level) ||
		    !cloog_constraint_involves(upper, level-1)) {
		insert_modulo_guard(upper,
				cloog_constraint_invalid(), i, next, infos);
		cloog_constraint_clear(upper);
	    }
	    cloog_constraint_release(upper);
	} else if (cloog_constraint_is_valid(upper =
		    cloog_constraint_set_defining_inequalities(constraints,
			      i, &lower, infos->names->nb_parameters))) {
	    if (!level || (nb_iter < level) ||
		    !cloog_constraint_involves(upper, level-1)) {
		insert_modulo_guard(upper, lower, i, next, infos);
		cloog_constraint_clear(upper);
		cloog_constraint_clear(lower);
	    }
	    cloog_constraint_release(upper);
	    cloog_constraint_release(lower);
	}
    }
}


/**
 * insert_guard function:
 * This function inserts a guard in the clast.
 * A guard on an element (level) is :
 * -> the conjunction of all the existing constraints where the coefficient of
 *    this element is 0 if the element is an iterator,
 * -> the conjunction of all the existing constraints if the element isn't an
 *    iterator.
 * For instance, considering these constraints and the element j:
 * -3*i +2*j -M >= 0
 *  2*i      +M >= 0
 * this function should return 'if (2*i+M>=0) {'.
 * - matrix is the polyhedron containing all the constraints,
 * - level is the column number of the element in matrix we want to use,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in matrix (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - November  3rd 2001: first version. 
 * - November 14th 2001: a lot of 'purifications'. 
 * - July     31th 2002: (debug) some guard parts are no more redundants. 
 * - August   12th 2002: polyhedra union ('or' conditions) are now supported.
 * - October  27th 2005: polyhedra union ('or' conditions) are no more supported
 *                       (the need came from loop_simplify that may result in
 *                       domain unions, now it should be fixed directly in
 *                       cloog_loop_simplify).
 */
static void insert_guard(CloogConstraintSet *constraints, int level,
			 struct clast_stmt ***next, CloogInfos *infos)
{ 
  int i, guarded, minmax=-1, nb_and = 0, nb_iter ;
  int total_dim;
  CloogConstraintSet *copy;
  CloogConstraint j, l;
  struct clast_guard *g;

  if (constraints == NULL)
    return;

    copy = cloog_constraint_set_copy(constraints);

    insert_extra_modulo_guards(copy, level, next, infos);
  
    total_dim = cloog_constraint_set_total_dimension(constraints);
    g = new_clast_guard(2 * total_dim);

    /* Well, it looks complicated because I wanted to have a particular, more
     * readable, ordering, obviously this function may be far much simpler !
     */
    nb_iter = cloog_constraint_set_n_iterators(constraints,
						infos->names->nb_parameters);
 
    nb_and = 0 ;
    /* We search for guard parts. */
    for (i = 1; i <= total_dim; i++)
      for (j = cloog_constraint_first(copy); cloog_constraint_is_valid(j);
	   j = cloog_constraint_next(j))
	if (cloog_constraint_involves(j, i-1) &&
	    (!level || (nb_iter < level) ||
	     !cloog_constraint_involves(j, level-1))) {
	  struct clast_expr *v;
	  struct clast_term *t;

	  v = cloog_constraint_variable_expr(j, i, infos->names);
	  g->eq[nb_and].LHS = &(t = new_clast_term(infos->state->one, v))->expr;
	  if (!level || cloog_constraint_is_equality(j)) {
	    /* put the "denominator" in the LHS */
	    cloog_constraint_coefficient_get(j, i-1, &t->val);
	    cloog_constraint_coefficient_set(j, i-1, infos->state->one);
	    if (cloog_int_is_neg(t->val)) {
	      cloog_int_neg(t->val, t->val);
	      cloog_constraint_coefficient_set(j, i-1, infos->state->negone);
	    }
	    if (level || cloog_constraint_is_equality(j))
	      g->eq[nb_and].sign = 0;
	    else if (cloog_constraint_is_lower_bound(j, i-1))
	      g->eq[nb_and].sign = 1;
	    else
	      g->eq[nb_and].sign = -1;
	    g->eq[nb_and].RHS = clast_bound_from_constraint(j, i, infos->names);
	  } else {
	    if (cloog_constraint_is_lower_bound(j, i-1)) {
		minmax = 1;
		g->eq[nb_and].sign = 1;
	    } else {
		minmax = 0;
		g->eq[nb_and].sign = -1;
	    }
	  
	    guarded = (nb_iter >= level) ? level : 0 ;
	    g->eq[nb_and].RHS = clast_minmax(copy,i,minmax,guarded,infos) ;
	  }
	  nb_and ++ ;

	  /* 'elimination' of the current constraint, this avoid to use one
	   * constraint more than once. The current line is always eliminated,
	   * and the next lines if they are in a min or a max.
	   */
	  cloog_constraint_clear(j);
	
	  if (minmax == -1)
	    continue;
	  l = cloog_constraint_copy(j);
	  for (l = cloog_constraint_next(l); cloog_constraint_is_valid(l);
	       l = cloog_constraint_next(l))
	    if (((minmax == 1) && cloog_constraint_is_lower_bound(l, i-1)) ||
		((minmax == 0) && cloog_constraint_is_upper_bound(l, i-1)))
	      cloog_constraint_clear(l);
	}
	cloog_constraint_set_free(copy);
  
  g->n = nb_and;
  if (nb_and) {
    clast_guard_sort(g);
    **next = &g->stmt;
    *next = &g->then;
  } else
    free_clast_stmt(&g->stmt);
  
  return;
}
 

/**
 * insert_modulo_guard:
 * This function inserts a modulo guard corresponding to an equality
 * or a pair of inequalities.
 * See insert_equation.
 * - matrix is the polyhedron containing all the constraints,
 * - upper and lower are the line numbers of the constraint in matrix
 *   we want to print; in particular, if we want to print an equality,
 *   then lower == -1 and upper is the row of the equality; if we want
 *   to print an inequality, then upper is the row of the upper bound
 *   and lower in the row of the lower bound
 * - level is the column number of the element in matrix we want to use,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in matrix (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 */
static void insert_modulo_guard(CloogConstraint upper,
				CloogConstraint lower, int level,
				struct clast_stmt ***next, CloogInfos *infos)
{
  int i, nb_elts = 0, len, len2, nb_iter, in_stride = 0, nb_par;
  struct cloog_vec *line_vector;
  cloog_int_t *line, val, bound;
  CloogConstraintSet *set;

  cloog_int_init(val);
  cloog_constraint_coefficient_get(upper, level-1, &val);
  if (cloog_int_is_one(val) || cloog_int_is_neg_one(val)) {
    cloog_int_clear(val);
    return;
  }

  len = cloog_constraint_total_dimension(upper) + 2;
  len2 = cloog_equal_total_dimension(infos->equal) + 2;
  nb_par = infos->names->nb_parameters;
  nb_iter = len - 2 - nb_par;

  cloog_int_init(bound);
  /* Check if would be emitting the redundant constraint mod(e,m) <= m-1 */
  if (cloog_constraint_is_valid(lower)) {
    cloog_constraint_constant_get(upper, &val);
    cloog_constraint_constant_get(lower, &bound);
    cloog_int_add(bound, val, bound);
    cloog_constraint_coefficient_get(lower, level-1, &val);
    cloog_int_sub_ui(val, val, 1);
    if (cloog_int_eq(val, bound)) {
      cloog_int_clear(val);
      cloog_int_clear(bound);
      return;
    }
  }

  set = cloog_constraint_set_for_reduction(upper, lower);
  set = cloog_constraint_set_reduce(set, level, infos->equal, nb_par, &bound);
  upper = cloog_constraint_first(set);
  if (!cloog_constraint_is_valid(upper)) {
    cloog_int_clear(val);
    cloog_int_clear(bound);
    cloog_constraint_set_free(set);
    return;
  }

  line_vector = cloog_vec_alloc(len);
  line = line_vector->p;
  cloog_constraint_copy_coefficients(upper, line+1);
  if (cloog_int_is_pos(line[level]))
    cloog_seq_neg(line+1, line+1, len-1);
  cloog_int_neg(line[level], line[level]);
  assert(cloog_int_is_pos(line[level]));

  nb_elts = 0;
  for (i = 1; i <= len-1; ++i) {
    if (i == level)
      continue;
    cloog_int_fdiv_r(line[i], line[i], line[level]);
    if (cloog_int_is_zero(line[i]))
      continue;
    if (i == len-1)
      continue;

    /* We need to know if an element of the equality has not to be printed
     * because of a stride that guarantees that this element can be divided by
     * the current coefficient. Because when there is a constant element, it
     * is included in the stride calculation (more exactly in the strided
     * iterator new lower bound: the 'offset') and we have not to print it.
     */
    if (i <= nb_iter && !cloog_constraint_is_valid(lower) &&
	cloog_int_is_divisible_by(infos->stride[i-1], line[level])) {
      in_stride = 1;
      continue;
    }

    nb_elts++;
  }

  if (nb_elts || (!cloog_int_is_zero(line[len-1]) && (!in_stride))) {
    struct clast_reduction *r;
    struct clast_expr *e;
    struct clast_guard *g;
    const char *name;

    r = new_clast_reduction(clast_red_sum, nb_elts+1);
    nb_elts = 0;

    /* First, the modulo guard : the iterators... */
    for (i=1;i<=nb_iter;i++) {
      if (i == level || cloog_int_is_zero(line[i]))
	continue;
      if (cloog_int_is_divisible_by(infos->stride[i-1], line[level]))
	continue;

      name = cloog_names_name_at_level(infos->names, i);

      r->elts[nb_elts++] = &new_clast_term(line[i],
				&new_clast_name(name)->expr)->expr;
    }

    /* ...the parameters... */
    for (i=nb_iter+1;i<=len-2;i++) {
      if (cloog_int_is_zero(line[i]))
	continue;

      name = infos->names->parameters[i-nb_iter-1] ;
      r->elts[nb_elts++] = &new_clast_term(line[i],
				&new_clast_name(name)->expr)->expr;
    }

    /* ...the constant. */
    if (!cloog_int_is_zero(line[len-1]))
      r->elts[nb_elts++] = &new_clast_term(line[len-1], NULL)->expr;

    /* our initial computation may have been an overestimate */
    r->n = nb_elts;

    e = &new_clast_binary(clast_bin_mod, &r->expr, line[level])->expr;
    g = new_clast_guard(1);
    if (!cloog_constraint_is_valid(lower)) {
      g->eq[0].LHS = e;
      cloog_int_set_si(val, 0);
      g->eq[0].RHS = &new_clast_term(val, NULL)->expr;
      g->eq[0].sign = 0;
    } else {
      g->eq[0].LHS = e;
      g->eq[0].RHS = &new_clast_term(bound, NULL)->expr;
      g->eq[0].sign = -1;
    }

    **next = &g->stmt;
    *next = &g->then;
  }

  cloog_constraint_release(upper);
  cloog_constraint_set_free(set);
  cloog_vec_free(line_vector);
  cloog_int_clear(val);
  cloog_int_clear(bound);
}


/**
 * insert_equation function:
 * This function inserts an equality 
 * constraint according to an element in the clast.
 * An equality can be preceded by a 'modulo guard'.
 * For instance, consider the constraint i -2*j = 0 and the
 * element j: pprint_equality should return 'if(i%2==0) { j = i/2 ;'.
 * - matrix is the polyhedron containing all the constraints,
 * - num is the line number of the constraint in matrix we want to print,
 * - level is the column number of the element in matrix we want to use,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in matrix (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - November 13th 2001: first version.
 * - June 26th 2003: simplification of the modulo guards (remove parts such as
 *                   modulo is 0, compare vivien or vivien2 with a previous
 *                   version for an idea).
 * - June 29th 2003: non-unit strides support.
 * - July 14th 2003: (debug) no more print the constant in the modulo guard when
 *                   it was previously included in a stride calculation.
 */
static void insert_equation(CloogConstraint upper, CloogConstraint lower,
			    int level, struct clast_stmt ***next, CloogInfos *infos)
{
  struct clast_expr *e;
  struct clast_assignment *ass;

  insert_modulo_guard(upper, lower, level, next, infos);

  if (cloog_constraint_is_valid(lower) ||
      !clast_equal_add(infos->equal, NULL, level, upper, infos))
  { /* Finally, the equality. */
		
    /* If we have to make a block by dimension, we start the block. Function
     * pprint knows if there is an equality, if this is the case, it checks
     * for the same following condition to close the brace.
     */
    if (infos->options->block) {
      struct clast_block *b = new_clast_block();
      **next = &b->stmt;
      *next = &b->body;
    }
		
    e = clast_bound_from_constraint(upper, level, infos->names);
    ass = new_clast_assignment(cloog_names_name_at_level(infos->names, level), e);

    **next = &ass->stmt;
    *next = &(**next)->next;
  }

  cloog_constraint_release(lower);
  cloog_constraint_release(upper);

  return;
}


/**
 * insert_for function:
 * This function inserts a for loop in the clast.
 * Returns 1 if the calling function should recurse into inner loops.
 *
 * A loop header according to an element is the conjunction of a minimum and a
 * maximum on a given element (they give the loop bounds).
 * For instance, considering these constraints and the element j:
 * i + j -9*M >= 0
 *    -j +5*M >= 0
 *     j -4*M >= 0
 * this function should return 'for (j=max(-i+9*M,4*M),j<=5*M;j++) {'.
 * If the given element is involved in modulo guards defined by
 * existentially quantified variables, then these guards should be
 * inserted inside the for loop.  However, the constraints involved
 * in this guard should not be used in determining the lower and upper
 * bound of the loop.  We therefore insert the guards first (which
 * removes the corresponding constraints from the constraint set)
 * and then reattach the guard inside the loop.
 * - constraints contains all constraints,
 * - level is the column number of the element in matrix we want to use,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in matrix (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 */
static int insert_for(CloogConstraintSet *constraints, int level,
		       struct clast_stmt ***next, CloogInfos *infos)
{
  const char *iterator;
  struct clast_expr *e1;
  struct clast_expr *e2;
  struct clast_assignment *ass;
  struct clast_stmt **old_next = *next;
  struct clast_stmt *guard;
  
  insert_extra_modulo_guards(constraints, 0, next, infos);
  guard = *old_next;

  iterator = cloog_names_name_at_level(infos->names, level);
  
  e1 = clast_minmax(constraints, level, 1, 0, infos);
  e2 = clast_minmax(constraints, level, 0, 0, infos);

  if (clast_expr_is_bigger_constant(e1, e2)) {
    free_clast_expr(e1);
    free_clast_expr(e2);
    return 0;
  }

  /* If min and max are not equal there is a 'for' else, there is a '='.
   * In the special case e1 = e2 = NULL, this is an infinite loop
   * so this is not a '='.
   */
  if (!clast_expr_equal(e1, e2) || !infos->options->otl || (!e1 && !e2)) {
    struct clast_for *f = new_clast_for(iterator, e1, e2, infos->stride[level-1]);
    *old_next = &f->stmt;
    if (guard)
	f->body = guard;
    else
	*next = &f->body;
  }
  else if (!clast_equal_add(infos->equal, constraints, level,
				cloog_constraint_invalid(), infos)) {
    if (infos->options->block) {
	struct clast_block *b = new_clast_block();
	*old_next = &b->stmt;
	if (guard)
	    b->body = guard;
	else
	    *next = &b->body;
    }
    ass = new_clast_assignment(iterator, e1);
    free_clast_expr(e2);
    *old_next = &ass->stmt;
    if (guard)
	ass->stmt.next = guard;
    else
	*next = &(**next)->next;
  } else {
    free_clast_expr(e1);
    free_clast_expr(e2);
  }

  return 1;    
}


/**
 * insert_block function:
 * This function inserts a statement block.
 * - block is the statement block,
 * - level is the number of loops enclosing the statement,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in domain (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - September 21th 2003: first version (pick from pprint function). 
 */
static void insert_block(CloogBlock *block, int level,
			 struct clast_stmt ***next, CloogInfos *infos)
{
    CloogStatement * statement ;
    struct clast_stmt *subs;
   
    if (!block)
	return;

    for (statement = block->statement; statement; statement = statement->next) {
	CloogStatement *s_next = statement->next;

	subs = clast_equal(level,infos);

	statement->next = NULL;
	**next = &new_clast_user_stmt(statement, subs)->stmt;
	statement->next = s_next;
	*next = &(**next)->next;
    }
}


/**
 * insert_loop function:
 * This function converts the content of a CloogLoop structure (loop) into a
 * clast_stmt (inserted at **next).
 * The iterator (level) of
 * the current loop is given by 'level': this is the column number of the
 * domain corresponding to the current loop iterator. The data of a loop are
 * written in this order:
 * 1. The guard of the loop, i.e. each constraint in the domain that does not
 *    depend on the iterator (when the entry in the column 'level' is 0).
 * 2. The iteration domain of the iterator, given by the constraints in the
 *    domain depending on the iterator, i.e.:
 *    * an equality if the iterator has only one value (possibly preceded by
 *      a guard verifying if this value is integral), *OR*
 *    * a loop from the minimum possible value of the iterator to the maximum
 *      possible value.
 * 3. The included statement block.
 * 4. The inner loops (recursive call).
 * 5. The following loops (recursive call).
 * - level is the recursion level or the iteration level that we are printing,
 * - the infos structure gives the user some options about code printing,
 *   the number of parameters in domain (nb_par), and the arrays of iterator
 *   names and parameters (iters and params). 
 **
 * - November   2nd 2001: first version. 
 * - March      6th 2003: infinite domain support. 
 * - April     19th 2003: (debug) NULL loop support. 
 * - June      29th 2003: non-unit strides support.
 * - April     28th 2005: (debug) level is level+equality when print statement!
 * - June      16th 2005: (debug) the N. Vasilache normalization step has been
 *                        added to avoid iteration duplication (see DaeGon Kim
 *                        bug in cloog_program_generate). Try vasilache.cloog
 *                        with and without the call to cloog_matrix_normalize,
 *                        using -f 8 -l 9 options for an idea.
 * - September 15th 2005: (debug) don't close equality braces when unnecessary.
 * - October   16th 2005: (debug) scalar value is saved for next loops.
 */
static void insert_loop(CloogLoop * loop, int level, int scalar,
			struct clast_stmt ***next, CloogInfos *infos)
{
    int equality=0, scalar_level;
    CloogConstraintSet *constraints, *temp;
    struct clast_stmt **top = *next;
    CloogConstraint i, j;
    int empty_loop = 0;

    /* It can happen that loop be NULL when an input polyhedron is empty. */
    if (loop == NULL)
	return;

    /* The constraints do not always have a shape that allows us to generate code from it,
    * thus we normalize it, we also simplify it with the equalities.
    */ 
    temp = cloog_domain_constraints(loop->domain);
    cloog_constraint_set_normalize(temp,level);
    constraints = cloog_constraint_set_simplify(temp,infos->equal,level,
				   infos->names->nb_parameters);
    cloog_constraint_set_free(temp);
    if (level)
      cloog_int_set(infos->stride[level-1], loop->stride);

    /* First of all we have to print the guard. */
    insert_guard(constraints,level, next, infos);

    if (level && cloog_constraint_set_contains_level(constraints, level,
					infos->names->nb_parameters)) {
	/* We scan all the constraints to know in which case we are :
	 * [[if] equation] or [for].
	 */
	if (cloog_constraint_is_valid(i =
		cloog_constraint_set_defining_equality(constraints, level))) {
	  insert_equation(i, cloog_constraint_invalid(), level, next, infos);
	  equality = 1 ;   
	} else if (cloog_constraint_is_valid(i =
		    cloog_constraint_set_defining_inequalities(constraints,
			      level, &j, infos->names->nb_parameters))) {
	    insert_equation(i, j, level, next, infos);
	} else
	    empty_loop = !insert_for(constraints, level, next, infos);
    }

    if (!empty_loop) {
	/* Finally, if there is an included statement block, print it. */
	insert_block(loop->block, level+equality, next, infos);

	/* Go to the next level. */
	if (loop->inner != NULL)
	    insert_loop(loop->inner, level+1,scalar, next, infos);
    }

    if (level)
      cloog_equal_del(infos->equal,level);
    cloog_constraint_set_free(constraints);

    /* Go to the next loop on the same level. */
    while (*top)
	top = &(*top)->next;
    if (loop->next != NULL)
	insert_loop(loop->next, level,scalar_level, &top,infos);
}


struct clast_stmt *cloog_clast_create(CloogProgram *program,
				      CloogOptions *options)
{
    CloogInfos *infos = ALLOC(CloogInfos);
    int i, nb_levels;
    struct clast_stmt *root = &new_clast_root(program->names)->stmt;
    struct clast_stmt **next = &root->next;

    infos->state      = options->state;
    infos->names    = program->names;
    infos->options  = options;
    infos->scaldims = program->scaldims;
    infos->nb_scattdims = program->nb_scattdims;

    /* Allocation for the array of strides, there is a +1 since the statement can
    * be included inside an external loop without iteration domain.
    */ 
    nb_levels = program->names->nb_scattering+program->names->nb_iterators+1;
    infos->stride = ALLOCN(cloog_int_t, nb_levels);
    for (i = 0; i < nb_levels; ++i)
	cloog_int_init(infos->stride[i]);

    infos->equal = cloog_equal_alloc(nb_levels,
			       nb_levels, program->names->nb_parameters);
	
    insert_loop(program->loop, 0, 0, &next, infos);

    cloog_equal_free(infos->equal);

    for (i = 0; i < nb_levels; ++i)
	cloog_int_clear(infos->stride[i]);
    free(infos->stride);
    free(infos);

    return root;
}
