#ifndef ISL_CONSTRAINT_H
#define ISL_CONSTRAINT_H

#include "isl_div.h"
#include "isl_set.h"

#if defined(__cplusplus)
extern "C" {
#endif

struct isl_constraint {
	int ref;
	struct isl_ctx *ctx;

	struct isl_basic_map	*bmap;
	isl_int			**line;
};

struct isl_constraint *isl_equality_alloc(struct isl_dim *dim);
struct isl_constraint *isl_inequality_alloc(struct isl_dim *dim);
struct isl_constraint *isl_basic_set_constraint(struct isl_basic_set *bset,
	isl_int **line);

struct isl_constraint *isl_constraint_cow(struct isl_constraint *c);
struct isl_constraint *isl_constraint_copy(struct isl_constraint *c);
struct isl_constraint *isl_constraint_free(struct isl_constraint *c);

struct isl_constraint *isl_basic_set_first_constraint(
	struct isl_basic_set *bset);
struct isl_constraint *isl_constraint_next(struct isl_constraint *c);
int isl_constraint_is_equal(struct isl_constraint *constraint1,
			    struct isl_constraint *constraint2);

struct isl_basic_map *isl_basic_map_add_constraint(
	struct isl_basic_map *bmap, struct isl_constraint *constraint);
struct isl_basic_set *isl_basic_set_add_constraint(
	struct isl_basic_set *bset, struct isl_constraint *constraint);

int isl_basic_set_has_defining_equality(
	struct isl_basic_set *bset, enum isl_dim_type type, int pos,
	struct isl_constraint **constraint);
int isl_basic_set_has_defining_inequalities(
	struct isl_basic_set *bset, enum isl_dim_type type, int pos,
	struct isl_constraint **lower,
	struct isl_constraint **upper);

int isl_constraint_dim(struct isl_constraint *constraint,
	enum isl_dim_type type);

void isl_constraint_get_constant(struct isl_constraint *constraint, isl_int *v);
void isl_constraint_get_coefficient(struct isl_constraint *constraint,
	enum isl_dim_type type, int pos, isl_int *v);
void isl_constraint_set_constant(struct isl_constraint *constraint, isl_int v);
void isl_constraint_set_coefficient(struct isl_constraint *constraint,
	enum isl_dim_type type, int pos, isl_int v);

struct isl_div *isl_constraint_div(struct isl_constraint *constraint, int pos);
struct isl_constraint *isl_constraint_add_div(struct isl_constraint *constraint,
	struct isl_div *div, int *pos);

void isl_constraint_clear(struct isl_constraint *constraint);
struct isl_constraint *isl_constraint_negate(struct isl_constraint *constraint);

int isl_constraint_is_equality(struct isl_constraint *constraint);

struct isl_basic_set *isl_basic_set_from_constraint(
	struct isl_constraint *constraint);

#if defined(__cplusplus)
}
#endif

#endif
