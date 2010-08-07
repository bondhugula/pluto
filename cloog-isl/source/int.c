#include <stdlib.h>
#include "../include/cloog/cloog.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n)*sizeof(type))

#if defined(CLOOG_INT_INT) || \
    defined(CLOOG_INT_LONG) || \
    defined(CLOOG_INT_LONG_LONG)

cloog_int_t cloog_gcd(cloog_int_t a, cloog_int_t b)
{
	while (a) {
		cloog_int_t t = b % a;
		b = a;
		a = t;
	}
	if (b < 0)
		b = -b;
	return b;
}

#endif

struct cloog_vec *cloog_vec_alloc(unsigned size)
{
	int i;
	struct cloog_vec *vec;

	vec = ALLOC(struct cloog_vec);
	if (!vec)
		return NULL;

	vec->p = ALLOCN(cloog_int_t, size);
	if (!vec->p)
		goto error;
	vec->size = size;

	for (i = 0; i < size; ++i)
		cloog_int_init(vec->p[i]);

	return vec;
error:
	free(vec);
	return NULL;
}

void cloog_vec_free(struct cloog_vec *vec)
{
	int i;

	if (!vec)
		return;

	for (i = 0; i < vec->size; ++i)
		cloog_int_clear(vec->p[i]);
	free(vec->p);
	free(vec);
}

void cloog_seq_neg(cloog_int_t *dst, cloog_int_t *src, unsigned len)
{
	int i;
	for (i = 0; i < len; ++i)
		cloog_int_neg(dst[i], src[i]);
}

void cloog_seq_combine(cloog_int_t *dst, cloog_int_t m1, cloog_int_t *src1,
			cloog_int_t m2, cloog_int_t *src2, unsigned len)
{
	int i;
	cloog_int_t tmp;

	cloog_int_init(tmp);
	for (i = 0; i < len; ++i) {
		cloog_int_mul(tmp, m1, src1[i]);
		cloog_int_addmul(tmp, m2, src2[i]);
		cloog_int_set(dst[i], tmp);
	}
	cloog_int_clear(tmp);
}
