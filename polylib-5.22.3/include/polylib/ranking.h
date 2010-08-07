/**
 * Tools to compute the ranking function of an iteration J: the number of
 * integer points in P that are lexicographically inferior to J 
 * @author B. Meister <meister@icps.u-strasbg.fr>
 * 6/2005
 * LSIIT-ICPS, UMR 7005 CNRS Université Louis Pasteur
 * HiPEAC Network
 */

#ifndef __BM_POLYLIB_RANKING_H__
#define __BM_POLYLIB_RANKING_H__
#include <polylib/polylib.h>

/*
 * Returns a list of polytopes needed to compute
 * the number of points in P that are lexicographically
 * smaller than a given point in D.
 * Only the first dim dimensions are taken into account
 * for computing the lexsmaller relation.
 * The remaining variables are assumed to be extra
 * existential/control variables.
 * When P == D, this is the conventional ranking function.
 * P and D are assumed to have the same parameter domain C.
 *
 * The first polyhedron in the list returned is the
 * updated context: a combination of D and C or an extended C.
 *
 * The order of the variables in the remaining polyhedra is
 * - first dim variables of P
 * - existential variables of P
 * - existential variables of D
 * - first dim variables of D
 * - the parameters
 */
Polyhedron *LexSmaller(Polyhedron *P, Polyhedron *D, unsigned dim,
			Polyhedron *C, unsigned MAXRAYS);

/*
 * Returns the number of points in P that are lexicographically
 * smaller than a given point in D.
 * Only the first dim dimensions are taken into account
 * for computing the lexsmaller relation.
 * The remaining variables are assumed to be extra
 * existential/control variables.
 * When P == D, this is the conventional ranking function.
 * P and D are assumed to have the same parameter domain C.
 * The variables in the Enumeration correspond to the first dim variables
 * in D followed by the parameters of D (the variables of C).
 */
Enumeration *Polyhedron_LexSmallerEnumerate(Polyhedron *P, Polyhedron *D, 
					    unsigned dim,
					    Polyhedron *C, unsigned MAXRAYS);

/*
 * Returns a function that assigns a unique number to each point in the
 * polytope P ranging from zero to (number of points in P)-1.
 * The order of the numbers corresponds to the lexicographical order.
 *
 * C is the parameter context of the polytope
 */
Enumeration *Polyhedron_Ranking(Polyhedron *P, Polyhedron *C, unsigned MAXRAYS);

#endif /* __BM_POLYLIB_RANKING_H__ */
