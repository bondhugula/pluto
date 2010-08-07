/**
 * $Id: matrix_permutations.h,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 * Permutations on matrices
 * Matrices are seen either as transformations (mtransformation) or as
 * polyhedra (mpolyhedron or Constraints).
 * @author B. Meister
 */

#ifndef __BM_MATRIX_PERMUTATIONS_H__
#define __BM_MATRIX_PERMUTATIONS_H__

#include<polylib/polylib.h>
#include<assert.h>

/* Permutations here are vectors that give the future position for each
   variable */

/* utility function : bit count */
unsigned int nb_bits(unsigned long long int x);

/* gives the inverse permutation vector of a permutation vector */
unsigned int * permutation_inverse(unsigned int * perm, unsigned int nb_elems);

/*
 * Given a linear tranformation on initial variables, and a variable
 * permutation, compute the tranformation for the permuted variables.  perm is
 * a vector giving the new "position of the k^th variable, k \in [1..n] we can
 * call it a "permutation vector" if you wish transf[x][y] ->
 * permuted[permutation(x)][permutation(y)]
 */
Matrix * mtransformation_permute(Matrix * transf, unsigned int * permutation);

/* permutes the variables of a matrix seen as a polyhedron */
Matrix * mpolyhedron_permute(Matrix * polyh, unsigned int * permutation);

/* permutes the variables of a matrix seen as a polyhedron */
void Constraints_permute(Matrix * C, unsigned int * perm, Matrix ** Cp);

/** Given a set of <i>equalities</i>, find a set of variables that can be
 * eliminated using these equalities.  The variables that we agree to eliminate
 * are in a zone of contiguous variables (or parameters).  <p>
 * Notes: 
 <ul>
 <li>brute force, surely enhanceable algorithm</li>
 <li>limited number of variables in the zone: limit = bitwidth of long long
 </ul>
 * @param Eqs the matrix of equalities.
 * @param start the rank of the first variable (inclusive) of the zone in Eqs
 * @param end the rank of the last variable (inclusive) of the zone
 * return a bitfield where bits set to one define the variables to eliminate
*/
unsigned long long int eliminable_vars(Matrix * Eqs, unsigned start, 
				       unsigned end);

/*
* find a valid permutation : for a set of m equations, find m variables that
* will be put at the beginning (to be eliminated) it must be possible to
* eliminate these variables : the submatrix built with their columns must be
* full-rank.  brute force method, that tests all the combinations until finding
* one which works.  
* <b>LIMITATIONS</b> : up to x-1 variables, where the long long
* format is x-1 bits (often 64 in year 2005).  
*/
unsigned int * find_a_permutation(Matrix * Eqs, unsigned int nb_parms);

/*
* compute the permutation of variables and parameters, according to some
* variables to keep.  put the variables not to be kept at the beginning, then
* the parameters and finally the variables to be kept.  strongly related to the
* function compress_to_full_dim2
*/
unsigned int * permutation_for_full_dim2(unsigned int * vars_to_keep, 
					 unsigned int nb_keep, 
					 unsigned int nb_vars_parms, 
					 unsigned int nb_parms);

#endif /*__BM_MATRIX_PERMUTATIONS_H__ */
