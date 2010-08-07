/** 
 * $Id: matrix_permutations.c,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 *
 * Permutations on matrices Matrices are seen either as transformations
 * (mtransformation) or as polyhedra (mpolyhedron)
 * @author B. Meister
 * LSIIT -ICPS 
 * UMR 7005 CNRS
 * Louis Pasteur University (ULP), Strasbourg, France 
 * 
 * Permutations are just indirection vectors: the k^th element of a permutation
 * vector is the position of the k^th variable in the permuted object.
 */

#include <stdlib.h>
#include <polylib/matrix_permutations.h>

/** utility function : bit count (i know, there are faster methods) */
unsigned int nb_bits(unsigned long long int x) {
  unsigned int i,n=0;
  unsigned long long int y=x;
  for (i=0; i< 64; i++) {
    n+=y%2;
    y>>=1;
  }
  return n;
}


/** Gives the inverse permutation vector of a permutation vector 
 * @param perm the permutation vector
 * @param 
*/
unsigned int * permutation_inverse(unsigned int * perm, unsigned int nb_elems) {
  int i;
  unsigned int * inv_perm = (unsigned int *)malloc(sizeof(unsigned int) * nb_elems);
  for (i=0; i< nb_elems; i++) inv_perm[perm[i]] = i;
  return inv_perm;
}
  

/**
 * Given a linear tranformation on initial variables, and a variable
 * permutation, computes the tranformation for the permuted variables.  perm is
 * a vector giving the new "position of the k^th variable, k \in [1..n] we can
 * call it a "permutation vector" if you wish transf[x][y] ->
 * permuted[permutation(x)][permutation(y)]
 */
Matrix * mtransformation_permute(Matrix * transf, unsigned int * permutation) {
  Matrix * permuted;
  unsigned int i,j;
  /* the transformation is supposed to be from Q^n to Q^n, so a square matrix. */
  assert(transf->NbRows==transf->NbColumns);
  permuted = Matrix_Alloc(transf->NbRows, transf->NbRows);
  for (i= 0; i< transf->NbRows; i++) {
    for (j= 0; j< transf->NbRows; j++) {
      value_assign(permuted->p[permutation[i]][permutation[j]], transf->p[i][j]);
    }
  }
  return permuted;
}


/** permutes the variables of the constraints of a polyhedron 
 * @param polyh the constraints of the polyhedron
 * @param permutation a permutation vector
*/
Matrix * mpolyhedron_permute(Matrix * polyh, unsigned int * permutation) {
  unsigned int i,j;
  Matrix * permuted = Matrix_Alloc(polyh->NbRows, polyh->NbColumns);
  for (i= 0; i< polyh->NbRows; i++) {
    value_assign(permuted->p[i][0], polyh->p[i][0]);
    for (j= 1; j< polyh->NbColumns; j++) {
      value_assign(permuted->p[i][permutation[j-1]+1], polyh->p[i][j]);
    }
  }
  return permuted;
}


/** permutes the variables of the constraints of a polyhedron 
 * @param C the original set of constraints
 * @param perm a permutation vector
 * @param Cp (returned) the set of constraints whose variables are
 * permuted. Allocated if set to NULL, assumed to be already allocated if not.
 */
void Constraints_permute(Matrix * C, unsigned int * perm, Matrix ** Cp) {
  unsigned int i,j;
  if ((*Cp)==NULL) {
    (*Cp) = Matrix_Alloc(C->NbRows, C->NbColumns);
  }
  else {
    assert((*Cp)->NbRows == C->NbRows && (*Cp)->NbColumns==C->NbColumns);
  }
  for (i= 0; i< C->NbRows; i++) {
    value_assign((*Cp)->p[i][0], C->p[i][0]);
    for (j= 1; j< C->NbColumns; j++) {
      value_assign((*Cp)->p[i][perm[j-1]+1], C->p[i][j]);
    }
  }
} /* Constraints_permute */


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
				      unsigned end) {
  unsigned long long int combination;
  unsigned int i,j,k;
  Matrix * M, * H, * Q, *U;
  Matrix * Square_Mat, *Eqs2;
  unsigned nb_vars = end - start + 1 ;
  Polyhedron * OverConstrained;

  assert (start>0 && end < Eqs->NbColumns-1);

  /* if the affine hull is overconstrained, return 0 */
  if (Eqs->NbRows >nb_vars) {
    /* FIXME: there is a magic maximum number of rays here */
    Eqs2 = Matrix_Copy(Eqs);
    OverConstrained = Constraints2Polyhedron(Eqs2,
					     Eqs->NbColumns*Eqs->NbColumns);
    Matrix_Free(Eqs2);
    if (emptyQ(OverConstrained)) {
      Polyhedron_Free(OverConstrained);
      return 0;
    }
    Polyhedron_Free(OverConstrained);
  }

  /* do not accept 0 = 0 equalities */
  for (i=0; i< Eqs->NbRows; i++) {
    assert (!Vector_IsZero(Eqs->p[i], Eqs->NbColumns));
  }
  
  Square_Mat= Matrix_Alloc(Eqs->NbRows, Eqs->NbRows);
  
  /* There are Eqs->NbRows variables to eliminate.
     Generate all the combinations of Eqs->NbRows variables (-> bits to 1 in
     the word "combination") among nb_vars WARNING : we assume here that we
     have not more than 64 variables.  You may convert it to use GNU MP to
     set it to an infinite number of bits 
  */
  for (combination = ((unsigned long long int) 1<<(Eqs->NbRows))-1;
       (combination < ((unsigned long long int) 1 << nb_vars)) ;
       combination++) {
    if (nb_bits(combination) == Eqs->NbRows) {
      k=0;
      /* 1- put the m colums in a square matrix */
      for (j=0; j< nb_vars; j++) {
	if ((combination>>j)%2) {
	  for (i=0; i< Eqs->NbRows; i++) {
	    value_assign(Square_Mat->p[i][k], Eqs->p[i][j+start]);
	  }
	  k++;
	}
      }
      /* 2- see if the matrix is full-row-rank */
      right_hermite(Square_Mat, &H, &Q, &U);
      Matrix_Free(Q);
      Matrix_Free(U);

      /* if it is full-row-rank, we have found a set of variables that can be
	 eliminated. */
      if ( value_notzero_p((H->p[Eqs->NbRows-1][Eqs->NbRows-1])) ) {
	Matrix_Free(Square_Mat);
	Matrix_Free(H);
	return combination;
      }
      Matrix_Free(H);
    }
  }
  Matrix_Free(Square_Mat);
  return (unsigned long long int) 0;
} /* eliminable_vars */



/** 
 * finds a valid permutation : for a set of m equations, find m variables that
 * will be put at the beginning (to be eliminated). 
 * Note: inherits the limited the number of variables from
 * <i>eliminable_vars</i>
 */
unsigned int * find_a_permutation(Matrix * Eqs, unsigned int nb_parms) {
  unsigned int i, j, k;
  int nb_vars = Eqs->NbColumns-nb_parms-2;
  unsigned long long int combination;
  unsigned int * permutation = (unsigned int *)malloc(sizeof(unsigned int) *
						      Eqs->NbColumns-1);

  /* 1- find a set of variables to eliminate */
  if ((combination = eliminable_vars(Eqs, 1, nb_vars)) == 0) {
    /* if it is impossible to eliminate enough variables, return error code */
    return NULL;
  }

  /* 2- make the permutation matrix
   *   a- deal with the variables */
  k=0;
  for (i=0; i< nb_vars; i++) {
    /* if the variable has to be eliminated, put them at the beginning */
    if (combination%2) {
      permutation[i] = k;
      k++;
    }
    /* if not, put the variables at the end */
    else permutation[i] = Eqs->NbRows+nb_parms+ i-k;
    combination>>=1;
  }
  /*  b- deal with the parameters */
  for (i=0; i< nb_parms; i++) {
    permutation[nb_vars+i] = Eqs->NbRows+i;
  }
  /*  c- deal with the constant */
  permutation[Eqs->NbColumns-2] = Eqs->NbColumns-2;
  
  return permutation;
} /* find_a_permutation */



/** computes the permutation of variables and parameters, according to some
 * variables to keep.  put the variables not to be kept at the beginning, then
 * the parameters and finally the variables to be kept.  strongly related to
 * the function compress_to_full_dim2
 */
unsigned int * permutation_for_full_dim2(unsigned int * vars_to_keep, 
					 unsigned int nb_keep, 
					 unsigned int nb_vars_parms, 
					 unsigned int nb_parms) {
  unsigned int * permutation = 
    (unsigned int*)malloc(sizeof(unsigned int) * nb_vars_parms+1);
  unsigned int i;
  int cur_keep =0, cur_go = 0;/*current number of variables to eliminate/keep*/
  for (i=0; i< nb_vars_parms - nb_parms; i++) {
    if (i==vars_to_keep[cur_keep]) {
      permutation[i] = nb_vars_parms-nb_keep+cur_keep;
      cur_keep++;
    }
    else {
      permutation[i] = cur_go;
      cur_go++;
    }
  }
  /* parameters are just left-shifted */
  for (i=0; i< nb_parms; i++)
    permutation[i+nb_vars_parms-nb_parms] = i+nb_vars_parms-nb_parms-nb_keep;

  /* contants stay where they are */
  permutation[nb_vars_parms] = nb_vars_parms;
  return permutation;
} /* permutation_for_full_dim2 */
