/**
 * $Id: compress_parms.h,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 * @author B. Meister 12/2003-2006
 * LSIIT -ICPS 
 * UMR 7005 CNRS
 * Louis Pasteur University (ULP), Strasbourg, France
 */
#ifndef __BM_COMPRESS_PARMS_H__
#define __BM_COMPRESS_PARMS_H__

#include "matrix_addon.h"
#include "matrix_permutations.h"
#include <assert.h>


/* ----- functions applying on equalities ----- */

/** 
 * Given a system of non-redundant equalities, looks if it has an integer
 * solution in the combined space, and if yes, returns one solution.
 */
void Equalities_integerSolution(Matrix * Eqs, Matrix ** sol);

/** 
 * Computes the validity lattice of a set of equalities. I.e., the lattice
 * induced on the last <tt>b</tt> variables by the equalities involving the
 * first <tt>a</tt> integer existential variables.
 */
void Equalities_validityLattice(Matrix * Eqs, int a, Matrix** vl);

/** 
 * Given an integer matrix B with m rows and integer m-vectors C and d,
 * computes the basis of the integer solutions to (BN+C) mod d = 0 (1).
 * This is an affine lattice (G): (N 1)^T= G(N' 1)^T, forall N' in Z^b.
 * If there is no solution, returns NULL.
*/
void Equalities_intModBasis(Matrix * B, Matrix * C, Matrix * d, Matrix ** imb);


/* ----- functions applying on constraints ----- */


/**
 * Eliminates all the equalities in a set of constraints and returns the set of
 * constraints defining a full-dimensional polyhedron, such that there is a
 * bijection between integer points of the original polyhedron and these of the
 * resulting (projected) polyhedron).
 */
void Constraints_fullDimensionize(Matrix ** M, Matrix ** C, Matrix ** VL, 
				  Matrix ** Eqs, Matrix ** ParmEqs, 
				  unsigned int ** elimVars, 
				  unsigned int ** elimParms,
				  int maxRays);

/* extracts equalities involving only parameters */
#define Constraints_removeParmEqs(a,b,c,d) Constraints_Remove_parm_eqs(a,b,c,d)
Matrix * Constraints_Remove_parm_eqs(Matrix ** M, Matrix ** Ctxt, 
				     int renderSpace, 
				     unsigned int ** elimParms);

/**
 * Eliminates the columns corresponding to a list of eliminated parameters.
 */
void Constraints_removeElimCols(Matrix * M, unsigned int nbVars, 
				unsigned int *elimParms, Matrix ** newM);


/* ----- function applying on a lattice ----- */

/**
 * Given a matrix that defines a full-dimensional affine lattice, returns the 
 * affine sub-lattice spanned in the k first dimensions.
 * Useful for instance when you only look for the parameters' validity lattice.
 */
void Lattice_extractSubLattice(Matrix * lat, unsigned int k, Matrix ** subLat);


/* ----- functions applying on a polyhedron ----- */


Polyhedron * Polyhedron_Remove_parm_eqs(Polyhedron ** P, Polyhedron ** C, 
				      int renderSpace, 
				      unsigned int ** elimParms, 
				      int maxRays);
#define Polyhedron_removeParmEqs(a,b,c,d,e) Polyhedron_Remove_parm_eqs(a,b,c,d,e)


/* ----- functions kept for backwards compatibility ----- */


/** 
 * given a full-row-rank nxm matrix M(made of row-vectors), 
 * computes the basis K (made of n-m column-vectors) of the integer kernel of M
 * so we have: M.K = 0
*/
Matrix * int_ker(Matrix * M);

/* given a matrix of m parameterized equations, compress the parameters and
 transform the variable space into a n-m space. */
Matrix * full_dimensionize(Matrix const * M, int nb_parms, 
			   Matrix ** Validity_Lattice);

/* Compute the overall period of the variables I for (MI) mod |d|,
 where M is a matrix and |d| a vector
 Produce a diagonal matrix S = (s_k) where s_k is the overall period of i_k */
Matrix * affine_periods(Matrix * M, Matrix * d);

/* given a matrix B' with m rows and m-vectors C' and d, computes the 
 basis of the integer solutions to (B'N+C') mod d = 0.
returns NULL if there is no integer solution */
Matrix * int_mod_basis(Matrix * Bp, Matrix * Cp, Matrix * d);

/* given a parameterized constraints matrix with m equalities, computes the
 compression matrix C such that there is an integer solution in the variables
 space for each value of N', with N = Cmp N' (N are the original parameters) */
Matrix * compress_parms(Matrix * E, int nb_parms);


#endif /* __BM_COMPRESS_PARMS_H__ */
