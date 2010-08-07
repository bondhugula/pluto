/** 
 * $Id: matrix_addon.h,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 * 
 * Polylib matrix addons
 * Mainly, deals with polyhedra represented in implicit form (set of
 * constraints).
 * @author Benoit Meister 
 */

#ifndef __BM_MATRIX_ADDON_H__
#define __BM_MATRIX_ADDON_H__

#include<polylib/polylib.h>
#include<assert.h>

/** Shortcut for Matrix_Print */
#define show_matrix(M) { printf(#M"= \n"); \
                         if (M!=NULL) { \
			 Matrix_Print(stderr,P_VALUE_FMT,(M));} \
                         else {printf("<NULL>\n");} \
                       } 

/** 
 * Allocates a matrix if it is null, or else asserts that it has at least a
 * certain size */
#define ensureMatrix(M, r, c) { if (M==NULL) M = Matrix_Alloc(r,c); \
                                else assert (M->NbRows>=r && M->NbColumns>=c); \
                              } 

/* Creates a view of the constraints of a polyhedron as a Matrix * */
Matrix * constraintsView(Polyhedron * P);

/* "Frees" a view of the constraints of a polyhedron */
void constraintsView_Free(Matrix * M);

/* splits a matrix of constraints M into a matrix of equalities Eqs and a
   matrix of inequalities Ineqs allocs the new matrices. */
void split_constraints(Matrix const * M, Matrix ** Eqs, Matrix **Ineqs);

/* returns the dim-dimensional identity matrix */
Matrix * Identity_Matrix(unsigned int dim);

/* given a n x n integer transformation matrix transf, compute its inverse M/g,
 where M is a nxn integer matrix.  g is a common denominator for elements of
 (transf^{-1})*/
void mtransformation_inverse(Matrix * transf, Matrix ** inv, Value * g);

/* simplifies a matrix seen as a polyhedron, by dividing its rows by the gcd of
their elements. */
void mpolyhedron_simplify(Matrix * polyh);

/* inflates a polyhedron (represented as a matrix) P, so that the apx of its
   Ehrhart Polynomial is an upper bound of the Ehrhart polynomial of P WARNING:
   this inflation is supposed to be applied on full-dimensional polyhedra. */
void mpolyhedron_inflate(Matrix * polyh, unsigned int nb_parms);

/* deflates a polyhedron (represented as a matrix) P, so that the apx of its
   Ehrhart Polynomial is a lower bound of the Ehrhart polynomial of P WARNING:
   this deflation is supposed to be applied on full-dimensional polyhedra. */
void mpolyhedron_deflate(Matrix * polyh, unsigned int nb_parms);

/* use an eliminator row to eliminate a variable in a victim row (without
changing the sign of the victim row -> important if it is an inequality).  */
void eliminate_var_with_constr(Matrix * Eliminator, 
			       unsigned int eliminator_row, Matrix * Victim, 
			       unsigned int victim_row, 
			       unsigned int var_to_elim);


/* ----- PARTIAL MAPPINGS ----- */

/* compresses the last vars/pars of the polyhedron M expressed as a polylib
   matrix
 - adresses the full-rank compressions only
 - modfies M */
void mpolyhedron_compress_last_vars(Matrix * M, Matrix * compression);
#define Constraints_compressLastVars(a, b) mpolyhedron_compress_last_vars(a, b)

/* uses a set of m equalities Eqs to eliminate m variables in the polyhedron.
 Ineqs represented as a matrix eliminates the m first variables 
- assumes that Eqs allows to eliminate the m equalities 
- modifies Ineqs */
unsigned int mpolyhedron_eliminate_first_variables(Matrix * Eqs, 
						   Matrix * Ineqs);
#define Constraints_eliminateFirstVars(a,b) mpolyhedron_eliminate_first_variables(a,b)

/** returns a contiguous submatrix of a matrix. */
void Matrix_subMatrix(Matrix * M, unsigned int sr, unsigned int sc, 
		      unsigned int nbR, unsigned int nbC, Matrix ** sub);

/**
 * Copies a contiguous submatrix of M1 into M2, at the indicated position.
 * M1 and M2 are assumed t be allocated already.
 */
void Matrix_copySubMatrix(Matrix *M1,
			  unsigned int sr1, unsigned int sc1,
			  unsigned int nbR, unsigned int nbC,
			  Matrix * M2,
			  unsigned int sr2, unsigned int sc2);

/** 
 * given a matrix M into -M
 */
void Matrix_oppose(Matrix * M);

#endif /* __BM_MATRIX_ADDON_H__ */
