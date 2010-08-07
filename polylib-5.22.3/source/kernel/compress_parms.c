/** 
 * $Id: compress_parms.c,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 *
 * The integer points in a parametric linear subspace of Q^n are generally
 * lying on a sub-lattice of Z^n.  
 * Functions here use and compute validity lattices, i.e. lattices induced on a
 * set of variables by such equalities involving another set of integer
 * variables.
 * @author B. Meister 12/2003-2006 meister@icps.u-strasbg.fr
 * LSIIT -ICPS 
 * UMR 7005 CNRS
 * Louis Pasteur University (ULP), Strasbourg, France 
*/

#include <stdlib.h>
#include <polylib/polylib.h>

/** 
 * debug flags (2 levels)
 */
#define dbgCompParm 0
#define dbgCompParmMore 0

#define dbgStart(a) if (dbgCompParmMore) { printf(" -- begin "); \
                                           printf(#a);        \
					   printf(" --\n"); }   \
                                           while(0)
#define dbgEnd(a) if (dbgCompParmMore) { printf(" -- end "); \
                                         printf(#a);      \
					 printf(" --\n"); } \
                                         while(0)


/** 
 * Given a full-row-rank nxm matrix M made of m row-vectors), computes the
 * basis K (made of n-m column-vectors) of the integer kernel of the rows of M
 * so we have: M.K = 0
*/
Matrix * int_ker(Matrix * M) {
  Matrix *U, *Q, *H, *H2, *K=NULL;
  int i, j, rk;

  if (dbgCompParm)
    show_matrix(M);
  /* eliminate redundant rows : UM = H*/
  right_hermite(M, &H, &Q, &U);
  for (rk=H->NbRows-1; (rk>=0) && Vector_IsZero(H->p[rk], H->NbColumns); rk--);
  rk++;
  if (dbgCompParmMore) {
    printf("rank = %d\n", rk);
  }
    
  /* there is a non-null kernel if and only if the dimension m of 
     the space spanned by the rows 
     is inferior to the number n of variables */
  if (M->NbColumns <= rk) {
    Matrix_Free(H);
    Matrix_Free(Q);
    Matrix_Free(U);
    K = Matrix_Alloc(M->NbColumns, 0);
    return K;
  }
  Matrix_Free(U); 
  Matrix_Free(Q);
  /* fool left_hermite  by giving NbRows =rank of M*/
  H->NbRows=rk;
  /* computes MU = [H 0] */
  left_hermite(H, &H2, &Q, &U); 
   if (dbgCompParmMore) {
    printf("-- Int. Kernel -- \n");
    show_matrix(M);
    printf(" = \n");
    show_matrix(H2);
    show_matrix(U); 
  }
  H->NbRows==M->NbRows;
  Matrix_Free(H);
  /* the Integer Kernel is made of the last n-rk columns of U */
  Matrix_subMatrix(U, 0, rk, U->NbRows, U->NbColumns, &K);

  /* clean up */
  Matrix_Free(H2);
  Matrix_Free(U);
  Matrix_Free(Q);
  return K;
} /* int_ker */


/** 
 * Computes the intersection of two linear lattices, whose base vectors are
 * respectively represented in A and B.
 * If I and/or Lb is set to NULL, then the matrix is allocated. 
 * Else, the matrix is assumed to be allocated already. 
 * I and Lb are rk x rk, where rk is the rank of A (or B).
 * @param A the full-row rank matrix whose column-vectors are the basis for the
 * first linear lattice.
 * @param B the matrix whose column-vectors are the basis for the second linear
 * lattice.
 * @param Lb the matrix such that B.Lb = I, where I is the intersection.
 * @return their intersection.
 */
static void linearInter(Matrix * A, Matrix * B, Matrix ** I, Matrix **Lb) {
  Matrix * AB=NULL;
  int rk = A->NbRows;
  int a = A->NbColumns;
  int b = B->NbColumns;
  int i,j, z=0;

  Matrix * H, *U, *Q;
  /* ensure that the spanning vectors are in the same space */
  assert(B->NbRows==rk);
  /* 1- build the matrix 
   * (A 0 1)
   * (0 B 1)
   */
  AB = Matrix_Alloc(2*rk, a+b+rk);
  Matrix_copySubMatrix(A, 0, 0, rk, a, AB, 0, 0);
  Matrix_copySubMatrix(B, 0, 0, rk, b, AB, rk, a);
  for (i=0; i< rk; i++) {
      value_set_si(AB->p[i][a+b+i], 1);
      value_set_si(AB->p[i+rk][a+b+i], 1);
  }
  if (dbgCompParm) {
    show_matrix(AB);
  }

  /* 2- Compute its left Hermite normal form. AB.U = [H 0] */
  left_hermite(AB, &H, &Q, &U);
  Matrix_Free(AB);
  Matrix_Free(Q);
  /* count the number of non-zero colums in H */ 
  for (z=H->NbColumns-1; value_zero_p(H->p[H->NbRows-1][z]); z--);
  z++;
  if (dbgCompParm) {
    show_matrix(H);
    printf("z=%d\n", z);
  }
  Matrix_Free(H);
  /* if you split U in 9 submatrices, you have: 
   * A.U_13 = -U_33
   * B.U_23 = -U_33,
   * where the nb of cols of U_{*3} equals the nb of zero-cols of H
   * U_33 is a (the smallest) combination of col-vectors of A and B at the same
   * time: their intersection.
  */
  Matrix_subMatrix(U, a+b, z, U->NbColumns, U->NbColumns, I);
  Matrix_subMatrix(U, a, z, a+b, U->NbColumns, Lb);
  if (dbgCompParm) {
    show_matrix(U);
  }
  Matrix_Free(U);
} /* linearInter */


/** 
 * Given a system of equalities, looks if it has an integer solution in the
 * combined space, and if yes, returns one solution.
 * <p>pre-condition: the equalities are full-row rank (without the constant
 * part)</p>
 * @param Eqs the system of equations (as constraints)
 * @param I a feasible integer solution if it exists, else NULL. Allocated if
 * initially set to NULL, else reused.
 */
void Equalities_integerSolution(Matrix * Eqs, Matrix **I) {
  Matrix * Hm, *H=NULL, *U, *Q, *M=NULL, *C=NULL, *Hi;
  Matrix *Ip;
  int i;
  Value mod;
  unsigned int rk;
  if (Eqs==NULL){
    if ((*I)!=NULL) Matrix_Free(*I);
    I = NULL;
    return;
  }
  /* we use: AI = C = (Ha 0).Q.I = (Ha 0)(I' 0)^T */
  /* with I = Qinv.I' = U.I'*/
  /* 1- compute I' = Hainv.(-C) */
  /* HYP: the equalities are full-row rank */
  rk = Eqs->NbRows;
  Matrix_subMatrix(Eqs, 0, 1, rk, Eqs->NbColumns-1, &M);
  left_hermite(M, &Hm, &Q, &U);
  Matrix_Free(M);
  Matrix_subMatrix(Hm, 0, 0, rk, rk, &H);
  if (dbgCompParmMore) {
    show_matrix(Hm);
    show_matrix(H);
    show_matrix(U);
  }
  Matrix_Free(Q);
  Matrix_Free(Hm);
  Matrix_subMatrix(Eqs, 0, Eqs->NbColumns-1, rk, Eqs->NbColumns, &C);
  Matrix_oppose(C);
  Hi = Matrix_Alloc(rk, rk+1);
  MatInverse(H, Hi);
  if (dbgCompParmMore) {
    show_matrix(C);
    show_matrix(Hi);
  }
  /* put the numerator of Hinv back into H */
  Matrix_subMatrix(Hi, 0, 0, rk, rk, &H);
  Ip = Matrix_Alloc(Eqs->NbColumns-2, 1);
  /* fool Matrix_Product on the size of Ip */
  Ip->NbRows = rk;
  Matrix_Product(H, C, Ip);
  Ip->NbRows = Eqs->NbColumns-2;
  Matrix_Free(H);
  Matrix_Free(C);
  value_init(mod);
  for (i=0; i< rk; i++) {
    /* if Hinv.C is not integer, return NULL (no solution) */
    value_pmodulus(mod, Ip->p[i][0], Hi->p[i][rk]);
    if (value_notzero_p(mod)) { 
      if ((*I)!=NULL) Matrix_Free(*I);
      value_clear(mod);
      Matrix_Free(U);
      Matrix_Free(Ip);
      Matrix_Free(Hi);
      I = NULL;
      return;
    }
    else {
      value_pdivision(Ip->p[i][0], Ip->p[i][0], Hi->p[i][rk]);
    }
  }
  /* fill the rest of I' with zeros */
  for (i=rk; i< Eqs->NbColumns-2; i++) {
    value_set_si(Ip->p[i][0], 0);
  }
  value_clear(mod);
  Matrix_Free(Hi);
  /* 2 - Compute the particular solution I = U.(I' 0) */
  ensureMatrix((*I), Eqs->NbColumns-2, 1);
  Matrix_Product(U, Ip, (*I));
  Matrix_Free(U);
  Matrix_Free(Ip);
  if (dbgCompParm) {
    show_matrix(*I);
  }
}


/** 
 * Computes the validity lattice of a set of equalities. I.e., the lattice
 * induced on the last <tt>b</tt> variables by the equalities involving the
 * first <tt>a</tt> integer existential variables.  The submatrix of Eqs that
 * concerns only the existential variables (so the first a columns) is assumed
 * to be full-row rank.
 * @param Eqs the equalities
 * @param a the number of existential integer variables, placed as first
 * variables
 * @param vl the (returned) validity lattice, in homogeneous form. It is
 * allocated if initially set to null, or reused if already allocated.
 */
void Equalities_validityLattice(Matrix * Eqs, int a, Matrix** vl) {
  unsigned int b = Eqs->NbColumns-2-a;
  unsigned int r = Eqs->NbRows;
  Matrix * A=NULL, * B=NULL, *I = NULL, *Lb=NULL, *sol=NULL;
  Matrix *H, *U, *Q;
  unsigned int i;

  if (dbgCompParm) {
    printf("Computing validity lattice induced by the %d first variables of:"
	   ,a);
    show_matrix(Eqs);
  }
  if (b==0) {
    ensureMatrix((*vl), 1, 1);
    value_set_si((*vl)->p[0][0], 1);
    return;
  }

  /* 1- check that there is an integer solution to the equalities */
  /* OPT: could change integerSolution's profile to allocate or not*/
  Equalities_integerSolution(Eqs, &sol);
  /* if there is no integer solution, there is no validity lattice */
  if (sol==NULL) {
    if ((*vl)!=NULL) Matrix_Free(*vl);
    return;
  }
  Matrix_subMatrix(Eqs, 0, 1, r, 1+a, &A);
  Matrix_subMatrix(Eqs, 0, 1+a, r, 1+a+b, &B);
  linearInter(A, B, &I, &Lb);
  Matrix_Free(A);
  Matrix_Free(B);
  Matrix_Free(I);
  if (dbgCompParm) {
    show_matrix(Lb);
  }
  
  /* 2- The linear part of the validity lattice is the left HNF of Lb */
  left_hermite(Lb, &H, &Q, &U);
  Matrix_Free(Lb);
  Matrix_Free(Q);
  Matrix_Free(U);

  /* 3- build the validity lattice */
  ensureMatrix((*vl), b+1, b+1);
  Matrix_copySubMatrix(H, 0, 0, b, b, (*vl), 0,0);
  Matrix_Free(H);
  for (i=0; i< b; i++) {
    value_assign((*vl)->p[i][b], sol->p[0][a+i]);
  }
  Matrix_Free(sol);
  Vector_Set((*vl)->p[b],0, b);
  value_set_si((*vl)->p[b][b], 1);
  
} /* validityLattice */


/**
 * Eliminate the columns corresponding to a list of eliminated parameters.
 * @param M the constraints matrix whose columns are to be removed
 * @param nbVars an offset to be added to the ranks of the variables to be
 * removed
 * @param elimParms the list of ranks of the variables to be removed
 * @param newM (output) the matrix without the removed columns
 */
void Constraints_removeElimCols(Matrix * M, unsigned int nbVars, 
			   unsigned int *elimParms, Matrix ** newM) {
  unsigned int i, j, k;
  if (elimParms[0]==0) {
    Matrix_clone(M, newM);
    return;
  }
  if ((*newM)==NULL) {
    (*newM) = Matrix_Alloc(M->NbRows, M->NbColumns - elimParms[0]);
  }
  else {
    assert ((*newM)->NbColumns==M->NbColumns - elimParms[0]);
  }
  for (i=0; i< M->NbRows; i++) {
    value_assign((*newM)->p[i][0], M->p[i][0]); /* kind of cstr */
    k=0;
    Vector_Copy(&(M->p[i][1]), &((*newM)->p[i][1]), nbVars);
    for (j=0; j< M->NbColumns-2-nbVars; j++) {
      if (j!=elimParms[k+1]) {
	value_assign((*newM)->p[i][j-k+nbVars+1], M->p[i][j+nbVars+1]);
      }
      else {
	k++;
      }
    }
    value_assign((*newM)->p[i][(*newM)->NbColumns-1], 
		 M->p[i][M->NbColumns-1]); /* cst part */
  }
} /* Constraints_removeElimCols */


/**
 * Eliminates all the equalities in a set of constraints and returns the set of
 * constraints defining a full-dimensional polyhedron, such that there is a
 * bijection between integer points of the original polyhedron and these of the
 * resulting (projected) polyhedron).
 * If VL is set to NULL, this funciton allocates it. Else, it assumes that
 * (*VL) points to a matrix of the right size.
 * <p> The following things are done: 
 * <ol>
 * <li> remove equalities involving only parameters, and remove as many
 *      parameters as there are such equalities. From that, the list of
 *      eliminated parameters <i>elimParms</i> is built.
 * <li> remove equalities that involve variables. This requires a compression
 *      of the parameters and of the other variables that are not eliminated.
 *      The affine compresson is represented by matrix VL (for <i>validity
 *      lattice</i>) and is such that (N I 1)^T = VL.(N' I' 1), where N', I'
 *      are integer (they are the parameters and variables after compression).
 *</ol>
 *</p>
 */
void Constraints_fullDimensionize(Matrix ** M, Matrix ** C, Matrix ** VL, 
				  Matrix ** Eqs, Matrix ** ParmEqs, 
				  unsigned int ** elimVars, 
				  unsigned int ** elimParms,
				  int maxRays) {
  unsigned int i, j;
  Matrix * A=NULL, *B=NULL;
  Matrix * Ineqs=NULL;
  unsigned int nbVars = (*M)->NbColumns - (*C)->NbColumns;
  unsigned int nbParms;
  int nbElimVars;
  Matrix * fullDim = NULL;

  /* variables for permutations */
  unsigned int * permutation;
  Matrix * permutedEqs=NULL, * permutedIneqs=NULL;
  
  /* 1- Eliminate the equalities involving only parameters. */
  (*ParmEqs) = Constraints_removeParmEqs(M, C, 0, elimParms);
  /* if the polyehdron is empty, return now. */
  if ((*M)->NbColumns==0) return;
  /* eliminate the columns corresponding to the eliminated parameters */
  if (elimParms[0]!=0) {
    Constraints_removeElimCols(*M, nbVars, (*elimParms), &A);
    Matrix_Free(*M);
    (*M) = A;
    Constraints_removeElimCols(*C, 0, (*elimParms), &B);
    Matrix_Free(*C);
    (*C) = B;
    if (dbgCompParm) {
      printf("After false parameter elimination: \n");
      show_matrix(*M);
      show_matrix(*C);
    }
  }
  nbParms = (*C)->NbColumns-2;

  /* 2- Eliminate the equalities involving variables */
  /*   a- extract the (remaining) equalities from the poyhedron */
  split_constraints((*M), Eqs, &Ineqs);
  nbElimVars = (*Eqs)->NbRows;
  /*    if the polyhedron is already full-dimensional, return */
  if ((*Eqs)->NbRows==0) {
    Matrix_identity(nbParms+1, VL);
    return;
  }
  /*   b- choose variables to be eliminated */
  permutation = find_a_permutation((*Eqs), nbParms);

  if (dbgCompParm) {
    printf("Permuting the vars/parms this way: [ ");
    for (i=0; i< (*Eqs)->NbColumns-2; i++) {
      printf("%d ", permutation[i]);
    }
    printf("]\n");
  }

  Constraints_permute((*Eqs), permutation, &permutedEqs);
  Equalities_validityLattice(permutedEqs, (*Eqs)->NbRows, VL);

  if (dbgCompParm) {
    printf("Validity lattice: ");
    show_matrix(*VL);
  }
  Constraints_compressLastVars(permutedEqs, (*VL));
  Constraints_permute(Ineqs, permutation, &permutedIneqs);
  if (dbgCompParmMore) {
    show_matrix(permutedIneqs);
    show_matrix(permutedEqs);
  }
  Matrix_Free(*Eqs);
  Matrix_Free(Ineqs);
  Constraints_compressLastVars(permutedIneqs, (*VL));
  if (dbgCompParm) {
    printf("After compression: ");
    show_matrix(permutedIneqs);
  }
  /*   c- eliminate the first variables */
  assert(Constraints_eliminateFirstVars(permutedEqs, permutedIneqs));
  if (dbgCompParmMore) {
    printf("After elimination of the variables: ");
    show_matrix(permutedIneqs);
  }

  /*   d- get rid of the first (zero) columns, 
       which are now useless, and put the parameters back at the end */
  fullDim = Matrix_Alloc(permutedIneqs->NbRows,
			 permutedIneqs->NbColumns-nbElimVars);
  for (i=0; i< permutedIneqs->NbRows; i++) {
    value_set_si(fullDim->p[i][0], 1);
    for (j=0; j< nbParms; j++) {
      value_assign(fullDim->p[i][j+fullDim->NbColumns-nbParms-1], 
		   permutedIneqs->p[i][j+nbElimVars+1]);
    }
    for (j=0; j< permutedIneqs->NbColumns-nbParms-2-nbElimVars; j++) {
      value_assign(fullDim->p[i][j+1], 
		   permutedIneqs->p[i][nbElimVars+nbParms+j+1]);
    }
    value_assign(fullDim->p[i][fullDim->NbColumns-1], 
		 permutedIneqs->p[i][permutedIneqs->NbColumns-1]);
  }
  Matrix_Free(permutedIneqs);

} /* Constraints_fullDimensionize */


/**
 * Given a matrix that defines a full-dimensional affine lattice, returns the 
 * affine sub-lattice spanned in the k first dimensions.
 * Useful for instance when you only look for the parameters' validity lattice.
 * @param lat the original full-dimensional lattice
 * @param subLat the sublattice
 */
void Lattice_extractSubLattice(Matrix * lat, unsigned int k, Matrix ** subLat) {
  Matrix * H, *Q, *U, *linLat = NULL;
  unsigned int i;
  dbgStart(Lattice_extractSubLattice);
  /* if the dimension is already good, just copy the initial lattice */
  if (k==lat->NbRows-1) {
    if (*subLat==NULL) {
      (*subLat) = Matrix_Copy(lat);
    }
    else {
      Matrix_copySubMatrix(lat, 0, 0, lat->NbRows, lat->NbColumns, (*subLat), 0, 0);
    }
    return;
  }
  assert(k<lat->NbRows-1);
  /* 1- Make the linear part of the lattice triangular to eliminate terms from 
     other dimensions */
  Matrix_subMatrix(lat, 0, 0, lat->NbRows, lat->NbColumns-1, &linLat);
  /* OPT: any integer column-vector elimination is ok indeed. */
  /* OPT: could test if the lattice is already in triangular form. */
  left_hermite(linLat, &H, &Q, &U);
  if (dbgCompParmMore) {
    show_matrix(H);
  }
  Matrix_Free(Q);
  Matrix_Free(U);
  Matrix_Free(linLat);
  /* if not allocated yet, allocate it */
  if (*subLat==NULL) {
    (*subLat) = Matrix_Alloc(k+1, k+1);
  }
  Matrix_copySubMatrix(H, 0, 0, k, k, (*subLat), 0, 0);
  Matrix_Free(H);
  Matrix_copySubMatrix(lat, 0, lat->NbColumns-1, k, 1, (*subLat), 0, k);
  for (i=0; i<k; i++) {
    value_set_si((*subLat)->p[k][i], 0);
  }
  value_set_si((*subLat)->p[k][k], 1);
  dbgEnd(Lattice_extractSubLattice);
} /* Lattice_extractSubLattice */


/** 
 * Computes the overall period of the variables I for (MI) mod |d|, where M is
 * a matrix and |d| a vector. Produce a diagonal matrix S = (s_k) where s_k is
 * the overall period of i_k 
 * @param M the set of affine functions of I (row-vectors)
 * @param d the column-vector representing the modulos
*/
Matrix * affine_periods(Matrix * M, Matrix * d) {
  Matrix * S;
  unsigned int i,j;
  Value tmp;
  Value * periods = (Value *)malloc(sizeof(Value) * M->NbColumns);
  value_init(tmp);
  for(i=0; i< M->NbColumns; i++) {
    value_init(periods[i]);
    value_set_si(periods[i], 1);
  }
  for (i=0; i<M->NbRows; i++) {
    for (j=0; j< M->NbColumns; j++) {
      Gcd(d->p[i][0], M->p[i][j], &tmp);
      value_division(tmp, d->p[i][0], tmp);
      Lcm3(periods[j], tmp, &(periods[j]));
     }
  }
  value_clear(tmp);

  /* 2- build S */
  S = Matrix_Alloc(M->NbColumns, M->NbColumns);
  for (i=0; i< M->NbColumns; i++) 
    for (j=0; j< M->NbColumns; j++)
      if (i==j) value_assign(S->p[i][j],periods[j]);
      else value_set_si(S->p[i][j], 0);

  /* 3- clean up */
  for(i=0; i< M->NbColumns; i++) value_clear(periods[i]);
  free(periods);
  return S;
} /* affine_periods */


/** 
 * Given an integer matrix B with m rows and integer m-vectors C and d,
 * computes the basis of the integer solutions to (BN+C) mod d = 0 (1).
 * This is an affine lattice (G): (N 1)^T= G(N' 1)^T, forall N' in Z^b.
 * If there is no solution, returns NULL.
 * @param B B, a (m x b) matrix
 * @param C C, a (m x 1) integer matrix
 * @param d d, a (1 x m) integer matrix
 * @param imb the affine (b+1)x(b+1) basis of solutions, in the homogeneous
 * form. Allocated if initially set to NULL, reused if not.
*/
void Equalities_intModBasis(Matrix * B, Matrix * C, Matrix * d, Matrix ** imb) {
  int b = B->NbColumns;
  /* FIXME: treat the case d=0 as a regular equality B_kN+C_k = 0: */
  /* OPT: could keep only equalities for which d>1 */
  int nbEqs = B->NbRows;
  unsigned int i;

  /* 1- buid the problem DI+BN+C = 0 */
  Matrix * eqs = Matrix_Alloc(nbEqs, nbEqs+b+1);
  for (i=0; i< nbEqs; i++) {
    value_assign(eqs->p[i][i], d->p[0][i]);
  }
  Matrix_copySubMatrix(B, 0, 0, nbEqs, b, eqs, 0, nbEqs);
  Matrix_copySubMatrix(C, 0, 0, nbEqs, 1, eqs, 0, nbEqs+b);

  /* 2- the solution is the validity lattice of the equalities */
  Equalities_validityLattice(eqs, nbEqs, imb);
  Matrix_Free(eqs);
} /* Equalities_intModBasis */


/** kept here for backwards compatiblity. Wrapper to Equalities_intModBasis() */
Matrix * int_mod_basis(Matrix * B, Matrix * C, Matrix * d) {
  Matrix * imb = NULL;
  Equalities_intModBasis(B, C, d, &imb);
  return imb;
} /* int_mod_basis */


/**
 * Given a parameterized constraints matrix with m equalities, computes the
 * compression matrix G such that there is an integer solution in the variables
 * space for each value of N', with N = G N' (N are the "nb_parms" parameters)
 * @param E a matrix of parametric equalities @param nb_parms the number of
 * parameters
 * <b>Note: </b>this function is mostly here for backwards
 * compatibility. Prefer the use of <tt>Equalities_validityLattice</tt>.
*/
Matrix * compress_parms(Matrix * E, int nbParms) {
  Matrix * vl=NULL;
  Equalities_validityLattice(E, E->NbColumns-2-nbParms, &vl);
  return vl;
}/* compress_parms */


/** Removes the equalities that involve only parameters, by eliminating some
 * parameters in the polyhedron's constraints and in the context.<p> 
 * <b>Updates M and Ctxt.</b>
 * @param M1 the polyhedron's constraints
 * @param Ctxt1 the constraints of the polyhedron's context
 * @param renderSpace tells if the returned equalities must be expressed in the
 * parameters space (renderSpace=0) or in the combined var/parms space
 * (renderSpace = 1)
 * @param elimParms the list of parameters that have been removed: an array
 * whose 1st element is the number of elements in the list.  (returned)
 * @return the system of equalities that involve only parameters.
 */
Matrix * Constraints_Remove_parm_eqs(Matrix ** M1, Matrix ** Ctxt1, 
				     int renderSpace, 
				     unsigned int ** elimParms) {
  int i, j, k, nbEqsParms =0;
  int nbEqsM, nbEqsCtxt, allZeros, nbTautoM = 0, nbTautoCtxt = 0;
  Matrix * M = (*M1);
  Matrix * Ctxt = (*Ctxt1);
  int nbVars = M->NbColumns-Ctxt->NbColumns;
  Matrix * Eqs;
  Matrix * EqsMTmp;
  
  /* 1- build the equality matrix(ces) */
  nbEqsM = 0;
  for (i=0; i< M->NbRows; i++) {
    k = First_Non_Zero(M->p[i], M->NbColumns);
    /* if it is a tautology, count it as such */
    if (k==-1) {
      nbTautoM++;
    }
    else {
      /* if it only involves parameters, count it */
      if (k>= nbVars+1) nbEqsM++;
    }
  }

  nbEqsCtxt = 0;
  for (i=0; i< Ctxt->NbRows; i++) {
    if (value_zero_p(Ctxt->p[i][0])) {
      if (First_Non_Zero(Ctxt->p[i], Ctxt->NbColumns)==-1) {
	nbTautoCtxt++;
      }
      else {
	nbEqsCtxt ++;
      }
    }
  }
  nbEqsParms = nbEqsM + nbEqsCtxt; 

  /* nothing to do in this case */
  if (nbEqsParms+nbTautoM+nbTautoCtxt==0) {
    (*elimParms) = (unsigned int*) malloc(sizeof(int));
    (*elimParms)[0] = 0;
    if (renderSpace==0) {
      return Matrix_Alloc(0,Ctxt->NbColumns);
    }
    else {
      return Matrix_Alloc(0,M->NbColumns);
    }
  }
  
  Eqs= Matrix_Alloc(nbEqsParms, Ctxt->NbColumns);
  EqsMTmp= Matrix_Alloc(nbEqsParms, M->NbColumns);
  
  /* copy equalities from the context */
  k = 0;
  for (i=0; i< Ctxt->NbRows; i++) {
    if (value_zero_p(Ctxt->p[i][0]) 
		     && First_Non_Zero(Ctxt->p[i], Ctxt->NbColumns)!=-1) {
      Vector_Copy(Ctxt->p[i], Eqs->p[k], Ctxt->NbColumns);
      Vector_Copy(Ctxt->p[i]+1, EqsMTmp->p[k]+nbVars+1, 
		  Ctxt->NbColumns-1);
      k++;
    }
  }
  for (i=0; i< M->NbRows; i++) {
    j=First_Non_Zero(M->p[i], M->NbColumns);
    /* copy equalities that involve only parameters from M */
    if (j>=nbVars+1) {
      Vector_Copy(M->p[i]+nbVars+1, Eqs->p[k]+1, Ctxt->NbColumns-1);
      Vector_Copy(M->p[i]+nbVars+1, EqsMTmp->p[k]+nbVars+1, 
		  Ctxt->NbColumns-1);
      /* mark these equalities for removal */
      value_set_si(M->p[i][0], 2);
      k++;
    }
    /* mark the all-zero equalities for removal */
    if (j==-1) {
      value_set_si(M->p[i][0], 2);
    }
  }

  /* 2- eliminate parameters until all equalities are used or until we find a
  contradiction (overconstrained system) */
  (*elimParms) = (unsigned int *) malloc((Eqs->NbRows+1) * sizeof(int));
  (*elimParms)[0] = 0;
  allZeros = 0;
  for (i=0; i< Eqs->NbRows; i++) {
    /* find a variable that can be eliminated */
    k = First_Non_Zero(Eqs->p[i], Eqs->NbColumns);
    if (k!=-1) { /* nothing special to do for tautologies */

      /* if there is a contradiction, return empty matrices */
      if (k==Eqs->NbColumns-1) {
	printf("Contradiction in %dth row of Eqs: ",k);
	show_matrix(Eqs);
	Matrix_Free(Eqs);
	Matrix_Free(EqsMTmp);
	(*M1) = Matrix_Alloc(0, M->NbColumns);
	Matrix_Free(M);
	(*Ctxt1) = Matrix_Alloc(0,Ctxt->NbColumns);
	Matrix_Free(Ctxt);
	free(*elimParms);
	(*elimParms) = (unsigned int *) malloc(sizeof(int));
	(*elimParms)[0] = 0;
	if (renderSpace==1) {
	  return Matrix_Alloc(0,(*M1)->NbColumns);
	}
	else {
	  return Matrix_Alloc(0,(*Ctxt1)->NbColumns);
	}
      }	
      /* if we have something we can eliminate, do it in 3 places:
	 Eqs, Ctxt, and M */
      else {
	k--; /* k is the rank of the variable, now */
	(*elimParms)[0]++;
	(*elimParms)[(*elimParms[0])]=k;
	for (j=0; j< Eqs->NbRows; j++) {
	  if (i!=j) {
	    eliminate_var_with_constr(Eqs, i, Eqs, j, k);
	    eliminate_var_with_constr(EqsMTmp, i, EqsMTmp, j, k+nbVars);
	  }
	}
	for (j=0; j< Ctxt->NbRows; j++) {
	  if (value_notzero_p(Ctxt->p[i][0])) {
	    eliminate_var_with_constr(Eqs, i, Ctxt, j, k);
	  }
	}
	for (j=0; j< M->NbRows; j++) {
	  if (value_cmp_si(M->p[i][0], 2)) {
	    eliminate_var_with_constr(EqsMTmp, i, M, j, k+nbVars);
	  }
	}
      }
    }
    /* if (k==-1): count the tautologies in Eqs to remove them later */
    else {
      allZeros++;
    }
  }
  
  /* elimParms may have been overallocated. Now we know how many parms have
     been eliminated so we can reallocate the right amount of memory. */
  if (!realloc((*elimParms), ((*elimParms)[0]+1)*sizeof(int))) {
    fprintf(stderr, "Constraints_Remove_parm_eqs > cannot realloc()");
  }

  Matrix_Free(EqsMTmp);

  /* 3- remove the "bad" equalities from the input matrices
     and copy the equalities involving only parameters */
  EqsMTmp = Matrix_Alloc(M->NbRows-nbEqsM-nbTautoM, M->NbColumns);
  k=0;
  for (i=0; i< M->NbRows; i++) {
    if (value_cmp_si(M->p[i][0], 2)) {
      Vector_Copy(M->p[i], EqsMTmp->p[k], M->NbColumns);
      k++;
    }
  }
  Matrix_Free(M);
  (*M1) = EqsMTmp;
  
  EqsMTmp = Matrix_Alloc(Ctxt->NbRows-nbEqsCtxt-nbTautoCtxt, Ctxt->NbColumns);
  k=0;
  for (i=0; i< Ctxt->NbRows; i++) {
    if (value_notzero_p(Ctxt->p[i][0])) {
      Vector_Copy(Ctxt->p[i], EqsMTmp->p[k], Ctxt->NbColumns);
      k++;
    }
  }
  Matrix_Free(Ctxt);
  (*Ctxt1) = EqsMTmp;
  
  if (renderSpace==0) {/* renderSpace=0: equalities in the parameter space */
    EqsMTmp = Matrix_Alloc(Eqs->NbRows-allZeros, Eqs->NbColumns);
    k=0;
    for (i=0; i<Eqs->NbRows; i++) {
      if (First_Non_Zero(Eqs->p[i], Eqs->NbColumns)!=-1) {
	Vector_Copy(Eqs->p[i], EqsMTmp->p[k], Eqs->NbColumns);
	k++;
      }
    }
  }
  else {/* renderSpace=1: equalities rendered in the combined space */
    EqsMTmp = Matrix_Alloc(Eqs->NbRows-allZeros, (*M1)->NbColumns);
    k=0;
    for (i=0; i<Eqs->NbRows; i++) {
      if (First_Non_Zero(Eqs->p[i], Eqs->NbColumns)!=-1) {
	Vector_Copy(Eqs->p[i], &(EqsMTmp->p[k][nbVars]), Eqs->NbColumns);
	k++;
      }
    }
  }
  Matrix_Free(Eqs);
  Eqs = EqsMTmp;

  return Eqs;
} /* Constraints_Remove_parm_eqs */


/** Removes equalities involving only parameters, but starting from a
 * Polyhedron and its context.
 * @param P the polyhedron
 * @param C P's context
 * @param renderSpace: 0 for the parameter space, =1 for the combined space.
 * @maxRays Polylib's usual <i>workspace</i>.
 */
Polyhedron * Polyhedron_Remove_parm_eqs(Polyhedron ** P, Polyhedron ** C, 
					int renderSpace, 
					unsigned int ** elimParms, 
					int maxRays) {
  Matrix * Eqs;
  Polyhedron * Peqs;
  Matrix * M = Polyhedron2Constraints((*P));
  Matrix * Ct = Polyhedron2Constraints((*C));

  /* if the Minkowski representation is not computed yet, do not compute it in
     Constraints2Polyhedron */
  if (F_ISSET((*P), POL_VALID | POL_INEQUALITIES) && 
      (F_ISSET((*C), POL_VALID | POL_INEQUALITIES))) {
    FL_INIT(maxRays, POL_NO_DUAL);
  }
    
  Eqs = Constraints_Remove_parm_eqs(&M, &Ct, renderSpace, elimParms);
  Peqs = Constraints2Polyhedron(Eqs, maxRays);
  Matrix_Free(Eqs);

  /* particular case: no equality involving only parms is found */
  if (Eqs->NbRows==0) {
    Matrix_Free(M);
    Matrix_Free(Ct);
    return Peqs;
  }
  Polyhedron_Free(*P);
  Polyhedron_Free(*C);
  (*P) = Constraints2Polyhedron(M, maxRays);
  (*C) = Constraints2Polyhedron(Ct, maxRays);
  Matrix_Free(M);
  Matrix_Free(Ct);
  return Peqs;
} /* Polyhedron_Remove_parm_eqs */


/**
 * Given a matrix with m parameterized equations, compress the nb_parms
 * parameters and n-m variables so that m variables are integer, and transform
 * the variable space into a n-m space by eliminating the m variables (using
 * the equalities) the variables to be eliminated are chosen automatically by
 * the function.
 * <b>Deprecated.</b> Try to use Constraints_fullDimensionize instead.
 * @param M the constraints 
 * @param the number of parameters
 * @param validityLattice the the integer lattice underlying the integer
 * solutions.
*/
Matrix * full_dimensionize(Matrix const * M, int nbParms, 
			   Matrix ** validityLattice) {
  Matrix * Eqs, * Ineqs;
  Matrix * permutedEqs, * permutedIneqs;
  Matrix * Full_Dim;
  Matrix * WVL; /* The Whole Validity Lattice (vars+parms) */
  unsigned int i,j;
  int nbElimVars;
  unsigned int * permutation, * permutationInv;
  /* 0- Split the equalities and inequalities from each other */
  split_constraints(M, &Eqs, &Ineqs);

  /* 1- if the polyhedron is already full-dimensional, return it */
  if (Eqs->NbRows==0) {
    Matrix_Free(Eqs);
    (*validityLattice) = Identity_Matrix(nbParms+1);
    return Ineqs;
  }
  nbElimVars = Eqs->NbRows;

  /* 2- put the vars to be eliminated at the first positions, 
     and compress the other vars/parms
     -> [ variables to eliminate / parameters / variables to keep ] */
  permutation = find_a_permutation(Eqs, nbParms);
  if (dbgCompParm) {
    printf("Permuting the vars/parms this way: [ ");
    for (i=0; i< Eqs->NbColumns; i++) {
      printf("%d ", permutation[i]);
    }
    printf("]\n");
  }
  permutedEqs = mpolyhedron_permute(Eqs, permutation);
  WVL = compress_parms(permutedEqs, Eqs->NbColumns-2-Eqs->NbRows);
  if (dbgCompParm) {
    printf("Whole validity lattice: ");
    show_matrix(WVL);
  }
  mpolyhedron_compress_last_vars(permutedEqs, WVL);
  permutedIneqs = mpolyhedron_permute(Ineqs, permutation);
  if (dbgCompParm) {
    show_matrix(permutedEqs);
  }
  Matrix_Free(Eqs);
  Matrix_Free(Ineqs);
  mpolyhedron_compress_last_vars(permutedIneqs, WVL);
  if (dbgCompParm) {
    printf("After compression: ");
    show_matrix(permutedIneqs);
  }
  /* 3- eliminate the first variables */
  if (!mpolyhedron_eliminate_first_variables(permutedEqs, permutedIneqs)) {
    fprintf(stderr,"full-dimensionize > variable elimination failed. \n"); 
    return NULL;
  }
  if (dbgCompParm) {
    printf("After elimination of the variables: ");
    show_matrix(permutedIneqs);
  }

  /* 4- get rid of the first (zero) columns, 
     which are now useless, and put the parameters back at the end */
  Full_Dim = Matrix_Alloc(permutedIneqs->NbRows,
			  permutedIneqs->NbColumns-nbElimVars);
  for (i=0; i< permutedIneqs->NbRows; i++) {
    value_set_si(Full_Dim->p[i][0], 1);
    for (j=0; j< nbParms; j++) 
      value_assign(Full_Dim->p[i][j+Full_Dim->NbColumns-nbParms-1], 
		   permutedIneqs->p[i][j+nbElimVars+1]);
    for (j=0; j< permutedIneqs->NbColumns-nbParms-2-nbElimVars; j++) 
      value_assign(Full_Dim->p[i][j+1], 
		   permutedIneqs->p[i][nbElimVars+nbParms+j+1]);
    value_assign(Full_Dim->p[i][Full_Dim->NbColumns-1], 
		 permutedIneqs->p[i][permutedIneqs->NbColumns-1]);
  }
  Matrix_Free(permutedIneqs);
  
  /* 5- Keep only the the validity lattice restricted to the parameters */
  *validityLattice = Matrix_Alloc(nbParms+1, nbParms+1);
  for (i=0; i< nbParms; i++) {
    for (j=0; j< nbParms; j++)
      value_assign((*validityLattice)->p[i][j], 
		   WVL->p[i][j]);
    value_assign((*validityLattice)->p[i][nbParms], 
		 WVL->p[i][WVL->NbColumns-1]);
  }
  for (j=0; j< nbParms; j++) 
    value_set_si((*validityLattice)->p[nbParms][j], 0);
  value_assign((*validityLattice)->p[nbParms][nbParms], 
	       WVL->p[WVL->NbColumns-1][WVL->NbColumns-1]);

  /* 6- Clean up */
  Matrix_Free(WVL);
  return Full_Dim;
} /* full_dimensionize */

#undef dbgCompParm
#undef dbgCompParmMore
