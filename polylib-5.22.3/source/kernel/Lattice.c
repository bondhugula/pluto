#include <stdlib.h>
#include <polylib/polylib.h>


typedef struct {
  int count;
  int *fac;
} factor;

static factor allfactors (int num);

/* 
 * Print the contents of a list of Lattices 'Head'  
 */
void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head) {
  
  LatticeUnion *temp;

  for(temp = Head; temp != NULL; temp = temp->next)
    Matrix_Print(fp,format,(Matrix *)temp->M);
  return;
} /* PrintLatticeUnion */

/* 
 * Free the memory allocated to a list of lattices 'Head' 
 */
void LatticeUnion_Free(LatticeUnion *Head) {

  LatticeUnion  *temp;

  while (Head != NULL) {
    temp = Head;
    Head = temp->next;
    Matrix_Free(temp->M);
    free(temp);
  }
  return;
} /* LatticeUnion_Free */

/* 
 * Allocate a heads for a list of Lattices
 */
LatticeUnion *LatticeUnion_Alloc(void) {

  LatticeUnion  *temp;

  temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
  temp->M=NULL;
  temp->next=NULL;
  return temp;
} /* LatticeUnion_Alloc */

/*
 * Given two Lattices 'A' and 'B', return True if they have the same affine 
 * part (the last column) otherwise return 'False'. 
 */
Bool sameAffinepart (Lattice *A, Lattice *B) {
  
  int i;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered SAMEAFFINEPART \n"); 
  fclose(fp);
#endif
  
  for (i = 0; i < A->NbRows; i ++)
    if (value_ne(A->p[i][A->NbColumns-1],B->p[i][B->NbColumns-1]))
      return False;
  return True;
} /* sameAffinepart */ 

/*
 * Return an empty lattice of dimension 'dimension-1'. An empty lattice is 
 * represented as [[0 0 ... 0] .... [0 ... 0][0 0.....0 1]]. 
 */ 
Lattice *EmptyLattice(int dimension) {

  Lattice *result;
  int i,j;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered NULLATTICE \n"); 
  fclose (fp);
#endif
  
  result = (Lattice *) Matrix_Alloc(dimension, dimension);
  for (i = 0; i < dimension; i ++)
    for (j = 0; j < dimension; j ++)
      value_set_si(result->p[i][j],0);
  value_set_si(result->p[i-1][i-1],1);
  return result;
} /* EmptyLattice */ 

/*
 * Return True if Lattice 'A' is empty, otherwise return False. 
 */
Bool isEmptyLattice (Lattice *A) {
  
  int i,j;  
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ISNULLATTICE \n"); 
  fclose(fp);
#endif
  
  for (i = 0; i < A->NbRows-1; i ++)
    for (j = 0; j < A->NbColumns-1; j ++)
      if(value_notzero_p(A->p[i][j])) {
	return False;
      }
  if (value_one_p(A->p[i][A->NbColumns-1])) {
    return True ;
  }
  return False ;
} /* isEmptyLaattice */ 

/*
 * Given a Lattice 'A', check whether it is linear or not, i.e. whether the 
 * affine part is NULL or not. If affine part is empty, it returns True other-
 * wise it returns False.
 */
Bool isLinear(Lattice  *A) {
  
  int i;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered ISLINEAR \n"); 
  fclose (fp);
#endif

  for (i = 0; i < A->NbRows-1; i ++)
    if (value_notzero_p(A->p[i][A->NbColumns-1])) {
      return False;
    }
  return True;
} /* isLinear */ 
          
/* 
 * Return the affine Hermite normal form of the affine lattice 'A'. The unique 
 * affine Hermite form if a lattice is stored in 'H' and the unimodular matrix 
 * corresponding to 'A = H*U' is stored in the matrix 'U'. 
 * Algorithm : 
 *            1) Check if the Lattice is Linear or not.
 *            2) If it is not Linear, then Homogenise the Lattice.
 *            3) Call Hermite.
 *            4) If the Lattice was Homogenised, the HNF H must be 
 *               Dehomogenised and also corresponding changes must
 *               be made to the Unimodular Matrix U.
 *            5) Return.
 */ 
void AffineHermite (Lattice *A, Lattice **H, Matrix **U) {
 
  Lattice *temp;
  Bool flag = True;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered AFFINEHERMITE \n"); 
  fclose (fp);
#endif

  if (isLinear(A) == False)
    temp = Homogenise(A,True); 
  else {
    flag = False ;
    temp = (Lattice *)Matrix_Copy(A);
  }  
  Hermite((Matrix *)temp,(Matrix **) H, U);
  if (flag == True) {
    Matrix_Free ((Matrix *) temp);
    temp = Homogenise(H[0],False);
    Matrix_Free((Matrix *) H[0]);
    H[0] = (Lattice *)Matrix_Copy(temp);
    Matrix_Free((Matrix *) temp);
    temp = Homogenise(U[0],False);
    Matrix_Free ((Matrix *) U[0]);
    U[0] = (Matrix *)Matrix_Copy(temp);
  }  
  Matrix_Free((Matrix *) temp);
  return;
} /* AffineHermite */

/*
 * Given a Polylib matrix 'A' that rerepresents an affine function, return the
 * affine Smith normal form 'Delta' of 'A' and unimodular matrices 'U' and 'V'
 * such that 'A = U*Delta*V'. 
 * Algorithm:
 *           (1) Homogenise the Lattice.
 *           (2) Call Smith
 *           (3) The Smith Normal Form Delta must be Dehomogenised and also 
 *               corresponding changes must be made to the Unimodular Matrices 
 *               U and V.
 *           4) Bring Delta into AffineSmith Form.
 */
void AffineSmith(Lattice *A, Lattice **U, Lattice **V, Lattice **Diag) {
 
  Lattice *temp;
  Lattice *Uinv;
  int i,j;
  Value sum, quo, rem;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered AFFINESMITH \n"); 
  fclose(fp);
#endif
  
  value_init(sum); 
  value_init(quo); value_init(rem);
  temp = Homogenise(A,True);  
  Smith((Matrix *)temp, (Matrix **)U, (Matrix **)V, (Matrix **)Diag);
  Matrix_Free((Matrix *)temp);
  
  temp = Homogenise (*U, False);
  Matrix_Free ((Matrix *) *U);
  *U = (Lattice *)Matrix_Copy ((Matrix *)temp);
  Matrix_Free ((Matrix *)temp);
  
  temp = Homogenise (*V, False);
  Matrix_Free ((Matrix *)*V);
  *V = (Lattice *) Matrix_Copy ((Matrix *)temp);
  Matrix_Free ((Matrix *)temp);
  
  temp = Homogenise (*Diag, False);
  Matrix_Free ((Matrix *)*Diag);
  *Diag = (Lattice *)Matrix_Copy ((Matrix *)temp);
  Matrix_Free ((Matrix *)temp);
  
  temp = (Lattice *) Matrix_Copy ((Matrix *) *U);
  Uinv = (Lattice *) Matrix_Alloc (U[0]->NbRows, U[0]->NbColumns);
  Matrix_Inverse( (Matrix *) temp, (Matrix *)  Uinv);
  Matrix_Free ((Matrix *) temp);
  
  for (i = 0; i < U[0]->NbRows-1; i ++) {
    value_set_si(sum,0);
    for(j = 0; j < U[0]->NbColumns-1; j ++) {
      value_addmul(sum, Uinv->p[i][j], U[0]->p[j][U[0]->NbColumns-1]);
    }
    value_assign(Diag[0]->p[i][j],sum);
  }
  Matrix_Free((Matrix *) Uinv);  
  for(i = 0; i < U[0]->NbRows-1; i ++) 
    value_set_si(U[0]->p[i][U[0]->NbColumns-1],0);
  for(i = 0; i < Diag[0]->NbRows-1; i ++) {
    value_division(quo,Diag[0]->p[i][Diag[0]->NbColumns-1],Diag[0]->p[i][i]);
    value_modulus(rem,Diag[0]->p[i][Diag[0]->NbColumns-1],Diag[0]->p[i][i]);
    
    fprintf(stdout," pourcent "); 
    value_print(stdout,VALUE_FMT,rem);
    fprintf(stdout," quotient ");
    value_print(stdout,VALUE_FMT,quo);
    fprintf(stdout," \n");
    
    /* Apparently the % operator is strange when sign are different */
    if(value_neg_p(rem)) {
      value_addto(rem,rem,Diag[0]->p[i][i]);
      value_decrement(quo,quo);
    };
    fprintf(stdout,"apres  pourcent "); 
    value_print(stdout,VALUE_FMT,rem);
    fprintf(stdout," quotient ");
    value_print(stdout,VALUE_FMT,quo);
    fprintf(stdout," \n");
    value_assign( Diag[0]->p[i][Diag[0]->NbColumns-1],rem);
    value_assign(V[0]->p[i][V[0]->NbColumns-1],quo);
  }  
  value_clear(sum); 
  value_clear(quo); value_clear(rem);
  return;
} /* AffineSmith */

/*
 * Given a lattice 'A' and a boolean variable 'Forward', homogenise the lattice
 * if 'Forward' is True, otherwise if 'Forward' is False, dehomogenise the 
 * lattice 'A'. 
 * Algorithm: 
 *            (1) If Forward == True
 *                Put the last row first. 
 *                Put the last columns first. 
 *            (2) Else 
 *                Put the first row last. 
 *                Put the first column last.
 *            (3) Return the result. 
 */
Lattice *Homogenise(Lattice *A, Bool Forward) {
  
  Lattice *result;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered HOMOGENISE \n"); 
  fclose(fp);
#endif
    
  result = (Lattice *)Matrix_Copy(A);  
  if (Forward == True ) { 
    PutColumnFirst((Matrix *)result, A->NbColumns-1);
    PutRowFirst((Matrix *)result, result->NbRows-1);
  }
  else  { 
    PutColumnLast((Matrix *)result,0);
    PutRowLast((Matrix *)result,0);
  }   
  return result;
} /* Homogenise */ 

/*
 * Given two lattices 'A' and 'B', verify if lattice 'A' is included in 'B' or
 * not. If 'A' is included in 'B' the 'A' intersection 'B', will be 'A'. So, 
 * compute 'A' intersection 'B' and check if it is the same as 'A'. 
 */
Bool LatticeIncludes(Lattice *A, Lattice *B) {
  
  Lattice *temp, *UA, *HA;
  Bool flag = False;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered LATTICE INCLUDES \n"); 
  fclose(fp);
#endif
  
  AffineHermite(A,&HA,&UA);  
  temp = LatticeIntersection(B,HA); 
  if (sameLattice(temp, HA) == True)
    flag = True;
  
  Matrix_Free((Matrix *)temp);
  Matrix_Free((Matrix *)UA);
  Matrix_Free((Matrix *)HA);
  return flag; 
} /* LatticeIncludes */

/*
 * Given two lattices 'A' and 'B', verify if 'A' and 'B' are the same lattice.
 * Algorithm: 
 *           The Affine Hermite form of two full dimensional matrices are 
 * unique. So, take the Affine Hermite form of both 'A' and 'B' and compare the
 * matrices. If they are equal, the function returns True, else it returns 
 * False. 
 */
Bool sameLattice(Lattice *A, Lattice *B) {
  
  Lattice *HA, *HB, *UA, *UB;
  int i,j;
  Bool result = True;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered SAME LATTICE \n"); 
  fclose(fp);
#endif
  
  AffineHermite(A, &HA, &UA);
  AffineHermite(B, &HB, &UB);
    
  for (i = 0 ; i < A->NbRows; i ++)
    for (j =  0; j < A->NbColumns; j ++)
      if (value_ne(HA->p[i][j],HB->p[i][j])) {	
	result = False; 
	break; 
      }  
  
  Matrix_Free ((Matrix *) HA); 
  Matrix_Free ((Matrix *) HB); 
  Matrix_Free ((Matrix *) UA); 
  Matrix_Free ((Matrix *) UB); 
  
  return result;
} /* sameLattice */ 

/*
 * Given a matrix 'A' and an integer 'dimension', do the following: 
 * If dimension < A->dimension), output a (dimension * dimension) submatrix of 
 * A. Otherwise the output matrix is [A 0][0 ID]. The order if the identity 
 * matrix is (dimension - A->dimension). The input matrix is not necessarily 
 * a Polylib matrix but the output is a polylib matrix. 
 */
Lattice *ChangeLatticeDimension(Lattice *A, int dimension) {
  
  int i, j;
  Lattice *Result ;
  
  Result = Matrix_Alloc(dimension, dimension);
  if(dimension <= A->NbRows) {
    for (i = 0; i < dimension; i ++)
      for (j = 0; j < dimension; j ++)
	value_assign(Result->p[i][j],A->p[i][j]);
    return Result;
  }  
  for (i = 0; i < A->NbRows; i ++)
    for (j = 0; j < A->NbRows; j ++)
      value_assign(Result->p[i][j],A->p[i][j]);
  
  for (i = A->NbRows; i < dimension; i ++)
    for (j = 0; j < dimension; j ++) {
	value_set_si(Result->p[i][j],0);
	value_set_si(Result->p[j][i],0);
    }  
  for (i = A->NbRows; i < dimension; i ++)
     value_set_si(Result->p[i][i],1);   
  return Result;
} /* ChangeLatticeDimension */

/* 
 * Given an affine lattice 'A', return a matrix of the linear part of the 
 * lattice.  
 */
Lattice *ExtractLinearPart(Lattice *A) {

  Lattice *Result;
  int i, j; 
  Result = (Lattice *) Matrix_Alloc(A->NbRows-1, A->NbColumns-1);
  for (i = 0; i < A->NbRows-1; i ++)
    for (j = 0; j < A->NbColumns-1; j ++)
      value_assign(Result->p[i][j],A->p[i][j]);  
  return Result;
} /* ExtractLinearPart */

static Matrix *MakeDioEqforInter(Matrix *A, Matrix *B);

/*
 * Given two lattices 'A' and 'B', return the intersection of the two lattcies.
 * The dimension of 'A' and 'B' should be the same. 
 * Algorithm:
 *           (1) Verify if the lattcies 'A' and 'B' have the same affine part. 
 *               If they have same affine part, then only their Linear parts 
 *               need to be intersected. If they don't have the same affine
 *               part then the affine part has to be taken into consideration. 
 *               For this, homogenise the lattices to get their Hermite Forms
 *               and then find their intersection.
 *
 *           (2) Step(2) involves, solving the Diophantine Equations in order 
 *               to extract the intersection of the Lattices. The Diophantine
 *               equations are formed taking into consideration whether the 
 *               affine part has to be included or not. 
 *
 *           (3) Solve the Diophantine equations. 
 *
 *           (4) Extract the necessary information from the result. 
 * 
 *           (5) If the lattices have different affine parts and they were 
 *               homogenised, the result is dehomogenised. 
 */ 
Lattice *LatticeIntersection(Lattice *X, Lattice *Y) {
  
  int i, j, exist;
  Lattice *result = NULL, *U = NULL ;
  Lattice *A = NULL, *B = NULL, *H = NULL;
  Matrix *fordio;
  Vector *X1 = NULL;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered LATTICEINTERSECTION \n"); 
  fclose(fp);
#endif

  if (X->NbRows != X->NbColumns) {
    fprintf(stderr, "\nIn LatticeIntersection : The Input Matrix X is a not a well defined Lattice\n");
    return EmptyLattice(X->NbRows);
  }
  
  if (Y->NbRows != Y->NbColumns) {
    fprintf (stderr, "\nIn LatticeIntersection : The Input Matrix Y is a not a well defined Lattice\n");
    return EmptyLattice(X->NbRows);
  }
  
  if (Y->NbRows != X->NbRows) {
    fprintf (stderr, "\nIn LatticeIntersection : The Input Lattices X and Y are of incompatible dimensions\n");
    return EmptyLattice(X->NbRows);
  }
  
  if (isinHnf(X))
    A = (Lattice *) Matrix_Copy(X);
  else {
    AffineHermite(X, &H, &U);
    A = (Lattice *)Matrix_Copy (H);
    Matrix_Free((Matrix *) H);
    Matrix_Free((Matrix *) U);
  }
  
  if (isinHnf(Y))
    B = (Lattice *)Matrix_Copy(Y);
  else {
    AffineHermite(Y, &H, &U);
    B = (Lattice *)Matrix_Copy (H);
    Matrix_Free((Matrix *) H);
    Matrix_Free((Matrix *) U);
  }
  
  if ((isEmptyLattice(A)) || (isEmptyLattice (B))) {
    result = EmptyLattice(X->NbRows);
    Matrix_Free ((Matrix *) A); 
    Matrix_Free ((Matrix *) B);
    return result;   
  }
  fordio = MakeDioEqforInter (A, B);
  Matrix_Free (A);
  Matrix_Free (B);
  exist = SolveDiophantine(fordio,(Matrix **) &U, &X1);
  if (exist < 0) { /* Intersection is NULL */    
    result = (EmptyLattice(X->NbRows)); 
    return result;
  }
  
  result = (Lattice *)Matrix_Alloc(X->NbRows, X->NbColumns);
  for (i = 0; i < result->NbRows-1; i ++)
    for (j = 0; j < result->NbColumns-1; j ++)
      value_assign(result->p[i][j],U->p[i][j]);
  
  for (i = 0; i < result->NbRows-1; i ++)
    value_assign(result->p[i][result->NbColumns-1],X1->p[i]); 
  for (i = 0; i < result->NbColumns-1; i ++)
    value_set_si(result->p[result->NbRows-1][i],0);
  value_set_si(result->p[result->NbRows-1][result->NbColumns-1],1);
  
  Matrix_Free((Matrix *) U);
  Vector_Free(X1);
  Matrix_Free(fordio);
  
  AffineHermite(result,&H,&U);  
  Matrix_Free((Matrix *)result);
  result = (Lattice *)Matrix_Copy(H); 
  
  Matrix_Free((Matrix *) H);
  Matrix_Free((Matrix *) U);
  
  /* Check whether the Lattice is NULL or not */
  
  if (isEmptyLattice (result)) {
    Matrix_Free ((Matrix *)result);
    return (EmptyLattice (X->NbRows));
  }  
  return result;
} /* LatticeIntersection */

static Matrix * MakeDioEqforInter (Lattice *A, Lattice *B) {
  
  Matrix *Dio ;
  int i,j;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered MAKEDIOEQFORINTER \n"); 
  fclose(fp);
#endif
  
 Dio = Matrix_Alloc(2*(A->NbRows-1) + 1, 3 * (A->NbColumns-1)+1);
 
 for (i = 0; i < Dio->NbRows; i ++)
   for (j = 0; j < Dio->NbColumns; j ++)
     value_set_si(Dio->p[i][j],0);
 
 for (i = 0; i < A->NbRows-1; i++) {
   value_set_si(Dio->p[i][i],1);
   value_set_si(Dio->p[i+A->NbRows-1][i],1);  
 } 
 for (i = 0; i < A->NbRows-1 ; i ++)
   for (j = 0; j < A->NbRows-1; j ++) {
     value_oppose(Dio->p[i][j+A->NbRows-1],A->p[i][j]);
     value_oppose(Dio->p[i+(A->NbRows-1)][j+2*(A->NbRows-1)],B->p[i][j]);
   }
 
 /* Adding the affine part */
 
 for (i = 0; i < A->NbColumns-1; i++) {
   value_oppose(Dio->p[i][Dio->NbColumns-1],A->p[i][A->NbColumns-1]);
   value_oppose(Dio->p[i+A->NbRows-1][Dio->NbColumns-1],B->p[i][A->NbColumns-1]) ;
 } 
 value_set_si(Dio->p[Dio->NbRows-1][Dio->NbColumns-1],1); 
 return Dio;
} /* MakeDioEqforInter */

static void AddLattice(LatticeUnion *,Matrix *,  Matrix *, int , int);
LatticeUnion *SplitLattice(Matrix *, Matrix *, Matrix *);



/*
 * The function is transforming a lattice X in a union of lattices based on a starting lattice Y.
 * Note1: If the intersection of X and Y lattices is empty the result is identic with the first argument (X) because no operation can be made.
 *Note2: The function is availabe only for simple Lattices and not for a union of Lattices.

 *       Step 1:  Find Intersection = LatticeIntersection (A, B).
 *       Step 2:  Extract the Linear Parts of the Lattices A and Intersection.
 *                (while dealing with Basis we only deal with the Linear Parts)
 *       Step 3:  Let M1 = Basis of A and M2 = Basis of B.
 *                Let B1 and B2 be the Basis of A and B respectively, 
 *                corresponding to the above Theorem.
 *                Then we Have B1 = M1 * U1 {a unimodular Matrix }
 *                and B2 = M2 * U2. M1 and M2 we know, they are the linear 
 *                parts we obtained in Step 2. Our Task is now to find U1 and
 *                U2. 
 *                We know that B1  * Delta = B2.
 *                i.e. M1 * U1 * Delta = M2 * U2
 *                or U1*Delta*U2Inverse = M1Inverse * M2.
 *                and Delta is the Diagonal Matrix which satisifies the 
 *                above properties (in the Theorem).
 *                So Delta is nothing but the Smith Normal Form of 
 *                M1Inverse * M2.
 *                So, first we have to find M1Inverse.
 *             
 *                This Step, involves finding the Inverse of the Matrix M1.
 *                We find the Inverse using the Polylib function 
 *                Matrix_Inverse. There is a catch here, the result of this
 *                function is an integral matrix, not necessarily the exact
 *                Inverse (since M1 need not be Unimodular), but a multiple
 *                of the actual inverse. The number by which we have to divide
 *                the matrix, is not obtained here as the input matrix is not
 *                a Polylib matrix { We input only the Linear part }. Later I
 *                give a way for finding that number.
 *
 *                M1Inverse = Matrix_Inverse ( M1 );
 *      
 *      Step 4 :  MtProduct = Matrix_Product (M1Inverse, M2);
 *      Step 5 :  SmithNormalFrom (MtProduct, Delta, U, V);
 *                U1 = U and U2Inverse = V.
 *      Step 6 :  Find U2 = Matrix_Inverse  (U2inverse). Here there is no prob
 *                as U1 and its inverse are unimodular.
 *      
 *      Step 7 :  Compute B1 = M1 * U1;
 *      Step 8 :  Compute B2 = M2 * U2;
 *      Step 9 :  Earlier when we computed M1Inverse, we knew that it was not
 *                the exact inverse but a multiple of it. Now we find the 
 *                number, such that ( M1Inverse / number ) would give us the 
 *                exact inverse of M1.
 *                We know that B1 * Delta = B2.
 *                Let k = B2[0][0] / B1[0][0].
 *                Let number = Delta[0][0]/k;
 *                This 'number' is the number we want.
 *                We Divide the matrix Delta by this number, to get the actual
 *                Delta such that B1 * Delta = B2.
 *     Step 10 :  Call Split Lattice (B1, B2, Delta ).
 *                This function returns the Union of Lattices in such a way 
 *                that B2 is at the Head of this List.
 *
 *If the intersection between X and Y is empty then the result is NULL.
 */


LatticeUnion *Lattice2LatticeUnion(Lattice *X,Lattice *Y)
{
  Lattice *B1 = NULL, *B2 = NULL, *newB1 = NULL, *newB2 = NULL, *Intersection=NULL;
  Matrix *U = NULL,*M1 = NULL, *M2 = NULL, *M1Inverse = NULL,*MtProduct = NULL;
  Matrix *Vinv, *V , *temp, *DiagMatrix ;

  LatticeUnion *Head = NULL, *tempHead = NULL;
  int i;
  Value k;
  

  Intersection = LatticeIntersection(X,Y);
  if (isEmptyLattice(Intersection) == True) {
    fprintf(stderr,"\nIn Lattice2LatticeUnion : The Input Lattices X and Y does not have any common part\n");
    return NULL;
  }  

  value_init(k);
  M1 = (Matrix *)ExtractLinearPart(X);
  M2 = (Matrix *)ExtractLinearPart(Intersection);

  M1Inverse = Matrix_Alloc(M1->NbRows,M1->NbColumns);
  temp = Matrix_Copy(M1);
  Matrix_Inverse(temp,M1Inverse);
  Matrix_Free(temp);

  MtProduct = Matrix_Alloc(M1->NbRows, M1->NbColumns);
  Matrix_Product(M1Inverse,M2,MtProduct) ;  
  Smith(MtProduct, &U, &Vinv, &DiagMatrix);  
  V = Matrix_Alloc(Vinv->NbRows,Vinv->NbColumns);
  Matrix_Inverse(Vinv, V);
  Matrix_Free(Vinv);  
  B1 = Matrix_Alloc(M1->NbRows, U->NbColumns);
  B2 = Matrix_Alloc(M2->NbRows, V->NbColumns);  
  Matrix_Product(M1, U, B1);
  Matrix_Product(M2, V, B2);
  Matrix_Free(M1);
  Matrix_Free(M2);  
  value_division(k,B2->p[0][0],B1->p[0][0]);
  value_division(k,DiagMatrix->p[0][0],k);  
  for (i = 0; i < DiagMatrix->NbRows; i++)
    value_division(DiagMatrix->p[i][i],DiagMatrix->p[i][i],k);
  newB1 = ChangeLatticeDimension(B1, B1->NbRows + 1); 
  Matrix_Free(B1);
  newB2 = ChangeLatticeDimension(B2, B2->NbRows +1);
  Matrix_Free(B2);  
  for(i = 0; i < newB1->NbRows - 1;i ++)
    value_assign(newB2->p[i][newB1->NbRows-1],Intersection->p[i][X->NbRows-1]);
  Head = SplitLattice(newB1,newB2,DiagMatrix); 
  Matrix_Free(newB1);
  Matrix_Free(DiagMatrix); 
  value_clear(k);
  return Head;
}



/**

***        Method : 
***               
**/
/* 
 * Return the Union of lattices that constitute the difference the lattices 
 * 'A' and 'B'. The dimensions of 'A' and 'B' should be the same. 
 * Note :
 *        Inorder to Find the Difference of Lattices, we make use of
 *        the following facts.
 *               
 * Theorem : Given Two Lattices L1 and L2, (L2 subset of L1) there exists a
 *           Basis B = {b1, b2,..bn} of L1 and integers {a1, a2...,an} such 
 *           that a1 divides a2, a2 divides a3 and so on and {a1b1, a2b2 ,...,
 *           .., anbn} is a Basis of L2. So given this theorem we can express 
 *           the Lattice L1 in terms of Union of Lattices Involving L2, such 
 *           that Lattice L1 = B1 = Union of (B2 + i1b1 + i2b2 + .. inbn) such
 *           that 0 <= i1 < a1; 0 <= i2 < a2; .......   0 <= in < an. We also
 *           know that A/B = A/(A Intersection B) and that (A Intersection B) 
 *           is a subset of A. So, Making use of these two facts, we find the 
 *           A/B. We Split The Lattice A in terms of Lattice (A Int B). From 
 *           this Union of Lattices Delete the Lattice (A Int B).
 *
 * Algorithm : 
 *
 *       Step 1:  Find Intersection = LatticeIntersection (A, B).
 *       Step 2:  Extract the Linear Parts of the Lattices A and Intersection.
 *                (while dealing with Basis we only deal with the Linear Parts)
 *       Step 3:  Let M1 = Basis of A and M2 = Basis of B.
 *                Let B1 and B2 be the Basis of A and B respectively, 
 *                corresponding to the above Theorem.
 *                Then we Have B1 = M1 * U1 {a unimodular Matrix }
 *                and B2 = M2 * U2. M1 and M2 we know, they are the linear 
 *                parts we obtained in Step 2. Our Task is now to find U1 and
 *                U2. 
 *                We know that B1  * Delta = B2.
 *                i.e. M1 * U1 * Delta = M2 * U2
 *                or U1*Delta*U2Inverse = M1Inverse * M2.
 *                and Delta is the Diagonal Matrix which satisifies the 
 *                above properties (in the Theorem).
 *                So Delta is nothing but the Smith Normal Form of 
 *                M1Inverse * M2.
 *                So, first we have to find M1Inverse.
 *             
 *                This Step, involves finding the Inverse of the Matrix M1.
 *                We find the Inverse using the Polylib function 
 *                Matrix_Inverse. There is a catch here, the result of this
 *                function is an integral matrix, not necessarily the exact
 *                Inverse (since M1 need not be Unimodular), but a multiple
 *                of the actual inverse. The number by which we have to divide
 *                the matrix, is not obtained here as the input matrix is not
 *                a Polylib matrix { We input only the Linear part }. Later I
 *                give a way for finding that number.
 *
 *                M1Inverse = Matrix_Inverse ( M1 );
 *      
 *      Step 4 :  MtProduct = Matrix_Product (M1Inverse, M2);
 *      Step 5 :  SmithNormalFrom (MtProduct, Delta, U, V);
 *                U1 = U and U2Inverse = V.
 *      Step 6 :  Find U2 = Matrix_Inverse  (U2inverse). Here there is no prob
 *                as U1 and its inverse are unimodular.
 *      
 *      Step 7 :  Compute B1 = M1 * U1;
 *      Step 8 :  Compute B2 = M2 * U2;
 *      Step 9 :  Earlier when we computed M1Inverse, we knew that it was not
 *                the exact inverse but a multiple of it. Now we find the 
 *                number, such that ( M1Inverse / number ) would give us the 
 *                exact inverse of M1.
 *                We know that B1 * Delta = B2.
 *                Let k = B2[0][0] / B1[0][0].
 *                Let number = Delta[0][0]/k;
 *                This 'number' is the number we want.
 *                We Divide the matrix Delta by this number, to get the actual
 *                Delta such that B1 * Delta = B2.
 *     Step 10 :  Call Split Lattice (B1, B2, Delta ).
 *                This function returns the Union of Lattices in such a way 
 *                that B2 is at the Head of this List.
 *     Step 11 :  To Remove B2 From the list of the Union of Lattices.
 *                Head = Head->next;
 *     Step 12 :  Free the Memory that is now not needed and return Head.
 *
 */
LatticeUnion *LatticeDifference(Lattice  *A,Lattice *B) {
 
  Lattice *Intersection = NULL;
  LatticeUnion *Head = NULL, *tempHead = NULL;
  Matrix *H , *U1 , *X, *Y ;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered LATTICEDIFFERENCE \n"); 
  fclose(fp);
#endif

  if (A->NbRows != A->NbColumns) { 
    fprintf(stderr,"\nIn LatticeDifference : The Input Matrix A is not a proper Lattice \n");
    return NULL;
  }
  
  if (B->NbRows != B->NbColumns) { 
    fprintf(stderr,"\nIn LatticeDifference : The Input Matrix B is not a proper Lattice \n");
    return NULL;
  }
  
  if (A->NbRows != B->NbRows) {
    fprintf(stderr,"\nIn Lattice Difference : The Input Lattices A and B have ");
    fprintf(stderr,"incompatible dimensions \n");
    return NULL;
  }

if (isinHnf (A) != True) {
    AffineHermite(A,&H,&U1);
    X = Matrix_Copy(H);    
    Matrix_Free(U1);
    Matrix_Free(H);
  }
  else
    X = Matrix_Copy(A);
  
  if (isinHnf(B) != True) {
    AffineHermite(B,&H,&U1);
    Y = Matrix_Copy(H);   
    Matrix_Free(H);
    Matrix_Free(U1);
  }
  else
    Y = Matrix_Copy(B);  
  if (isEmptyLattice(X)) {
    return NULL;  
  } 

  Head=Lattice2LatticeUnion(X,Y);

/* If the spliting operation can't be done the result is X. */

  if (Head == NULL) {
    Head = (LatticeUnion *)malloc(sizeof(LatticeUnion));
    Head->M = Matrix_Copy(X);
    Head->next = NULL;
    Matrix_Free(X);
    Matrix_Free(Y);
    return Head;
  } 

  tempHead = Head;
  Head = Head->next;  
  Matrix_Free (tempHead->M);
  tempHead->next = NULL; 
  free(tempHead);  

  if ((Head != NULL))
    Head = LatticeSimplify (Head);
  Matrix_Free (X);
  Matrix_Free (Y); 

  return Head;
} /* LatticeDifference */


/*
 * Given a Lattice 'B1' and a Lattice 'B2' and a Diagonal Matrix 'C' such that
 * 'B2' is a subset of 'B1' and C[0][0] divides C[1][1], C[1][1] divides C[2]
 * [2] and so on, output the list of matrices whose union is B1. The function
 * expresses the Lattice B1 in terms of B2 Unions of B1 = Union of {B2 + i0b0 +
 * i1b1 + .... + inbn} where 0 <= i0 < C[0][0]; 0 <= i1 < C[1][1] and so on and
 * {b0 ... bn} are the columns of Lattice B1. The list is so formed that the 
 * Lattice B2 is the Head of the list. 
 */     
LatticeUnion *SplitLattice(Lattice *B1, Lattice *B2, Matrix *C) {
  
  int i;
  
  LatticeUnion *Head = NULL;  
  Head = (LatticeUnion *)malloc(sizeof(LatticeUnion));
  Head->M = (Lattice *)B2;
  Head->next = NULL;  
  for (i = 0; i < C->NbRows ; i++)
    AddLattice(Head,B1,B2,VALUE_TO_INT(C->p[i][i]),i);  
  return Head;
} /* SplitLattice */

/*
 * Given lattices 'B1' and 'B2', an integer 'NumofTimes', a column number 
 * 'Colnumber' and a pointer to a list of lattices, the function does the 
 * following :-
 * For every lattice in the list, it adds a set of lattices such that the 
 * affine part of the new lattices is greater than the original lattice by 0 to
 * NumofTimes-1 * {the (ColumnNumber)-th column of B1}. 
 * Note : 
 * Three pointers are defined to point at various points of the list. They are:
 * Head   -> It always points to the head of the list. 
 * tail   -> It always points to the last element in the list. 
 * marker -> It points to the element, which is the last element of the Input 
 *           list.
 */ 
static void AddLattice (LatticeUnion *Head, Matrix  *B1,  Matrix *B2, int NumofTimes, int Colnumber) {
  
  LatticeUnion *temp, *tail, *marker;
  int i,j;
  Value tmp;
  
  value_init(tmp);
  tail =  Head;  
  while (tail->next != NULL)
    tail = tail->next;  
  marker = tail;
  
  for(temp = Head; temp != NULL; temp=temp->next) {
    for (i = 1; i < NumofTimes; i++) { 
      Lattice *tempMatrix, *H, *U;
      
      tempMatrix = (Lattice *)Matrix_Copy(temp->M);	  
      for (j = 0; j < B2->NbRows; j++) {
	value_set_si(tmp,i);
	value_addmul(tempMatrix->p[j][B2->NbColumns-1], tmp, B1->p[j][Colnumber]);
      }
      tail->next = (LatticeUnion *)malloc(sizeof(LatticeUnion)); 
      AffineHermite(tempMatrix,&H,&U);
      Matrix_Free((Matrix *)tempMatrix);
      Matrix_Free(U);
      tail->next->M = H;
      tail->next->next=NULL;
      tail = tail->next;
    }
    if (temp == marker)
      break;
  }  
  value_clear(tmp);
  return;
} /* AddLattice */

/* 
 * Given a polyhedron 'A', store the Hermite basis 'B' and return the true 
 * dimension of the polyhedron 'A'. 
 * Algorithm : 
 *
 *             1) First we find all the vertices of the Polyhedron A.
 *                Now suppose the vertices are [v1, v2...vn], then 
 *                a particular set of vectors governing the space of A are 
 *                given by [v1-v2, v1-v3, ... v1-vn] (let us say V).
 *                So we initially calculate these vectors.
 *             2) Then there are the rays and lines which contribute to the
 *                space in which A is going to lie.
 *                So we append to the rays and lines. So now we get a matrix
 *                {These are the rows} [ V ] [l1] [l2]...[lk] 
 *                where l1 to lk are either rays or lines of the Polyhedron A.
 *             3) The above matrix is the set of vectors which determine
 *                the space in which A is going to lie.
 *                Using this matrix we find a Basis which is such that
 *                the first 'm' columns of it determine the space of A. 
 *             4) But we also have to ensure that in the last 'n-m' 
 *                coordinates the Polyhedron is '0', this is done by
 *                taking the image by B(inv) of A and finding the remaining
 *                equalities, and composing it with the matrix B, so as
 *                to get a new matrix which is the actual Hermite Basis of
 *                the Polyhedron.
 */
int FindHermiteBasisofDomain(Polyhedron *A, Matrix **B) {
  
  int i, j;
  Matrix *temp,*temp1, *tempinv, *Newmat ;
  Matrix *vert, *rays, *result;
  Polyhedron *Image;
  int rank, equcount ;
  int noofvertices = 0, noofrays = 0;
  int vercount , raycount;
  Value lcm, fact;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered FINDHERMITEBASISOFDOMAIN \n"); 
  fclose(fp);
#endif
  
  POL_ENSURE_FACETS(A);
  POL_ENSURE_VERTICES(A);

  /* Checking is empty */  
  if (emptyQ(A)) {
      B[0] = Identity(A->Dimension+1);
      return(-1);
   }

  value_init(lcm); value_init(fact);
  value_set_si(lcm,1);

  /* Finding the Vertices */  
  for (i = 0; i < A->NbRays; i++)
    if ((value_notzero_p(A->Ray[i][0])) && value_notzero_p(A->Ray[i][A->Dimension+1]))
      noofvertices++;
    else 
      noofrays ++;
  
  vert = Matrix_Alloc(noofvertices,A->Dimension+1);
  rays = Matrix_Alloc(noofrays,A->Dimension);
  vercount = 0;
  raycount = 0;
  
  for(i = 0; i < A->NbRays; i++) {
    if ((value_notzero_p(A->Ray[i][0])) && value_notzero_p(A->Ray[i][A->Dimension+1])) {
      for(j = 1; j < A->Dimension+2; j++) 
	value_assign(vert->p[vercount][j-1],A->Ray[i][j]);
      Lcm3(lcm, A->Ray[i][j-1], &lcm);
      vercount++;
    }
    else {
      for (j = 1; j < A->Dimension+1; j++)
	value_assign(rays->p[raycount][j-1],A->Ray[i][j]);
      raycount++;	
    }
  }
  
  /* Multiplying the rows by the lcm */
  for(i = 0; i < vert->NbRows; i ++) {    
    value_division(fact,lcm,vert->p[i][vert->NbColumns-1]);
    for (j = 0; j < vert->NbColumns-1; j++)
      value_multiply(vert->p[i][j],vert->p[i][j],fact);
  }
  
  /* Drop the Last Columns */
  temp = RemoveColumn(vert,vert->NbColumns-1);
  Matrix_Free(vert);
  
  /* Getting the Vectors */
  vert = Matrix_Alloc(temp->NbRows-1, temp->NbColumns);
  for (i = 1; i < temp->NbRows; i++)
    for (j = 0; j < temp->NbColumns ; j++)
      value_subtract(vert->p[i-1][j],temp->p[0][j],temp->p[i][j]);

  Matrix_Free(temp);
  
  /* Add the Rays and Lines */
  /* Combined Matrix */  
  result = Matrix_Alloc(vert->NbRows+rays->NbRows, vert->NbColumns);
  for (i = 0; i < vert->NbRows; i++)
    for (j = 0 ;j < result->NbColumns ; j++)
      value_assign(result->p[i][j],vert->p[i][j]);
  
  for (; i<result->NbRows; i++)
    for (j = 0; j < result->NbColumns; j++)
      value_assign(result->p[i][j],rays->p[i-vert->NbRows][j]);

  Matrix_Free(vert);
  Matrix_Free(rays);

  rank = findHermiteBasis(result, &temp);
  temp1 = ChangeLatticeDimension(temp,temp->NbRows+1);

  Matrix_Free(result);
  Matrix_Free(temp);
  
  /* Adding the Affine Part to take care of the Equalities */  
  temp = Matrix_Copy(temp1);
  tempinv = Matrix_Alloc(temp->NbRows,temp->NbColumns);
  Matrix_Inverse(temp,tempinv); 
  Matrix_Free(temp);
  Image = DomainImage(A,tempinv,MAXNOOFRAYS);  
  Matrix_Free(tempinv);
  Newmat = Matrix_Alloc(temp1->NbRows,temp1->NbColumns);
  for(i = 0; i < rank ; i++)
    for(j = 0; j < Newmat->NbColumns ; j++) 
      value_set_si(Newmat->p[i][j],0);  
  for(i = 0; i < rank; i++)
    value_set_si(Newmat->p[i][i],1);  
  equcount = 0;  
  for (i = 0; i < Image->NbConstraints; i ++)
    if (value_zero_p(Image->Constraint[i][0])) {
      for (j = 1; j<Image->Dimension+2; j ++)
	value_assign(Newmat->p[rank+equcount][j-1],Image->Constraint[i][j]);
      ++equcount ;
    } 
  Domain_Free(Image);
  for (i = 0; i < Newmat->NbColumns-1; i++)
    value_set_si(Newmat->p[Newmat->NbRows-1][i],0);
  value_set_si(Newmat->p[Newmat->NbRows-1][Newmat->NbColumns-1],1);
  temp = Matrix_Alloc(Newmat->NbRows, Newmat->NbColumns);
  Matrix_Inverse(Newmat,temp);
  Matrix_Free(Newmat);
  B[0] = Matrix_Alloc(temp1->NbRows,temp->NbColumns);

  Matrix_Product(temp1,temp,B[0]);  
  Matrix_Free(temp1);
  Matrix_Free(temp);
  value_clear(lcm);
  value_clear(fact);
  return rank;
} /* FindHermiteBasisofDomain */ 

/*
 * Return the image of a lattice 'A' by the invertible, affine, rational 
 * function 'M'. 
 */
Lattice *LatticeImage(Lattice *A, Matrix *M) {
 
  Lattice *Img, *temp, *Minv;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered LATTICEIMAGE \n"); 
  fclose(fp);
#endif

  if ((A->NbRows != M->NbRows) || (M->NbRows != M->NbColumns))
    return (EmptyLattice (A->NbRows));
  
  if (value_one_p(M->p[M->NbRows-1][M->NbColumns-1])) {
    Img = Matrix_Alloc ( M->NbRows, A->NbColumns );
    Matrix_Product (M,A,Img);
    return Img;
  } 
  temp = Matrix_Copy(M);
  Minv = Matrix_Alloc(temp->NbColumns, temp->NbRows);
  Matrix_Inverse(temp, Minv);
  Matrix_Free(temp);
  
  Img = LatticePreimage(A, Minv);
  Matrix_Free (Minv);
  return Img;
} /* LatticeImage */

/* 
 * Return the preimage of a lattice 'L' by an affine, rational function 'G'.
 * Algorithm: 
 *           (1) Prepare Diophantine equation :
 *               [Gl -Ll][x y] = [Ga -La]{"l-linear, a-affine"}
 *           (2) Solve the Diophantine equations.
 *           (3) If there is solution to the Diophantine eq., extract the 
 *               general solution and the particular solution of x and that 
 *               forms the preimage of 'L' by 'G'. 
 */
Lattice *LatticePreimage(Lattice *L, Matrix *G) {
 
  Matrix *Dio, *U ;
  Lattice *Result;
  Vector *X;
  int i,j;
  int rank;
  Value divisor, tmp;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered LATTICEPREIMAGE \n"); 
  fclose(fp);
#endif
  
  /* Check for the validity of the function */  
  if (G->NbRows != L->NbRows) {
    fprintf (stderr, "\nIn LatticePreimage: Incompatible types of Lattice and the function\n");
    return (EmptyLattice(G->NbColumns));
  }
  
  value_init(divisor); value_init(tmp);
  
  /* Making Diophantine Equations [g -L] */
  value_assign(divisor,G->p[G->NbRows-1][G->NbColumns-1]);
  Dio = Matrix_Alloc(G->NbRows, G->NbColumns+L->NbColumns-1);
  for (i = 0; i < G->NbRows-1; i++)
    for (j = 0; j < G->NbColumns-1; j++)
      value_assign(Dio->p[i][j],G->p[i][j]);
  
  for (i = 0;i < G->NbRows-1; i++)
    for (j = 0; j < L->NbColumns-1; j++) {
      value_multiply(tmp,divisor,L->p[i][j]);
      value_oppose(Dio->p[i][j+G->NbColumns-1],tmp);
    }
  
  for (i = 0; i < Dio->NbRows-1; i++) {
    value_multiply(tmp,divisor,L->p[i][L->NbColumns-1]);
    value_subtract(tmp,G->p[i][G->NbColumns-1],tmp);
    value_assign(Dio->p[i][Dio->NbColumns-1],tmp);
  }
  for (i = 0; i < Dio->NbColumns-1; i++)
    value_set_si(Dio->p[Dio->NbRows-1][i],0);
  
  value_set_si(Dio->p[Dio->NbRows-1][Dio->NbColumns-1],1); 
  rank = SolveDiophantine(Dio, &U, &X);
  
  if (rank == -1)
    Result = EmptyLattice(G->NbColumns);
  else {
    Result = Matrix_Alloc (G->NbColumns, G->NbColumns);
    for (i = 0; i < Result->NbRows-1; i++)
      for (j = 0; j < Result->NbColumns-1; j++)
	value_assign(Result->p[i][j],U->p[i][j]);
    
    for (i = 0; i < Result->NbRows-1; i ++)
      value_assign(Result->p[i][Result->NbColumns-1],X->p[i]);  
    Matrix_Free (U);
    Vector_Free (X);   
    for (i = 0; i < Result->NbColumns-1; i ++)
      value_set_si(Result->p[Result->NbRows-1][i],0);
    value_set_si(Result->p[i][i],1);
  } 
  Matrix_Free(Dio);
  value_clear(divisor);
  value_clear(tmp);
  return Result;
} /* LatticePreimage */
 
/*
 * Return True if the matrix 'm' is a valid lattice, otherwise return False. 
 * Note: A valid lattice has the last row as [0 0 0 ... 1]. 
 */ 
Bool IsLattice(Matrix *m) {
  
  int i;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered ISLATTICE \n"); 
  fclose (fp);
#endif
  
  /* Is it necessary to check if the lattice
     is fulldimensional or not here only? */
  
  if (m->NbRows != m->NbColumns)
    return False;
  
  for (i = 0; i < m->NbColumns-1; i++)
    if (value_notzero_p(m->p[m->NbRows-1][i])) 
      return False ;
  if (value_notone_p(m->p[i][i])) 
    return False;
  return True ;
} /* IsLattice */ 
 
/*
 *  Check whether the matrix 'm' is full row-rank or not. 
 */ 
Bool isfulldim(Matrix *m) {
  
  Matrix *h, *u ;
  int i ;
  
  /* 
     res = Hermite (m, &h, &u);
     if (res != m->NbRows)
     return False ;
  */
  
  Hermite(m, &h, &u);
  for (i = 0; i < h->NbRows; i ++)
    if (value_zero_p(h->p[i][i])) {
      Matrix_Free (h);
      Matrix_Free (u);
      return False;
    }
  Matrix_Free (h);
  Matrix_Free (u);
  return True;
} /* isfulldim */
     
/* 
 * This function takes as input a lattice list in which the lattices have the 
 * same linear part, and almost the same affinepart, i.e. if A and B are two 
 * of the lattices in the above lattice list and [a1, .. , an] and [b1 .. bn]
 * are the affineparts of A and B respectively, then for 0 < i < n ai = bi and
 * 'an' may not be equal to 'bn'. These are not the affine parts in the n-th
 * dimension, but the lattices have been tranformed such that the value of the 
 * elment in the dimension on which we are simplifying is in the last row and
 * also the lattices are in a sorted order. 
 *              This function also takes as input the dimension along which we
 * are simplifying and takes the diagonal element of the lattice along that
 * dimension and tries to find out the factors of that element and sees if the
 * list of lattices can be simplified using these factors. The output of this 
 * function is the list of lattices in the simplified form and a flag to indic-
 * ate whether any form of simplification was actually done or not.  
 */
static Bool Simplify(LatticeUnion **InputList, LatticeUnion **ResultList, int dim) {
  
  int i;
  LatticeUnion *prev, *temp;
  factor allfac;
  Bool retval = False;
  int width;
  Value cnt, aux, k, fac, num, tmp, foobar;
  
  if ((*InputList == NULL) || (InputList[0]->next == NULL))
    return False ;

  value_init(aux); value_init(cnt); 
  value_init(k); value_init(fac);
  value_init(num); value_init(tmp);
  value_init(foobar);

  width = InputList[0]->M->NbRows-1; 
  allfac = allfactors(VALUE_TO_INT(InputList[0]->M->p[dim][dim]));
  value_set_si(cnt,0);
  for (temp = InputList[0]; temp != NULL; temp = temp->next)
    value_increment(cnt,cnt);  
  for(i = 0; i < allfac.count; i++) {
    value_set_si(foobar,allfac.fac[i]);
    value_division(aux,InputList[0]->M->p[dim][dim],foobar);
    if(value_ge(cnt,aux))
      break;
  }
  if (i == allfac.count) {
    value_clear(cnt); value_clear(aux);   
    value_clear(k); value_clear(fac);
    value_clear(num); value_clear(tmp);
    value_clear(foobar);
    return False; 
  }
  for (; i < allfac.count; i++) {    
    Bool Present = False;
    value_set_si(k,0);
    
    if (*InputList == NULL) {
      value_clear(cnt); value_clear(aux);
      value_clear(k); value_clear(fac);
      value_clear(num); value_clear(tmp);
      value_clear(foobar);
      return retval;
    }
    value_set_si(foobar,allfac.fac[i]);
    value_division(num,InputList[0]->M->p[dim][dim],foobar);
    while (value_lt(k,foobar)) {          
      Present = False;
      value_assign(fac,k); 
      for (temp = *InputList; temp != NULL; temp = temp->next) {
	if (value_eq(temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1],fac)) {
	  value_set_si(foobar,allfac.fac[i]);
	  value_addto(fac,fac,foobar);
	  if (value_ge(fac,(*InputList)->M->p[dim][dim])) {
	    Present = True;
	    break;
	  }
	}	
	if (value_gt(temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1],fac))
	  break;
      }      
      if (Present == True) {	
	retval = True;
	if (*ResultList == NULL)
	  *ResultList = temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
	else {
	  for (temp = *ResultList; temp->next != NULL; temp = temp->next);
	  temp->next = (LatticeUnion *) malloc (sizeof (LatticeUnion));
	  temp = temp->next;
	}
	temp->M = Matrix_Copy(InputList[0]->M); 
	temp->next = NULL;
	value_set_si(foobar,allfac.fac[i]);
	value_assign(temp->M->p[dim][dim],foobar);
	value_assign(temp->M->p[dim][width],k);
	value_set_si(temp->M->p[width][width],1);
	
	/* Deleting the Lattices from the curlist */
	value_assign(tmp,k);
	prev = NULL;
	temp = InputList[0];
	while (temp != NULL) {
	  if (value_eq(temp->M->p[width][width],tmp)) {
	    if (temp == InputList[0]) {
	      prev = temp;
	      temp = InputList [0] = temp->next;
	      Matrix_Free(prev->M);
	      free(prev);	  
	    }
	    else {
	      prev->next = temp->next;
	      Matrix_Free(temp->M);
	      free(temp);
	      temp = prev->next;
	    }
	    value_set_si(foobar,allfac.fac[i]);
	    value_addto(tmp,tmp,foobar); 
	  }
	  else {
	    prev = temp;
	    temp = temp->next;
	  }
	}   
      } 
      value_increment(k,k);
    }    
  }
  value_clear(cnt); value_clear(aux);
  value_clear(k); value_clear(fac);
  value_clear(num); value_clear(tmp);
  value_clear(foobar);
  return retval;          
} /* Simplify */     
 
/*
 * This function is used in the qsort function in sorting the lattices. Given 
 * two lattices 'A' and 'B', both in HNF, where A = [ [a11 0], [a21, a22, 0] .
 * .... [an1, .., ann] ] and B = [ [b11 0], [b21, b22, 0] ..[bn1, .., bnn] ],
 * then A < B, if there exists a pair <i,j> such that [aij < bij] and for every
 * other pair <i1, j1>, 0<=i1<i, 0<=j1<j [ai1j1 = bi1j1]. 
 */
static int LinearPartCompare(const void *A, const void *B) {
 
  Lattice **L1, **L2;
  int i, j;
  
  L1 = (Lattice **) A; 
  L2 = (Lattice **) B;
  
  for (i = 0;  i < L1[0]->NbRows-1; i++)
    for (j = 0; j <= i ; j++) {
      if (value_gt(L1[0]->p[i][j],L2[0]->p[i][j]))
	return 1;      
      if (value_lt(L1[0]->p[i][j],L2[0]->p[i][j]))
	return -1;
    }   
  return 0;
} /* LinearPartCompare */
  
/*
 * This function takes as input a List of Lattices and sorts them on the basis
 * of their Linear parts. It sorts in place, as a result of which the input 
 * list is modified to the sorted order.
 */ 
static void LinearPartSort (LatticeUnion *Head) {

  int  cnt;
  Lattice **Latlist;
  LatticeUnion *temp ;
  
  cnt = 0;
  for (temp = Head; temp != NULL; temp = temp->next)
    cnt ++;
 
  Latlist = (Lattice **) malloc ( sizeof (Lattice *) * cnt);
  
  cnt = 0; 
  for (temp = Head; temp != NULL; temp = temp->next)
    Latlist[cnt++] = temp->M; 
  
  qsort(Latlist, cnt, sizeof(Lattice *), LinearPartCompare);
  
  cnt = 0;  
  for (temp = Head; temp != NULL; temp = temp->next)
    temp->M = Latlist[cnt++];
  
 free (Latlist);
 return;
} /* LinearPartSort */

/*
 * This function is used in 'AfiinePartSort' in sorting the lattices with the
 * same linear part. GIven two lattices 'A' and 'B' with affineparts [a1 .. an]
 * and [b1 ... bn], then A < B if for some 0 < i <= n, ai < bi and for 0 < i1 <
 * i, ai1 = bi1. 
 */ 
static int AffinePartCompare(const void *A, const void *B) {
  
  int i;
  Lattice **L1, **L2; 

  L1 = (Lattice **)A;
  L2 = (Lattice **)B;
  
  for (i = 0; i < L1[0]->NbRows; i++) {
    if (value_gt(L1[0]->p[i][L1[0]->NbColumns-1],L2[0]->p[i][L1[0]->NbColumns-1]))
      return 1;  
    
    if (value_lt(L1[0]->p[i][L1[0]->NbColumns-1],L2[0]->p[i][L1[0]->NbColumns-1]))
      return -1;  
  }  
  return 0 ;
} /* AffinePartCompare */

/* 
 * This function takes a list of lattices with the same linear part and sorts
 * them on the basis of their affine part. The sorting is done in place.
 */
static void AffinePartSort (LatticeUnion *List) {
  
  int cnt;
  Lattice **LatList;
  LatticeUnion *tmp;
  
  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next)
    cnt ++;  
  
  LatList = (Lattice **) malloc (sizeof(Lattice *) * cnt);
  
  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next)
    LatList[cnt++] = tmp->M;
  
  qsort(LatList,cnt, sizeof (Lattice *), AffinePartCompare);
  
  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next) 
    tmp->M = LatList[cnt++];
  return;
} /* AffinePartSort */

static Bool AlmostSameAffinePart(LatticeUnion *A, LatticeUnion *B) {
  
  int i;
  
  if ((A == NULL) || (B == NULL))
    return False;

  for (i = 0; i < A->M->NbRows-1; i ++)
    if (value_ne(A->M->p[i][A->M->NbColumns-1],B->M->p[i][A->M->NbColumns-1]))
      return False;   
  return True;
} /* AlmostSameAffinePart */

/* 
 * This function takes a list of lattices having the same linear part and tries
 * to simplify these lattices. This may not be the only way of simplifying the
 * lattices. The function returns a list of partially simplified lattices and 
 * also a flag to tell whether any simplification was performed at all. 
 */
static Bool AffinePartSimplify(LatticeUnion *curlist, LatticeUnion **newlist) {
  
  int i;
  Value aux;
  LatticeUnion *temp, *curr, *next;
  LatticeUnion *nextlist;
  Bool change = False, chng;
  
  if (curlist == NULL) 
    return False;
  
  if (curlist->next == NULL) {
    curlist->next = newlist[0];
    newlist[0] = curlist;
    return False ;
  }
  
  value_init(aux);
  for (i = 0; i < curlist->M->NbRows - 1; i ++) {
    
    /* Interchanging the elements of the Affine part for easy computation
       of the sort (using qsort) */
    
    for (temp = curlist; temp != NULL; temp = temp->next) {      
      value_assign(aux,temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1]);
      value_assign(temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1],temp->M->p[i][temp->M->NbColumns-1]);
      value_assign(temp->M->p[i][temp->M->NbColumns-1],aux);
    }     
    AffinePartSort(curlist);
    nextlist = NULL;
    curr = curlist;
    while (curr != NULL) {
      next = curr->next;
      if (!AlmostSameAffinePart(curr, next)) {
	curr->next = NULL;
	chng = Simplify(&curlist, newlist, i);
	if (nextlist == NULL)
	  nextlist = curlist;
	else {
	  LatticeUnion *tmp;
	  for (tmp = nextlist; tmp->next; tmp=tmp->next);
	  tmp->next = curlist;
	}
	change = (Bool)(change | chng);
	curlist = next;
      }
      curr = next;
    }  
    curlist = nextlist;
    
    /* Interchanging the elements of the Affine part for easy computation
       of the sort (using qsort) */
    
    for(temp = curlist; temp != NULL; temp = temp->next) {
      value_assign(aux,temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1]);
      value_assign(temp->M->p[temp->M->NbRows-1][temp->M->NbColumns-1],temp->M->p[i][temp->M->NbColumns-1]);
      value_assign(temp->M->p[i][temp->M->NbColumns-1],aux);
    }     
    if (curlist == NULL)
      break;
  }
  if ( *newlist == NULL)
    *newlist = nextlist;
  else {
    for (curr = *newlist; curr->next != NULL; curr = curr->next);
    curr->next = nextlist;
  }  
  value_clear(aux);
  return change;
} /* AffinePartSimplify */
     
static Bool SameLinearPart(LatticeUnion *A, LatticeUnion *B) {
  
  int i, j;  
  if ((A == NULL) || (B ==NULL))
    return False;
  for (i = 0; i < A->M->NbRows-1; i++)
    for (j = 0; j <= i; j++)
      if (value_ne(A->M->p[i][j],B->M->p[i][j]))
	return False;
  
  return True;
} /* SameLinearPart */

/* 
 * Given a union of lattices, return a simplified list of lattices. 
 */ 
LatticeUnion *LatticeSimplify(LatticeUnion *latlist) {
  
  LatticeUnion  *curlist, *nextlist;
  LatticeUnion *curr, *next;
  Bool change = True, chng;
  
  curlist = latlist;
  while (change == True) {
    change = False;
    LinearPartSort(curlist);
    curr = curlist;
    nextlist = NULL;
    while(curr != NULL) {
      next = curr->next;
      if (!SameLinearPart(curr, next)) {
	curr->next = NULL; 
	chng = AffinePartSimplify(curlist, &nextlist);
	change = (Bool)(change | chng);
	curlist = next; 
      }
      curr = next; 
    }
    curlist = nextlist; 
  }
  return curlist;
} /* LatticeSimplify */

int intcompare (const void *a, const void *b) {

  int *i, *j;
  
  i = (int *) a;
  j = (int *) b;
  if (*i > *j)
    return 1;
  if (*i < *j)
    return -1;
  return 0;
} /* intcompare */

static int polylib_sqrt(int i);
static factor allfactors (int num) {
 
  int i,j, tmp;
  int noofelmts = 1;
  int *list, *newlist;
  int count;
  factor result;
  
  list = (int *)malloc(sizeof (int));
  list[0] = 1;
  
  tmp = num;
  for (i = 2; i <= polylib_sqrt(tmp); i++) {
    if ((tmp % i) == 0) {
      if (noofelmts == 0) {
	list = (int *) malloc (sizeof (int));
	list[0] = i;
	noofelmts = 1;
      }
      else {
	newlist = (int *) malloc (sizeof (int) * 2 * noofelmts + 1);
	for (j = 0; j < noofelmts; j++)
	  newlist[j] = list[j] ;
	newlist[j] = i;
	for (j = 0; j < noofelmts; j++)
	  newlist[j+noofelmts+1] = i * list[j];
	free (list);
	list = newlist;
	noofelmts= 2*noofelmts+1;
      }
      tmp = tmp / i;
      i = 1;
    } 
  }
  
  if ((tmp != 0) && (tmp != num)) {
    newlist = (int *) malloc (sizeof (int) * 2 * noofelmts + 1);
    for (j = 0; j < noofelmts; j ++)
      newlist[j] = list[j] ;
    newlist[j] = tmp;
    for (j = 0; j < noofelmts; j ++)
      newlist[j+noofelmts+1] = tmp * list[j];
    free (list);
    list = newlist;
    noofelmts= 2*noofelmts+1;
  }  
  qsort (list, noofelmts, sizeof(int), intcompare);
  count = 1;
  for (i = 1; i < noofelmts; i ++)
    if (list[i] != list[i-1])
      list[count++] = list[i]; 
  if (list[count-1] == num)
    count --;
  
  result.fac = (int *) malloc (sizeof (int) * count);
  result.count = count;
  for (i = 0; i < count; i ++)
    result.fac[i] = list[i];
  free (list); 
  return result;
} /* allfactors */

static int polylib_sqrt (int i) {
  
  int j;
  j = 0;
  i = i > 0 ? i : -i;
  
  while (1) {
    if ((j * j) > i)
      break;
    else
      j ++;
  }
  return (j-1);
} /* polylib_sqrt */







