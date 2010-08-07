#include <stdlib.h>
#include <polylib/polylib.h>

static void RearrangeMatforSolveDio(Matrix *M);

/*
 *  Solve Diophantine Equations :
 *        This function takes as input a system of equations in the form
 *          Ax + C = 0 and finds the solution for it, if it exists
 *       
 *        Input : The matrix form the system of the equations Ax + C = 0
 *                 ( a pointer to a Matrix. )
 *                 A pointer to the pointer, where the matrix U 
 *                  corresponding to the free variables of the equation
 *                  is stored.
 *                 A pointer to the pointer of a vector is a solution to T.
 *
 *
 *        Output : The above matrix U and the vector T.
 *
 *        Algorithm :
 *                    Given an integral matrix A, we can split it such that
 *                    A = HU, where H is in HNF (lowr triangular)
 *                      and U is unimodular.
 *                    So Ax = c -> HUx = c -> Ht = c ( where Ux = t).
 *                       Solving for Ht = c is easy.
 *                       Using 't' we find x = U(inverse) * t.
 * 
 *        Steps :
 *                   1) For the above algorithm to work correctly to
 *                      need the condition that the first 'rank' rows are
 *                      the rows which contribute to the rank of the matrix.
 *                      So first we copy Input into a matrix 'A' and 
 *                      rearrange the rows of A (if required) such that
 *                      the first rank rows contribute to the rank.
 *                   2) Extract A and C from the matrix 'A'. A = n * l matrix.
 *                   3) Find the Hermite normal form of the matrix A.
 *                       ( the matrices the lower tri. H and the unimod U).
 *                   4) Using H, find the values of T one by one.
 *                      Here we use a sort of Gaussian elimination to find
 *                      the solution. You have a lower triangular matrix
 *                      and a vector, 
 *                      [ [a11, 0], [a21, a22, 0] ...,[arank1...a rankrank 0]]
 *                       and the solution vector [t1.. tn] and the vector 
 *                      [ c1, c2 .. cl], now as we are traversing down the
 *                      rows one by one, we will have all the information 
 *                      needed to calculate the next 't'.
 *                      
 *                      That is to say, when you want to calculate t2, 
 *                      you would have already calculated the value of t1
 *                      and similarly if you are calculating t3, you will 
 *                      need t1 and t2 which will be available by that time.
 *                      So, we apply a sort of Gaussian Elimination inorder
 *                      to find the vector T.
 *
 *                   5) After finding t_rank, the remaining (l-rank) t's are
 *                      made equal to zero, and we verify, if these values
 *                      agree with the remaining (n-rank) rows of A.
 *
 *                   6) If a solution exists, find the values of X using 
 *                        U (inverse) * T.
 */

int SolveDiophantine(Matrix *M, Matrix **U, Vector **X) {
  
  int i, j, k1, k2, min, rank;
  Matrix *A, *temp, *hermi, *unimod,  *unimodinv ;
  Value *C; /* temp storage for the vector C */
  Value *T; /* storage for the vector t */
  Value sum, tmp;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered SOLVEDIOPHANTINE\n"); 
  fclose(fp);
#endif

  value_init(sum); value_init(tmp);
  
  /* Ensuring that the first rank row of A contribute to the rank*/ 
  A = Matrix_Copy(M);
  RearrangeMatforSolveDio(A);
  temp = Matrix_Alloc(A->NbRows-1, A->NbColumns-1);
  
  /* Copying A into temp, ignoring the Homogeneous part */ 
  for (i = 0; i < A->NbRows -1; i++)
    for (j = 0; j < A->NbColumns-1; j++)
      value_assign(temp->p[i][j],A->p[i][j]);
  
  /* Copying C into a temp, ignoring the Homogeneous part */ 
  C = (Value *) malloc (sizeof(Value) * (A->NbRows-1));
  k1 = A->NbRows-1;
  
  for (i = 0; i < k1; i++) {
    value_init(C[i]);
    value_oppose(C[i],A->p[i][A->NbColumns-1]);
  }
  Matrix_Free (A); 
  
  /* Finding the HNF of temp */  
  Hermite(temp, &hermi, &unimod);
  
  /* Testing for existence of a Solution */
  
  min=(hermi->NbRows <= hermi->NbColumns ) ? hermi->NbRows : hermi->NbColumns ;
  rank = 0;
  for (i = 0; i < min ; i++) {
    if (value_notzero_p(hermi->p[i][i]))
      rank ++;
    else
      break ;
  }
  
  /* Solving the Equation using Gaussian Elimination*/
  
  T = (Value *) malloc(sizeof(Value) * temp->NbColumns);
  k2 = temp->NbColumns;
  for(i=0;i< k2; i++) 
    value_init(T[i]);

  for (i = 0; i < rank ; i++) {
    value_set_si(sum,0);
    for (j = 0; j < i; j++) {
      value_addmul(sum, T[j], hermi->p[i][j]);
    } 
    value_subtract(tmp,C[i],sum);
    value_modulus(tmp,tmp,hermi->p[i][i]);
    if (value_notzero_p(tmp)) { /* no solution to the equation */
      *U = Matrix_Alloc(0,0);
      *X = Vector_Alloc (0);
      value_clear(sum); value_clear(tmp);
      for (i = 0; i < k1; i++) 
	value_clear(C[i]);
      for (i = 0; i < k2; i++) 
	value_clear(T[i]);
      free(C);
      free(T);
      return (-1);
    };
    value_subtract(tmp,C[i],sum);
    value_division(T[i],tmp,hermi->p[i][i]);
  }
  
  /** Case when rank < Number of Columns; **/
  
  for (i = rank; i < hermi->NbColumns; i++)
    value_set_si(T[i],0);
  
  /** Solved the equtions **/
  /** When rank < hermi->NbRows; Verifying whether the solution agrees 
      with the remaining n-rank rows as well. **/
  
  for (i = rank; i < hermi->NbRows; i++) {
    value_set_si(sum,0);
    for (j = 0; j < hermi->NbColumns; j++) {
      value_addmul(sum, T[j], hermi->p[i][j]);
    }  
    if (value_ne(sum,C[i])) {
      *U = Matrix_Alloc(0,0);
      *X = Vector_Alloc (0);
      value_clear(sum); value_clear(tmp);
      for (i = 0; i < k1; i++) 
	value_clear(C[i]);
      for (i = 0; i < k2; i++) 
	value_clear(T[i]);
      free(C);
      free(T);
      return (-1);
    }
  }     
  unimodinv = Matrix_Alloc(unimod->NbRows, unimod->NbColumns);
  Matrix_Inverse(unimod, unimodinv);
  Matrix_Free(unimod);
  *X = Vector_Alloc(M->NbColumns-1);
  
  if (rank == hermi->NbColumns)
    *U = Matrix_Alloc(0,0);
  else { /* Extracting the General solution form U(inverse) */
    
    *U = Matrix_Alloc(hermi->NbColumns, hermi->NbColumns-rank);   
    for (i = 0; i < U[0]->NbRows; i++)
      for (j = 0; j < U[0]->NbColumns; j++)
	value_assign(U[0]->p[i][j],unimodinv->p[i][j+rank]);
  }
  
  for (i = 0; i < unimodinv->NbRows; i++) { 
    
    /* Calculating the vector X = Uinv * T */
    value_set_si(sum,0);
    for (j = 0; j < unimodinv->NbColumns; j++) {
      value_addmul(sum, unimodinv->p[i][j], T[j]);
    }  
    value_assign(X[0]->p[i],sum);
  }
  
  /*
    for (i = rank; i < A->NbColumns; i ++)
    X[0]->p[i] = 0;
  */
  Matrix_Free (unimodinv);
  Matrix_Free (hermi);
  Matrix_Free (temp);
  value_clear(sum); value_clear(tmp);
  for (i = 0; i < k1; i++) 
    value_clear(C[i]);
  for (i = 0; i < k2; i++) 
    value_clear(T[i]);
  free(C);
  free(T);
  return (rank);  
} /* SolveDiophantine */

/*
 * Rearrange :
 *            This function takes as input a matrix M (pointer to it)
 *            and it returns the tranformed matrix M, such that the first
 *            'rank' rows of the new matrix M are the ones which contribute
 *            to the rank of the matrix M. 
 *          
 *            1) For a start we try to put all the zero rows at the end.
 *            2) Then cur = 1st row of the remaining matrix.
 *            3) nextrow = 2ndrow of M.
 *            4) temp = cur + nextrow
 *            5) If (rank(temp) == temp->NbRows.) {cur = temp;nextrow ++}
 *            6) Else (Exchange the nextrow of M with the currentlastrow.
 *                     and currentlastrow --).
 *            7) Repeat steps 4,5,6 till it is no longer possible.
 *            
 */
static void RearrangeMatforSolveDio(Matrix *M) {
  
  int i, j, curend, curRow, min, rank=1;
  Bool add = True;
  Matrix *A, *L, *H, *U;
  
  /* Though I could have used the Lattice function
	Extract Linear Part, I chose not to use it so that
	this function can be independent of Lattice Operations */

  L = Matrix_Alloc(M->NbRows-1,M->NbColumns-1);
  for (i = 0; i < L->NbRows; i++)
    for (j = 0; j < L->NbColumns; j++)
      value_assign(L->p[i][j],M->p[i][j]);
  
  /* Putting the zero rows at the end */
  curend = L->NbRows-1;
  for (i = 0; i < curend; i++) {
    for (j = 0; j < L->NbColumns; j++) 
      if (value_notzero_p(L->p[i][j]))
	break;
    if (j == L->NbColumns) {
      ExchangeRows(M,i,curend);
      curend --;
    }
  } 
  
  /* Trying to put the redundant rows at the end */
  
  if (curend > 0) { /* there are some useful rows */
    
    Matrix *temp;
    A = Matrix_Alloc(1, L->NbColumns); 
    
    for (i = 0; i <L->NbColumns; i++) 
      value_assign(A->p[0][i],L->p[0][i]);     
    curRow = 1;    
    while (add == True ) {
      temp= AddANullRow(A);
      for (i = 0;i <A->NbColumns; i++)
	value_assign(temp->p[curRow][i],L->p[curRow][i]);      
      Hermite(temp, &H, &U);
      for (i = 0; i < H->NbRows; i++)
	if (value_zero_p(H->p[i][i]))
	  break;
      if (i != H->NbRows) {
	ExchangeRows(M, curRow, curend);
	curend --;
      }
      else {
	curRow ++;
	rank ++;
	Matrix_Free (A);
	A = Matrix_Copy (temp);
	Matrix_Free (temp);
      }      
      Matrix_Free (H);
      Matrix_Free (U);
      min = (curend >= L->NbColumns) ? L->NbColumns : curend ;
      if (rank==min || curRow >= curend)
	break;
    }
    Matrix_Free (A);
  }
  Matrix_Free (L);
  return;
} /* RearrangeMatforSolveDio */
