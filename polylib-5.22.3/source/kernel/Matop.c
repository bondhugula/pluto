#include <stdlib.h>
#include <polylib/polylib.h>

/* computes c = lcm(a,b) using Gcd(a,b,&c) */
void Lcm3(Value a, Value b, Value *c)
{
  Value tmp;

  if (value_zero_p(a)) {
    value_assign(*c, b);
    return;
  }
  if (value_zero_p(b)) {
    value_assign(*c, a);
    return;
  }
  value_init(tmp);
  value_multiply(tmp, a, b);
  value_absolute(tmp, tmp);
  Gcd(a,b,c);
  value_division(*c, tmp, *c);
  value_clear(tmp);
}

/*
 * Return the lcm of 'i' and 'j' 
 */ 
Value *Lcm (Value i, Value j)
{
  Value *tmp;
  
  tmp = (Value *) malloc (sizeof(Value));
  value_init(*tmp);
  Lcm3(i, j, tmp);
  return tmp;
} /* Lcm */

/* 
 * Return an identity matrix of size 'size'. 
 */
Matrix *Identity(unsigned size)
{
  unsigned i;
  Matrix *A;
  
  A = Matrix_Alloc(size, size);
  for (i = 0; i < size; i++)  
    value_set_si(A->p[i][i], 1);
  return A;
} /* Identity */

/*
 * Exchange the rows 'Row1' and 'Row2' of the matrix 'M'
 */
void ExchangeRows(Matrix *M, int Row1, int Row2) {
  
  int i;
  Value temp;
  
  value_init(temp);
  for (i = 0; i < (int)M->NbColumns;  i++) {
    value_assign(temp,M->p[Row1][i]);
    value_assign(M->p[Row1][i],M->p[Row2][i]);
    value_assign(M->p[Row2][i],temp);
  }  
  value_clear(temp);
  return ;
} /* ExchangeRows */

/*
 * Exchange the columns 'Column1' and 'Column2' of the matrix 'M'
 */
void ExchangeColumns(Matrix *M, int Column1, int Column2) {
  
  int i;
  
  for (i = 0; i < M->NbRows; i++)
    value_swap(M->p[i][Column1],M->p[i][Column2]);

  return;
} /* ExchangeColumns */

/* 
 * Return the Transpose of a matrix 'A' 
 */
Matrix *Transpose (Matrix *A) {
  
  int i,j;
  Matrix *transpose;
  
  transpose = Matrix_Alloc (A->NbColumns, A->NbRows);
  for (i = 0; i < (int)A->NbRows; i++)
    for (j = 0; j < (int)A->NbColumns; j++)
      value_assign(transpose->p[j][i],A->p[i][j]);
  return transpose;
} /* Transpose */

/*
 * Return a copy of the contents of a matrix 'Src'.
 */
Matrix *Matrix_Copy(Matrix const *Src )
{
  Matrix *Dst;
  unsigned i, j;
  
  Dst = Matrix_Alloc(Src->NbRows, Src->NbColumns);
  
  for (i = 0; i < Src->NbRows; i++)
    for (j = 0; j < Src->NbColumns; j++)
      value_assign(Dst->p[i][j],Src->p[i][j]);
  return Dst;
} /* Matrix_copy */

/* 
 * Test if the input matrix is integral or not. 
 */
Bool isIntegral (Matrix *A) {
  
  unsigned i, j;
  Value divisor, tmp;
  
  value_init(divisor); value_init(tmp);
  value_assign(divisor,A->p[A->NbRows-1][A->NbColumns-1]);

  for (i = 0; i < A->NbRows ; i++)
    for (j = 0; j < A->NbColumns ; j++) {
      value_modulus(tmp,A->p[i][j],divisor);
      if (value_notzero_p(tmp)) {
	value_clear(divisor); value_clear(tmp);
	  return False;  
      }	  
    }
  value_clear(divisor); value_clear(tmp);
  return True ;
} /* isIntegral */

/*
 * Check if the matrix 'A' is in Hermite normal form or not. 
 */
Bool isinHnf(Matrix *A) {
 
  Matrix *temp ;
  unsigned i, j ;
  Value rem;
  
  value_init(rem); 
  temp = Homogenise(A,True) ;
  for (i = 0; i < temp->NbRows; i++) {
    value_assign(rem,temp->p[i][i]);
    for (j = 0; j < i ; j++)
      if (value_ge(temp->p[i][j],rem)) {
	Matrix_Free(temp);
	value_clear(rem);
	return False ;
      }
    for (j = i+1; j < temp->NbColumns ; j++)
      if (value_notzero_p(temp->p[i][j])) {
	Matrix_Free(temp);
	value_clear(rem);
	return False ;
      }
  }
  value_clear(rem);
  return True ;
} /* isinHnf */

/*
 * Remove the row 'Rownumber' and place it at the end of the matrix 'X'  
 */ 
void PutRowLast (Matrix *X, int Rownumber) {
  
  int i, j ;
  Value temp;
 
  if (Rownumber == X->NbRows-1)
    return;
  
  value_init(temp);
  for (j = 0; j < X->NbColumns; j++) {
    value_assign(temp,X->p[Rownumber][j]);
    for (i = Rownumber; i < X->NbRows-1; i++)
      value_assign(X->p[i][j],X->p[i+1][j]);
    value_assign(X->p[i][j],temp);
  }
  value_clear(temp);
  return;
} /* PutRowLast */

/* 
 * Remove the row 'Rownumber' and place it at the begining of the matrix 'X' 
 */ 
void PutRowFirst(Matrix *X, int Rownumber) {
  
  int i, j ;
  Value temp;

  value_init(temp);
  for (j = 0; j < X->NbColumns; j++) {
    value_assign(temp,X->p[Rownumber][j]);
    for (i = Rownumber; i > 0; i--)
      value_assign(X->p[i][j],X->p[i-1][j]);
    value_assign(X->p[i][j],temp);
  }
  value_clear(temp);
  return;
} /* PutRowFirst */

/*
 * Remove the column 'Columnnumber' and place it at the begining of the matrix
 * 'X'
 */
void PutColumnFirst (Matrix *X, int Columnnumber) {
  
  int i, j ;
  Value temp;

  value_init(temp);
  for (i = 0; i < X->NbRows; i ++) {
    value_assign(temp,X->p[i][Columnnumber]);
    for (j = Columnnumber; j > 0; j --)
      value_assign(X->p[i][j],X->p[i][j-1]);
    value_assign(X->p[i][0],temp);
  }
  value_clear(temp);
  return ;
} /* PutColumnFirst */

/*
 * Remove the column 'Columnnumber' and place it at the end of the matrix 'X'
 */
void PutColumnLast (Matrix *X, int Columnnumber) {
  
  int i, j ;
  Value temp;

  value_init(temp);
  for (i = 0; i < X->NbRows; i++) {
    value_assign(temp,X->p[i][Columnnumber]);
    for (j = Columnnumber; j < X->NbColumns-1; j++)
      value_assign(X->p[i][j],X->p[i][j+1]);
    value_assign(X->p[i][X->NbColumns-1],temp);
  }
  value_clear(temp);
  return ;
} /* PutColumnLast */

/*
 * Add a row of zeros at the end of the matrix 'M' and return the new matrix 
 */
Matrix *AddANullRow (Matrix *M) {
  
  int i,j;
  Matrix *Result;
  
  Result = Matrix_Alloc(M->NbRows+1,M->NbColumns);
  for (i = 0;i < M->NbRows; i++)
    for (j = 0; j < M->NbColumns ; j++)
      value_assign(Result->p[i][j],M->p[i][j]);
  for (j = 0; j < M->NbColumns; j++)
    value_set_si(Result->p[i][j],0);  
  return Result;
} /* AddANullRow */

/* 
 * Add a column of zeros at the end of the matrix 'M' and return the new 
 * matrix
 */   
Matrix *AddANullColumn(Matrix *M) {
  
  int i,j;
  Matrix *Result;
  
  Result = Matrix_Alloc(M->NbRows, M->NbColumns+1);
  for (i = 0;i < M->NbRows; i++)
    for (j = 0; j < M->NbColumns ; j++)
      value_assign(Result->p[i][j],M->p[i][j]);
  for (i = 0; i < M->NbRows; i++)
    value_set_si(Result->p[i][M->NbColumns],0);
  return Result;
} /* AddANullColumn */

/*
 * Remove a row 'Rownumber' from matrix 'M' and return the new matrix.
 */
Matrix *RemoveRow(Matrix *M, int Rownumber) {
  
  Matrix *Result;
  int i;
  
  Result = Matrix_Alloc(M->NbRows-1, M->NbColumns);
  
  for (i = 0; i < Rownumber; i++)
    Vector_Copy(M->p[i], Result->p[i], M->NbColumns);
  for ( ; i < Result->NbRows; i++)
    Vector_Copy(M->p[i+1], Result->p[i], M->NbColumns);

  return Result;
} /* RemoveRow */

/*
 * Remove NumColumns columns starting at column number 'FirstColumnnumber' from matrix 'M'
 * and return the new matrix.
 */
Matrix *RemoveNColumns (Matrix *M, int FirstColumnnumber, int NumColumns) {
  
  Matrix *Result;
  int i;
  
  Result = Matrix_Alloc (M->NbRows, M->NbColumns-NumColumns);
  
  for (i = 0; i < Result->NbRows; i++) {
    Vector_Copy(M->p[i], Result->p[i], FirstColumnnumber);
    Vector_Copy(M->p[i]+FirstColumnnumber+NumColumns, Result->p[i]+FirstColumnnumber, 
		M->NbColumns-NumColumns-FirstColumnnumber);
  }
  return Result;
} /* RemoveColumn */

/*
 * Remove a column 'Columnnumber' from matrix 'M' and return the new matrix.
 */
Matrix *RemoveColumn (Matrix *M, int Columnnumber) {
  
  Matrix *Result;
  int i;
  
  Result = Matrix_Alloc (M->NbRows, M->NbColumns-1);
  
  for (i = 0; i < Result->NbRows; i++) {
    Vector_Copy(M->p[i], Result->p[i], Columnnumber);
    Vector_Copy(M->p[i]+Columnnumber+1, Result->p[i]+Columnnumber, 
		M->NbColumns-1-Columnnumber);
  }
  return Result;
} /* RemoveColumn */

/* 
 * Given a Matrix M of dimesnion n * l and rank l1, find a unimodular matrix
 * 'Result' such that the Vector Space spanned by M is the subset of the vector
 * Space spanned by the first l1 Rows of Result. The function returns the rank
 * l1 and the Matrix Result.  
 */
int findHermiteBasis(Matrix *M, Matrix **Result) {
  
  int i, j;
  Matrix *C, *curMat, *temp1, *temp2;
  Matrix *H, *U;
  Vector *V;
  int dim, curDim, curVect,rank;
  
  if (M->NbRows == 0) {
    Result[0] = Identity (M->NbColumns);
    return 0;
  }
  
  if (M->NbRows <= M->NbColumns) {
    Hermite(M, &H, &U);
    
    for (i = 0; i < H->NbRows; i++)
      if (value_zero_p(H->p[i][i]))
	break;
    
    if (i == H->NbRows) {
      Result[0] = Transpose(U);
      Matrix_Free(H);
      Matrix_Free(U);
      return(i);
    }
    Matrix_Free (H);
    Matrix_Free (U);
  } 
  
  /* Eliminating the Zero Rows  */
  
  C = Matrix_Copy (M);
  for (i = 0; i < C->NbRows; i++) {
    for (j = 0; j < C->NbColumns; j++)
      if (value_notzero_p(C->p[i][j]))
	break;
    if (j == C->NbColumns) {
      Matrix *temp;
      temp = RemoveRow(C, i); 
      Matrix_Free(C);
      C = Matrix_Copy(temp);
      Matrix_Free(temp);
      i --;
    }
  }
  
  /* Eliminating the Redundant Rows */
  
  curDim = 1;
  curVect = 1;
  dim = C->NbColumns;
  
  curMat = Matrix_Alloc(1,C->NbColumns);
  for (i = 0; i < C->NbColumns; i ++)
    value_assign(curMat->p[0][i],C->p[0][i]); 
  
  while((curVect < C->NbRows) && (curDim < dim)) {
    Matrix *temp;
    temp = AddANullRow(curMat);
    for (i = 0; i < C->NbColumns; i++)
      value_assign(temp->p[temp->NbRows-1][i],C->p[curVect][i]);
    
    temp1 = AddANullRow(temp);
    temp2 = AddANullColumn(temp1);
    rank = SolveDiophantine(temp2, &U, &V);
    if (rank == temp->NbRows) {
      Matrix_Free(curMat);
      curMat = Matrix_Copy(temp);
      curDim ++;
    }    
    curVect ++;    
    Matrix_Free (U);
    Vector_Free (V);
    Matrix_Free (temp1);
    Matrix_Free (temp);
    Matrix_Free (temp2);
  }
  Matrix_Free(C);
  
  Hermite(curMat, &H, &U);
  rank = curMat->NbRows;
  Matrix_Free(curMat);
  
  Result[0] = Transpose (U);
  Matrix_Free (H);
  Matrix_Free (U);
  return (rank);
} /* findHermiteBasis */




