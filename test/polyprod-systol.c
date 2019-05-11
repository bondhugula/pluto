CONSTANT n;

DO i = 0, n DO j = 0,
   n A[i, j] = A[i, j + 1] B[i, j] = B[i - 1, j] C[i, j] =
       C[i - 1, j + 1] + A[i, j] * B[i, j] END DO END DO
