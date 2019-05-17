CONSTANT n;

DO i = 0, n DO j = 0, n C[i, j] = C[i - 1, j + 1] + A[i] * B[j] END DO END DO
