

/* A is an UPPER triangular array that is a NOT-UNIT triangular */

for (i=0; i<=N-1; i++)
  {
    X[i] = A[i][i]*X[i];
    for (j=i+1; j<=N-1; j++)
      X[i] = X[i] + A[i][j]*X[j];
  }

/* A is an UPPER triangular array that is a UNIT triangular */

for (i=0; i<=N-1; i++)
  for (j=i+1; j<=N-1; j++)
    X[i] = X[i] + A[i][j]*X[j];


