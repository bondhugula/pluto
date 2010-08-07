

for (i=0; i<=N-1; i++)
  for (j=0; j<=i; j++)
    A[i][j] = A[i][j] + alpha*x[i]*x[j];
