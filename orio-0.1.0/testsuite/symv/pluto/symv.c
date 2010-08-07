

/* pluto start (N,alpha,beta) */

for (j=0; j<=N-1; j++)
  for (i=0; i<=j-1; i++)
    y[i] = beta*y[i] + alpha*A[i][j]*x[j];

for (i=0; i<=N-1; i++)
  for (j=i; j<=N-1; j++)
    y[i] = beta*y[i] + alpha*A[i][j]*x[j];

/* pluto end */
