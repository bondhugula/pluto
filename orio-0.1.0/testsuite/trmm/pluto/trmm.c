
/* pluto start (N,alpha) */

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    for (k=i+1; k<=N-1; k++)
      B[i][j] = B[i][j] + alpha*A[i][k]*B[k][j];

/* pluto end */

