
/* pluto start (N,alpha) */

for (i=0; i<=N-1; i++)
    for (j=i; j<=N-1; j++)
        A[i][j] = A[i][j] +  alpha*x[i]*x[j];

/* pluto end */
