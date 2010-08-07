
/* pluto start (N) */

for (i=0; i<=N-1; i++) 
  for (j=i+1; j<=N-1; j++) 
    X[i] = X[i] + A[i][j]*X[j]; 

/* pluto end */

