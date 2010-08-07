

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    B[i][j] = A[i][j] + u1[i]*v1[j] + u2[i]*v2[j];

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    x[i] = x[i] + beta* B[j][i]*y[j];


for (i=0; i<=N-1; i++)
  x[i] = x[i] + z[i];

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    w[i] = w[i] + alpha* B[i][j]*x[j];

