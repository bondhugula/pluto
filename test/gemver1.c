constant N, alpha, beta;

for (i = 0; i < N; i++)
  for (j = 0; j < N; j++)
    A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

for (i = 0; i < N; i++)
  for (j = 0; j < N; j++)
    x[i] = x[i] + beta * A[j][i] * y[j];
