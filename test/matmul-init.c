#pragma scop
for (i = 0; i < n; i++)
  for (j = 0; j < n; j++) {
    C[i][j] = 0;
    for (k = 0; k < n; k++)
      C[i][j] = C[i][j] + A[i][k] * B[k][j];
  }
#pragma endscop
