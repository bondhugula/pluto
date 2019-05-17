
#pragma scop
for (i = 0; i < N; i++)
  for (j = 1; j < N; j++)
    a[i][j] = a[j][i] + a[i][i];
#pragma endscop
