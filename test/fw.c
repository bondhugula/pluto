
#pragma scop
for (k = 1; k <= n; k++)
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      a[i][j] = a[i][j] - a[i][k] * a[k][j];
#pragma endscop
