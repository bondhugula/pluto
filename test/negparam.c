// CHECK: Output written

#pragma scop
for (i = 4 - n; i <= n + 2; i++)
  for (j = 4 - n; j <= n + 2; j++)
    a[i][j] = a[i - 1][j] + 2;

for (i = 4 - n; i <= n + 2; i++)
  for (j = 4 - n; j <= n + 2; j++)
    a[i][j] = a[i][j] + 1;

#pragma endscop
