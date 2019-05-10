#pragma scop
for (i = 0; i < n; i++)
  A[i] = ina[i];

for (i = 0; i < n; i++)
  C[i] = inc[i];

for (i = 0; i < n; i++)
  for (j = 0; j < n; j++)
    C[i] = C[i] + A[j];
#pragma endscop
