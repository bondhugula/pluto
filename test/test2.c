#pragma scop

for (i = 0; i < n; i++)
  A[i] = ina[i];

for (i = 0; i < n; i++)
  B[i] = inb[i];

for (i = 0; i < n; i++)
  for (j = 0; j < n; j++)
    B[i] = B[i] + A[j];

for (i = 0; i < n; i++)
  outb[i] = B[i];
#pragma endscop
