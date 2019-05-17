// CHECK: Output written
#pragma scop
for (i = 0; i < N; i++)
  for (j = 0; j < N; j++) {
    a[i] = 0;
    a[j] = 0;
  }
#pragma endscop
