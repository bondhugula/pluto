// CHECK: Output written

#pragma scop
for (i = 0; i < N; i++) {
  y[i] = 0;
  for (j = 0; j < N; j++) {
    y[i] = y[i] + a[i][j] * x[j];
  }
}
#pragma endscop
