// CHECK: Output written
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    y[i] = y[i] + a[i][j] * x[j];
  }
}

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    z[i] = z[i] + b[i][j] * y[j];
  }
}
#pragma endscop
