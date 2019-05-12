// CHECK: Output written
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    a[2 * i] = a[2 * i] + b[j];
  }
  for (j = 1; j < N; j++) {
    a[2 * i + 1] = a[2 * i + 1] + b[j];
  }
}
#pragma endscop
