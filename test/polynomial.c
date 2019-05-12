// Hyperplane i+j is parallel.
// CHECK: T(S2): (1, i+j, i)
#pragma scop
for (i = 0; i < 2 * n; i++) {
  c[i] = 0;
}

for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    c[i + j] = c[i + j] + a[i] * b[j];
  }
}
#pragma endscop
