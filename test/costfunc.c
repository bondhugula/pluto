// This validates Pluto's cost function. (i+j) leads to a constant amount of
// communication/boundary misses, i.e., u = 0.
// CHECK: T(S1): (i+j, i)
#pragma scop

for (i = 0; i < N; i++) {
  for (j = 1; j < N; j++) {
    a[i][j] = a[j][i] + a[i][j - 1];
  }
}

#pragma endscop
