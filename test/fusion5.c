// CHECK: Output written

#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    B[i][j] = 2 * A[i][j];
  }
}

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    C[i][j] = 3 * B[i - 1][N - j];
  }
}
#pragma endscop
