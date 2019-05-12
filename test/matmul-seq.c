// CHECK: Output written

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    C[i][j] = 0;
  }
}

#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++) {
      C[i][j] = C[i][j] + A[i][k] * B[k][j];
    }
  }
}

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    for (k = 0; k < N; k++) {
      D[i][j] = D[i][j] + C[i][k] * E[k][j];
    }
  }
}
#pragma endscop
