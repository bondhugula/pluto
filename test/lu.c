
#pragma scop
for (k = 0; k < N; k++) {
  for (j = k + 1; j < N; j++) {
    a[k][j] = a[k][j] / a[k][k];
  }
  for (i = k + 1; i < N; i++) {
    for (j = k + 1; j < N; j++) {
      a[i][j] = a[i][j] - a[i][k] * a[k][j];
    }
  }
}
#pragma endscop
