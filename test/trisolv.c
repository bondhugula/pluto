#pragma scop
for (i = 0; i < N; i++) {
  x[i] = c[i];
  for (j = 0; j <= i - 1; j++) {
    x[i] = x[i] - a[i][j] * x[j];
  }
  x[i] = x[i] / a[i][i];
}

#pragma endscop
