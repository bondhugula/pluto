constant N, z;

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    a[i][j] = a[i + z][j] + a[i][j - 1];
  }
}
