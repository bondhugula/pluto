
#pragma scop
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    a[i][j] = a[i][j] + c[i][j];
  }
}
#pragma endscop
