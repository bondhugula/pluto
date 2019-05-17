
#pragma scop
for (i = 0; i < size; i++) {
  for (j = 0; j < size; j++) {
    c[i][j] = a[i][j] + b[i][j];
  }
}
#pragma endscop
