#pragma scop

for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    a[i][j] = a[i - 2][j] + a[i - 1][j + 1];
  }
}
/* pluto end */
#pragma endscop
