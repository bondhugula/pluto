
#pragma scop
for (i = 1; i <= N; i++) {
  for (j = 1; j <= N; j++) {
    s = s + a[i][j];
  }
}
#pragma endscop
