/* Pluto doesn't support negative number in transformation matrices - so can't
 * do anything on this; however, Pluto+ can. */
#pragma scop

for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
    x[i][j] = x[i][j] + y[i - 1][j];
    y[i][j] = x[i][j - 1] * y[i][j];
  }
}
#pragma endscop
