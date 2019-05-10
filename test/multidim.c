
#pragma scop

/* pluto start (N) */

for (i = 1; i <= N; i++) {
  for (j = 1; j <= N; j++) {
    s = s + a[i][j];
  }
}

/* pluto end */
#pragma endscop
