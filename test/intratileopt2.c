// CHECK: T(S1): (i, j)
// CHECK: [pluto] After intra-tile optimize
// CHECK: T(S1): (j, i)
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < M; j++) {
    d[j][i] = d[j][i] + a[i][j];
  }
}
#pragma endscop
