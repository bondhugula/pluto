// More accesses with spatial locality with (i, j) order.
//
// CHECK: T(S1): (i, j)
// CHECK: [pluto] After intra-tile optimize
// CHECK: T(S1): (j, i)
//
// When tiling, only the intra-tile loop order is changed for locality here.
// TILE-PARALLEL: T(S1): (i, j)
// TILE-PARALLEL: [Pluto] After tiling:
// TILE-PARALLEL: T(S1): (i/32, j/32, i, j)
// TILE-PARALLEL: [pluto] After intra-tile optimize
// TILE-PARALLEL: T(S1): (i/32, j/32, j, i)
//
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < M; j++) {
    d[j][i] = d[j][i] + a[i][j];
  }
}
#pragma endscop
