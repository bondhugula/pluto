// Outer parallelism + intra-tile locality: while the (i, j) loop order is
// better for coarse-grained parallelism, the (j, i) is better for locality.
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
    a[i] += b[j]*c[j][i];
  }
}
#pragma endscop
