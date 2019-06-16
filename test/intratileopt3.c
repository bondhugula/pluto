//
// Per Pluto intra-tile loop order cost function, accesses with temporal reuse
// are more important than spatial reuse.
//
// CHECK: T(S1): (j, i)
//
// For tiling and parallelization as well as for intra-tile locality, it's
// profitable to have the j loop outside.
// TILE-PARALLEL: T(S1): (j, i)
// TILE-PARALLEL: [Pluto] After tiling:
// TILE-PARALLEL: T(S1): (j/32, i/32, j, i)
// TILE-PARALLEL: [pluto] After intra-tile optimize
// TILE-PARALLEL: T(S1): (j/32, i/32, j, i)
//
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < M; j++) {
    d[j] += a[i][j];
  }
}
#pragma endscop
