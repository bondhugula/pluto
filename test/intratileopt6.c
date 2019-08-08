// CHECK: T(S1): (i, j, 0, 1, 0)
// CHECK: T(S2): (i, j, 1, 0, k)
// CHECK: T(S3): (i, j, 0, 0, 0)
// CHECK: T(S4): (i, j+k, 1, 1, k)
// TILE-PARALLEL: After tiling:
// TILE-PARALLEL: T(S1): (0, i/32, j/32, i, j, 1, 0, 0, 0)
// TILE-PARALLEL: T(S2): (1, i/32, j/32, 0, k/32, i, j, 0, k)
// TILE-PARALLEL: T(S3): (0, i/32, j/32, i, j, 0, 0, 0, 0)
// TILE-PARALLEL: T(S4): (1, i/32, k/32, 1, j/32, i, k, 1, j)
// After intra-tile optimize
// TILE-PARALLEL: T(S1): (0, i/32, j/32, i, j, 1, 0, 0, 0)
// TILE-PARALLEL: T(S2): (1, i/32, j/32, 0, k/32, i, k, 0, j)
// TILE-PARALLEL: T(S3): (0, i/32, j/32, i, j, 0, 0, 0, 0)
// TILE-PARALLEL: T(S4): (1, i/32, k/32, 1, j/32, i, k, 1, j)
// CHECK: Output written

/* This 2mm kernel is taken from Polybench to test intra tile optimization.
 * With maxfuse, the intra-tile optimization has to move hyperplane j at
 * level 6 to the innermost level only for the statement S2. This test case
 * checks for intra tile optimizations for statements that (1)lie in the same
 * band (2) Are distributed in the inter tile space (3) intra-tile iterator for
 * only of the statemetent has to be moved to the innermost level. In general,
 * intra tile opt, should be able to move intra tile iterators at different
 * levels to the innermost level for statements in the band that are distributed
 * in the inter-tile space. */

#pragma scop
/* D := alpha*A*B*C + beta*D */
for (i = 0; i < _PB_NI; i++)
  for (j = 0; j < _PB_NJ; j++) {
    tmp[i][j] = 0.0;
    for (k = 0; k < _PB_NK; ++k)
      tmp[i][j] += A[i][k] * B[k][j];
  }
for (i = 0; i < _PB_NI; i++)
  for (j = 0; j < _PB_NL; j++) {
    D[i][j] *= 1.0;
    for (k = 0; k < _PB_NJ; ++k)
      D[i][j] += tmp[i][k] * C[k][j];
  }
#pragma endscop

