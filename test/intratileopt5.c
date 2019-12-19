// CHECK: T(S1): (0, i, j, 1, 0)
// CHECK: T(S2): (1, i, j, 0, k)
// CHECK: T(S3): (0, i, j, 0, 0)
// CHECK: T(S4): (1, i, k, 1, j)
// TILE-PARALLEL: After tiling:
// TILE-PARALLEL: T(S1): (0, i/32, j/32, i, j, 1, 0, 0, 0)
// TILE-PARALLEL: T(S2): (1, i/32, j/32, 0, k/32, i, j, 0, k)
// TILE-PARALLEL: T(S3): (0, i/32, j/32, i, j, 0, 0, 0, 0)
// TILE-PARALLEL: T(S4): (1, i/32, k/32, 1, j/32, i, k, 1, j)
// After intra-tile optimize
// TILE-PARALLEL: T(S1): (0, i/32, j/32, i, j, 1, 0, 0, 0)
// TILE-PARALLEL: T(S2): (1, i/32, j/32, 0, k/32, i, 0, k, j)
// TILE-PARALLEL: T(S3): (0, i/32, j/32, i, j, 0, 0, 0, 0)
// TILE-PARALLEL: T(S4): (1, i/32, k/32, 1, j/32, i, k, 1, j)
// CHECK: Output written

/* This 2mm kernel is taken from Polybench to test intra tile optimization.
 * Intra-tile optimization has to move hyperplane j at level 6 to the innermost
 * level only for the statement S2. This test case checks for intra tile
 * optimizations for statements that lie in the same band and require different
 * intra-tile permutations for statements that do not share intra-tile loops.
 * When loop nest is not tiled, loops in the innermost permutable band should
 * not be moved acorss scalar hyperplanes, even if outer and inner permutable
 * bands are the same. */

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
