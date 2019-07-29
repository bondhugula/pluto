// CHECK: T(S1): (i, j, 0, 1, 0)
// CHECK: T(S2): (i, j, 1, 0, k)
// CHECK: T(S3): (i, j, 0, 0, 0)
// CHECK: T(S4): (i, j+k, 1, 1, k)
// CHECK-TILE: After tiling:
// CHECK-TILE: T(S1): (i/32, j/32, 0, 1, 0, i, j, 0, 1, 0)
// CHECK-TILE: T(S2): (i/32, j/32, 1, 0, k/32, i, j, 1, 0, k)
// CHECK-TILE: T(S3): (i/32, j/32, 0, 0, 0, i, j, 0, 0, 0)
// CHECK-TILE: T(S4): (i/32, (j+k)/32, 1, 1, k/32, i, j+k, 1, 1, k)
// After intra-tile optimize
// CHECK-TILE: T(S1): (i/32, j/32, 0, 1, 0, i, j, 0, 1, 0)
// CHECK-TILE: T(S2): (i/32, j/32, 1, 0, k/32, i, k, 1, 0, j)
// CHECK-TILE: T(S3): (i/32, j/32, 0, 0, 0, i, j, 0, 0, 0)
// CHECK-TILE: T(S4): (i/32, (j+k)/32, 1, 1, k/32, i, k, 1, 1, j+k)
// CHECK-NOTILE: T(S1): (i, j, 0, 1, 0)
// CHECK-NOTILE: T(S2): (i, j, 1, 0, k)
// CHECK-NOTILE: T(S3): (i, j, 0, 0, 0)
// CHECK-NOTILE: T(S4): (i, j+k, 1, 1, k)
// CHECK: Output written

/* This is a 2mm kernel is taken from Polybench to test intra tile optimization.
 * With maxfuse, the intra-tile optimization has to move hyperplane (j+k) at
 * level 6 to the innermost level for the statement S4. This test case thus
 * tests for intra tile optimization to move loops across scalar hyperplanes. If
 * the loop nest is not tiled, then it will be incorrect to move a loop across
 * scalar hyperplane in the loop nest. The key here is that in case of a tiled
 * loop nest, the distribution would be made in the tile space (in the outermost
 * permutable band) and not within a tile. */
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

