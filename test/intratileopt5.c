// CHECK: T(S1): (i, j, 0, 1, 0)
// CHECK: T(S2): (i, j, 1, 0, k)
// CHECK: T(S3): (i, j, 0, 0, 0)
// CHECK: T(S4): (i, j+k, 1, 1, k)
// After tiling:
// CHECK: T(S1): (i/32, j/32, 0, 1, 0, i, j, 0, 1, 0)
// CHECK: T(S2): (i/32, j/32, 1, 0, k/32, i, j, 1, 0, k)
// CHECK: T(S3): (i/32, j/32, 0, 0, 0, i, j, 0, 0, 0)
// CHECK: T(S4): (i/32, (j+k)/32, 1, 1, k/32, i, j+k, 1, 1, k)
// After intra-tile optimize
// CHECK: T(S1): (i/32, j/32, 0, 1, 0, i, j, 0, 1, 0)
// CHECK: T(S2): (i/32, j/32, 1, 0, k/32, i, k, 1, 0, j)
// CHECK: T(S3): (i/32, j/32, 0, 0, 0, i, j, 0, 0, 0)
// CHECK: T(S4): (i/32, (j+k)/32, 1, 1, k/32, i, k, 1, 1, j+k)
// CHECK : Output Written

/* This is a 2mm kernel is taken from Polybench to test intra tile optimization.
 * With maxfuse, the intra-tile optimization has to move hyperplane (j+k) at
 * level 6 to the innermost level for the statement S4. This test case thus
 * tests for intra tile optimization to move loops across scalar hyperplanes. */
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

