/* Even with typed fusion, these loop nests should not be distributed at
 * level 2. The test case looks for loop fusion with typed fuse, when all the
 * dimensions have been completely coloured. It does not matter whether S2 is
 * seqential or parallel. The statements S1 and S2 should not be distributed.
 * Fusion here enables wavefront parallelism.
 */

// TYPED-FUSE-CHECK: T(S1): (i, i+j, 0)
// TYPED-FUSE-CHECK: T(S2): (i, 2i+1, 1)
// [Pluto] After tiling:
// TYPED-FUSE-CHECK: T(S1): (i/32, (i+j)/32, i, i+j, 0)
// TYPED-FUSE-CHECK: T(S2): (i/32, zT2, i, 2i+1, 1)
// [Pluto] After tile scheduling:
// TYPED-FUSE-CHECK: T(S1): (i/32+(i+j)/32, (i+j)/32, i, i+j, 0)
// TYPED-FUSE-CHECK: T(S2): (i/32+zT2, zT2, i, 2i+1, 1)
// CHECK: Output written

#pragma scop
for (i = 1; i < N; i++) {
  for (j = 1; j < N; j++) {
    A[i][j] = A[i - 1][j] + A[i][j - 1]; // S1
  }
  A[i][i] = SQRT_FUN(A[i][i - 1]); // S2
}
#pragma endscop
