/* With typed fusion, these loop nests have to be distributed. The test case
 * looks for loop distribution when the two sccs are connected. */

// TYPED-FUSE-CHECK: T(S1): (0, i)
// TYPED-FUSE-CHECK: T(S2): (1, i)
// CHECK: Output written

#pragma scop
for (i = 2; i < N; i++) {
  A[i] = A[i - 1] + A[i - 2];
}

for (i = 0; i < N; i++) {
  B[i] = A[i];
}
#pragma endscop
