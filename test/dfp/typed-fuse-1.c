// TYPED-FUSE-CHECK: T(S1): (1, i)
// TYPED-FUSE-CHECK: T(S2): (0, i)
// CHECK: Output written

/* With typed fusion, these loop nests have to be distributed. The test case
 * looks for loop distribution when the two sccs are not connected. */
#pragma scop
for (i = 0; i < N; i++) {
  s += A[i];
}
for (i = 0; i < N; i++) {
  B[i] = i;
}
#pragma endscop
