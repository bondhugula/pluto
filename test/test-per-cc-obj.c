// CC-OBJ-CHECK: T(S1): (0, i, k, 2)
// CC-OBJ-CHECK: T(S2): (0, i, k, 3)
// CC-OBJ-CHECK: T(S3): (1, i+1, 0, 0)
// CC-OBJ-CHECK: T(S4): (1, i, 0, 1)
// CHECK: After intra-tile optimize
// CC-OBJ-CHECK: T(S1): (0, k, i, 2)
// CC-OBJ-CHECK: T(S2): (0, k, i, 3)
// CC-OBJ-CHECK: T(S3): (1, i+1, 0, 0)
// CC-OBJ-CHECK: T(S4): (1, i, 0, 1)
// CHECK: Output written

/* Had per CC objective not been used, then Pluto would have found the solution
 * T(S3): (1, i, 0, 0), T(S4): (1, i, 0, 1), which is not parallel for
 * statements S3 and S4 */
#pragma scop
for (k = 0; k < N; k++) {
  for (i = 0; i < N; i++) {
    A[k][i] = i;
  }
  for (i = 0; i < N; i++) {
    B[k][i] = A[k][i - 1] + A[N - k][i];
  }
}
for (i = 0; i < N; i++) {
  C[i] = B[0][i];
}
for (i = 0; i < N; i++) {
  D[i] = C[i - 1];
}
#pragma endscop
