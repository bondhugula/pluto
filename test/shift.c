/* From Darte's paper on medium-grained parallelism detection.
 * S1 will be shifted forward by 1 (relative to S2) and appear after
 * S2 in the transformed body. */
// CHECK: T(S1): (i, j, 1)
// CHECK: T(S2): (i, j+1, 0)
// CHECK: Output written
#pragma scop
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    a[i][j] = b[i][j - 1] + a[i][j - 1];
    b[i][j] = a[i - 1][j + 1] + b[i - 1][j];
  }
}
#pragma endscop
