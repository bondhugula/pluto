// A test case with a different parameter for each nest -- this tests the
// bounding function in a particular way: for the fused nest, the bounding
// function will have to be for example N + M + L.
//
// CHECK: T(S1): (i, 0)
// CHECK: T(S2): (i, 1)
// CHECK: T(S3): (i, 2)
// CHECK: Output written
#pragma scop
for (i = 0; i < N; i++) {
  a[i] = b[i];
}
for (i = 0; i < M; i++) {
  a[i] = b[i];
}
for (i = 0; i < L; i++) {
  a[i] = b[i];
}
#pragma endscop
