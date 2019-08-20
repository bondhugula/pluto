/* Tests for loop distribution when there are no dependences between the
 * statements and one of the statements is a scalar. More precisely, must
 * distribute edges must not be added and has must distrite routine should not
 * check for such edges in this case. */

// CHECK: T(S1): (i)
// CHECK: T(S2): (0)
// [Pluto] After tiling:
// CHECK: T(S1): (i)
// CHECK: T(S2): (0)

#pragma scop
for (i = 0; i < N; i++) {
  A[i] = A[i - 1];
}
s = B[i];
#pragma endscop

