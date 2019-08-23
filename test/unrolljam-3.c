// CHECK: for (t1=0;t1<=n-1;t1++) {
// CHECK:   for (t2=0;t2<=n-1;t2++) {
// CHECK:     S1(t1,t2);
// CHECK:   }
// CHECK: }
// CHECK: /* End of CLooG code */

/* The loop nest can not be unroll jammed as the domain is not permutable. */
#pragma scop
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    s = s + b[i][j] + 1;
  }
}
#pragma endscop
