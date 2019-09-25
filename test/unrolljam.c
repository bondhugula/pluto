// CHECK:for (t3=32*t1;t3<=min(n-1,32*t1+31);t3++) {
// CHECK:    lbv=32*t2;
// CHECK:    ubv=min(n-1,32*t2+31);
// CHECK:#pragma ivdep
// CHECK:#pragma vector always
// CHECK:    for (t4=lbv;t4<=ubv;t4++) {
// CHECK:        S1(t1,t2,t3,t4);
// CHECK:    }
// CHECK:}
// CHECK: /* End of CLooG code */

/* This test case tests for profitability of unroll jamming a loop nest. The
 * loop nest does not have any dependences and hence unroll jamming is not
 * considered to be profitable. */
#pragma scop
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    a[i][j] = b[i][j] + 1;
  }
}
#pragma endscop
