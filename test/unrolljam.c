// CHECK:   /* CLoog Unroll jammed loop */
// CHECK:for (t3=32*t1;t3<=(min(n-1,32*t1+31))-3;t3+=4) {
// CHECK:    lbv=32*t2;
// CHECK:    ubv=min(n-1,32*t2+31);
// CHECK:#pragma ivdep
// CHECK:#pragma vector always
// CHECK:    for (t4=lbv;t4<=ubv;t4++) {
// CHECK:        S1(t1,t2,t3,t4);
// CHECK:        S1(t1,t2,(t3+1),t4);
// CHECK:        S1(t1,t2,(t3+2),t4);
// CHECK:        S1(t1,t2,(t3+3),t4);
// CHECK:    }
// CHECK:}
// CHECK:for (;t3<=min(n-1,32*t1+31);t3++) {
// CHECK:    lbv=32*t2;
// CHECK:    ubv=min(n-1,32*t2+31);
// CHECK:#pragma ivdep
// CHECK:#pragma vector always
// CHECK:    for (t4=lbv;t4<=ubv;t4++) {
// CHECK:        S1(t1,t2,t3,t4);
// CHECK:    }
// CHECK:}
// CHECK: /* End of CLooG code */

/* This test case tests for unroll jamming a loop nest. It checks if the loop
 * bounds, loop strides and the epilogue of the unroll jammed loop nest are
 * updated correctly. It also checks if the epilogue is vectorized when there is
 * a vectorizable loop in a unroll jammed loop nest. */
#pragma scop
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++) {
    a[i][j] = b[i][j] + 1;
  }
}
#pragma endscop
