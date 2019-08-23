// CHECK: for (t3=32*t1;t3<=min(n-1,32*t1+31);t3++) {
// CHECK:     lbv=max(32*t2,t3);
// CHECK:     ubv=min(n-1,32*t2+31);
// CHECK: #pragma ivdep
// CHECK: #pragma vector always
// CHECK:     for (t4=lbv;t4<=ubv;t4++) {
// CHECK:         S1(t1,t2,t3,t4);
// CHECK:     }
// CHECK: }
// CHECK: /* End of CLooG code */

/* The loop i can not be unroll jammed as the domain is not rectangular. */
#pragma scop
for (i = 0; i < n; i++) {
  for (j = i; j < n; j++) {
    a[i][j] = b[i][j] + 1;
  }
}
#pragma endscop
