// CHECK: for (t4=32*t1;t4<=(min(M-1,32*t1+31))-1;t4+=2) {
// CHECK:   for (t5=32*t3;t5<=(min(K-1,32*t3+31))-1;t5+=2) {
// CHECK:     lbv=32*t2;
// CHECK:     ubv=min(N-1,32*t2+31);
// CHECK:     #pragma ivdep
// CHECK:     #pragma vector always
// CHECK:     for (t6=lbv;t6<=ubv;t6++) {
// CHECK:       S1(t1,t2,t3,t4,t6,t5);
// CHECK:       S1(t1,t2,t3,(t4+1),t6,t5);
// CHECK:       S1(t1,t2,t3,t4,t6,(t5+1));
// CHECK:       S1(t1,t2,t3,(t4+1),t6,(t5+1));
// CHECK:     }
// CHECK:   }
// CHECK:   for (;t5<=min(K-1,32*t3+31);t5++) {
// CHECK:     lbv=32*t2;
// CHECK:     ubv=min(N-1,32*t2+31);
// CHECK:     #pragma ivdep
// CHECK:     #pragma vector always
// CHECK:     for (t6=lbv;t6<=ubv;t6++) {
// CHECK:       S1(t1,t2,t3,t4,t6,t5);
// CHECK:       S1(t1,t2,t3,(t4+1),t6,t5);
// CHECK:     }
// CHECK:   }
// CHECK: }
// CHECK: for (;t4<=min(M-1,32*t1+31);t4++) {
// CHECK:   for (t5=32*t3;t5<=min(K-1,32*t3+31);t5++) {
// CHECK:     lbv=32*t2;
// CHECK:     ubv=min(N-1,32*t2+31);
// CHECK:     #pragma ivdep
// CHECK:     #pragma vector always
// CHECK:     for (t6=lbv;t6<=ubv;t6++) {
// CHECK:       S1(t1,t2,t3,t4,t6,t5);
// CHECK:     }
// CHECK:   }
// CHECK: }
/* Checks for multi-loop unroll jam after intra tile optimize. The i and k loops
 * are unroll and jammed by a factor of 2. */
#pragma scop
for (i = 0; i < M; i++)
  for (j = 0; j < N; j++)
    for (k = 0; k < K; k++)
      C[i][j] = C[i][j] + A[i][k] * B[k][j];
#pragma endscop
