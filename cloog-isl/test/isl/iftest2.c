/* Generated from ../../../git/cloog/test/iftest2.cloog by CLooG 0.14.0-162-g1e599e0 gmp bits in 0.00s. */
if ((M >= 1) && (N >= 1)) {
  for (i=1;i<=N;i++) {
    for (j=1;j<=M;j++) {
      if (i <= 2*M) {
        S1(i,j);
      }
      if (i >= max(M,2*M+1)) {
        S1(i,j);
      }
    }
  }
}
