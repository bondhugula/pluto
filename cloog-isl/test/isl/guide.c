/* Generated from ../../../git/cloog/test/guide.cloog by CLooG 0.14.0-162-g1e599e0 gmp bits in 0.00s. */
if (N >= 1) {
  for (i=1;i<=N;i++) {
    if (i <= 2*M) {
      S1(i);
    }
    if (i >= max(M,2*M+1)) {
      S1(i);
    }
  }
  for (i=N+1;i<=2*N;i++) {
    S2(i);
  }
}
