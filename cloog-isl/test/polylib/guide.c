/* Generated from ../../../git/cloog/test/guide.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.00s. */
if (N >= 1) {
  for (i=1;i<=N;i++) {
    if (i >= M) {
      S1(i) ;
    }
    if (i <= min(2*M,M-1)) {
      S1(i) ;
    }
  }
  for (i=N+1;i<=2*N;i++) {
    S2(i) ;
  }
}
