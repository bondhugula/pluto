/* Generated from ../../../git/cloog/test/infinite2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.00s. */
for (i=1;i<=N;i++) {
  S1(i) ;
  for (j=1;j<=M;j++) {
    S2(i,j) ;
  }
}
for (i=N+1;;i++) {
  S1(i) ;
}
