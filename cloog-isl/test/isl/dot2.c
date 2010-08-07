/* Generated from ../../../git/cloog/test/dot2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
for (i=1;i<=min(M,N);i++) {
  S1(i) ;
  for (j=1;j<=M;j++) {
    S2(i,j) ;
  }
}
for (i=N+1;i<=M;i++) {
  S1(i) ;
}
for (i=M+1;i<=N;i++) {
  for (j=1;j<=M;j++) {
    S2(i,j) ;
  }
}
