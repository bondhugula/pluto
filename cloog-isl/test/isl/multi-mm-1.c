/* Generated from ../../../git/cloog/test/multi-mm-1.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
for (i=0;i<=N;i++) {
  for (j=0;j<=i;j++) {
    S1(i,j) ;
    S2(i,j) ;
  }
}
for (i=N+1;i<=M;i++) {
  for (j=0;j<=N;j++) {
    S1(i,j) ;
    S2(i,j) ;
  }
  for (j=N+1;j<=i;j++) {
    S1(i,j) ;
  }
}
