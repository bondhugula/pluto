/* Generated from ../../../git/cloog/test/min-1-1.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if ((M >= 0) && (N >= 1)) {
  for (i=1;i<=N;i++) {
    for (j=0;j<=min(min(M,i),-i+N);j++) {
      S1(i,j) ;
    }
  }
}
