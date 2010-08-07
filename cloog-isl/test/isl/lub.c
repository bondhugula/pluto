/* Generated from ../../../git/cloog/test/lub.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
if (M >= 2) {
  for (i=1;i<=M-1;i++) {
    for (j=i+1;j<=M;j++) {
      S1(i,j) ;
      for (k=i+1;k<=M;k++) {
        S2(i,j,k) ;
        S3(i,j,k) ;
      }
      S4(i,j) ;
    }
  }
}
