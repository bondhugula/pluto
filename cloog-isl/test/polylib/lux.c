/* Generated from ../../../git/cloog/test/lux.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
if (M >= 2) {
  for (l=2;l<=M;l++) {
    S1(1,1,M,l) ;
  }
  for (i=2;i<=M-1;i++) {
    for (j=1;j<=i-1;j++) {
      for (k=j+1;k<=M;k++) {
        S2(i,j,k,k,i) ;
      }
    }
    for (l=i+1;l<=M;l++) {
      S1(i,i,M,l) ;
    }
  }
  for (j=1;j<=M-1;j++) {
    for (k=j+1;k<=M;k++) {
      S2(M,j,k,k,M) ;
    }
  }
}
