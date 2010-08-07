/* Generated from ../../../git/cloog/test/esced.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if (m >= 1) {
  if (n >= 1) {
    for (i=1;i<=m;i++) {
      S1(i) ;
      for (j=1;j<=n;j++) {
        S2(i,j) ;
      }
    }
  }
  if (n <= 0) {
    for (i=1;i<=m;i++) {
      S1(i) ;
    }
  }
}
