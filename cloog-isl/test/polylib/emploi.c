/* Generated from ../../../git/cloog/test/emploi.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if (n >= 1) {
  if (m >= 1) {
    for (i=1;i<=n;i++) {
      if (i >= m) {
        S1(i) ;
      }
      if (i <= min(2*m,m-1)) {
        S1(i) ;
      }
      for (j=1;j<=m;j++) {
        S2(i,j) ;
      }
    }
  }
  if (m <= 0) {
    for (i=1;i<=n;i++) {
      S1(i) ;
    }
  }
}
