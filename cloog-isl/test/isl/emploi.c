/* Generated from ../../../git/cloog/test/emploi.cloog by CLooG 0.14.0-162-g1e599e0 gmp bits in 0.00s. */
if (n >= 1) {
  if (m >= 1) {
    for (i=1;i<=n;i++) {
      if (i <= 2*m) {
        S1(i);
      }
      if (i >= max(m,2*m+1)) {
        S1(i);
      }
      for (j=1;j<=m;j++) {
        S2(i,j);
      }
    }
  }
  if (m <= 0) {
    for (i=1;i<=n;i++) {
      S1(i);
    }
  }
}
