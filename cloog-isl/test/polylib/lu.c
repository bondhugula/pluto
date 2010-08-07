/* Generated from ../../../git/cloog/test/lu.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if (n >= 2) {
  for (j=2;j<=n;j++) {
    S1(1,j) ;
  }
  for (c1=2;c1<=n-1;c1++) {
    for (c2=2;c2<=n-1;c2++) {
      for (i=1;i<=min(c1-1,c2-1);i++) {
        S2(i,c2,c1) ;
      }
    }
    for (i=1;i<=c1-1;i++) {
      S2(i,n,c1) ;
    }
    for (j=c1+1;j<=n;j++) {
      S1(c1,j) ;
    }
  }
  for (c2=2;c2<=n;c2++) {
    for (i=1;i<=c2-1;i++) {
      S2(i,c2,n) ;
    }
  }
}
