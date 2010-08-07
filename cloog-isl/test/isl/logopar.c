/* Generated from ../../../git/cloog/test/logopar.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.04s. */
for (j=0;j<=m;j++) {
  S1(1,j) ;
}
if (m >= n+1) {
  for (i=2;i<=n;i++) {
    for (j=0;j<=i-2;j++) {
      S2(i,j) ;
    }
    for (j=i-1;j<=n;j++) {
      S1(i,j) ;
      S2(i,j) ;
    }
    for (j=n+1;j<=m;j++) {
      S1(i,j) ;
    }
  }
}
if (m == n) {
  for (i=2;i<=m;i++) {
    for (j=0;j<=i-2;j++) {
      S2(i,j) ;
    }
    for (j=i-1;j<=m;j++) {
      S1(i,j) ;
      S2(i,j) ;
    }
  }
}
for (i=n+1;i<=m+1;i++) {
  for (j=i-1;j<=m;j++) {
    S1(i,j) ;
  }
}
