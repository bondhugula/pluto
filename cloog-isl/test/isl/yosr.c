/* Generated from ../../../git/cloog/test/yosr.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if (n >= 2) {
  for (j=2;j<=n;j++) {
    S1(1,j) ;
  }
  for (proc=2;proc<=n-1;proc++) {
    for (i=1;i<=proc-1;i++) {
      for (j=i+1;j<=n;j++) {
        S2(i,j,proc) ;
      }
    }
    for (j=proc+1;j<=n;j++) {
      S1(proc,j) ;
    }
  }
  for (i=1;i<=n-1;i++) {
    for (j=i+1;j<=n;j++) {
      S2(i,j,n) ;
    }
  }
}
