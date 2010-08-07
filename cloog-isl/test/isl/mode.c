/* Generated from ../../../git/cloog/test/mode.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
if (M >= 0) {
  for (i=0;i<=min(M,N-1);i++) {
    for (j=0;j<=i;j++) {
      S1(i,j) ;
      S2(i,j) ;
    }
    for (j=i+1;j<=N;j++) {
      S2(i,j) ;
    }
  }
  if ((M >= N) && (N >= 0)) {
    for (j=0;j<=N;j++) {
      S1(N,j) ;
      S2(N,j) ;
    }
  }
  if (N >= 0) {
    for (i=N+1;i<=M;i++) {
      for (j=0;j<=N;j++) {
        S1(i,j) ;
        S2(i,j) ;
      }
      for (j=N+1;j<=i;j++) {
        S1(i,j) ;
      }
    }
  }
  if (N <= -1) {
    for (i=0;i<=M;i++) {
      for (j=0;j<=i;j++) {
        S1(i,j) ;
      }
    }
  }
}
