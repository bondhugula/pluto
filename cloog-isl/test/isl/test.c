/* Generated from ../../../git/cloog/test/test.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.12s. */
for (i=1;i<=2;i++) {
  for (j=1;j<=M;j++) {
    S1(i,j) ;
  }
}
for (i=3;i<=M-1;i++) {
  for (j=1;j<=i-1;j++) {
    S1(i,j) ;
  }
  S1(i,i) ;
  S2(i,i) ;
  for (j=i+1;j<=M;j++) {
    S1(i,j) ;
  }
}
for (j=1;j<=M-1;j++) {
  S1(M,j) ;
}
S1(M,M) ;
S2(M,M) ;
for (i=M+1;i<=N;i++) {
  for (j=1;j<=M;j++) {
    S1(i,j) ;
  }
  S2(i,i) ;
}
