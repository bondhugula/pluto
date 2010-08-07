/* Generated from ../../../git/cloog/test/lineality-1-2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
S1(1,1) ;
S2(1,1) ;
for (j=2;j<=M;j++) {
  S1(1,j) ;
}
for (i=2;i<=M-1;i++) {
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
