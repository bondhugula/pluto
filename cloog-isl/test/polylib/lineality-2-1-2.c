/* Generated from ../../../git/cloog/test/lineality-2-1-2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
for (i=1;i<=M-3;i++) {
  for (j=1;j<=i+1;j++) {
    S1(i,j) ;
  }
  S1(i,i+2) ;
  S2(i,i+2) ;
  for (j=i+3;j<=M;j++) {
    S1(i,j) ;
  }
}
if (M >= 3) {
  for (j=1;j<=M-1;j++) {
    S1(M-2,j) ;
  }
  S1(M-2,M) ;
  S2(M-2,M) ;
}
for (i=M-1;i<=M;i++) {
  for (j=1;j<=M;j++) {
    S1(i,j) ;
  }
  S2(i,i+2) ;
}
