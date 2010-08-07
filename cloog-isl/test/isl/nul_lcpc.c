/* Generated from ../../../git/cloog/test/nul_lcpc.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.03s. */
for (i=1;i<=6;i+=2) {
  for (j=1;j<=i;j++) {
    S1(i,(i-1)/2,j) ;
    S2(i,(i-1)/2,j) ;
  }
  for (j=i+1;j<=p;j++) {
    S1(i,(i-1)/2,j) ;
  }
}
for (i=7;i<=m;i+=2) {
  for (j=1;j<=p;j++) {
    S1(i,(i-1)/2,j) ;
  }
}
