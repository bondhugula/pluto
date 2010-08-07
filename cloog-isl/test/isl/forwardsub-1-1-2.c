/* Generated from ../../../git/cloog/test/forwardsub-1-1-2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
S3(1,1) ;
S1(2,1) ;
S4(2,2) ;
for (i=3;i<=M;i++) {
  S1(i,1) ;
  for (j=2;j<=i-1;j++) {
    S2(i,j) ;
  }
  S4(i,i) ;
}
