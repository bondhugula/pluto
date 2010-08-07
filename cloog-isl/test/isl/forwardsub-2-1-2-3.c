/* Generated from ../../../git/cloog/test/forwardsub-2-1-2-3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
S3(1,0) ;
for (k=2;k<=M;k++) {
  S1(1,1,k) ;
}
for (i=2;i<=M-1;i++) {
  S4(i,0) ;
  for (k=i+1;k<=M;k++) {
    S2(i,1,k) ;
  }
}
S4(M,0) ;
