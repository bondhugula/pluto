/* Generated from ../../../git/cloog/test/forwardsub-3-1-2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
S3(2,1) ;
S1(3,1) ;
S1(4,1) ;
S4(4,2) ;
for (i=5;i<=M+1;i++) {
  S1(i,1) ;
  for (j=2;j<=floord(i-1,2);j++) {
    S2(i,j) ;
  }
  if (i%2 == 0) {
    S4(i,i/2) ;
  }
}
for (i=M+2;i<=2*M-1;i++) {
  for (j=i-M;j<=floord(i-1,2);j++) {
    S2(i,j) ;
  }
  if (i%2 == 0) {
    S4(i,i/2) ;
  }
}
S4(2*M,M) ;
