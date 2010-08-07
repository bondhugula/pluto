/* Generated from ../../../git/cloog/test/./non_optimal/youcef.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
for (i=0;i<=3;i++) {
  S1(i,i) ;
  S2(i,i) ;
  for (j=i+1;j<=4;j++) {
    S2(i,j) ;
  }
  S2(i,5) ;
  S3(i,5) ;
}
S1(4,4) ;
S2(4,4) ;
S2(4,5) ;
S3(4,5) ;
S1(5,5) ;
S2(5,5) ;
S3(5,5) ;
