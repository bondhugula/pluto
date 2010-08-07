/* Generated from ../../../git/cloog/test/reservoir/bastoul3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.00s. */
for (i=3;i<=9;i++) {
  for (j=max(1,i-6);j<=min(3,i-2);j++) {
    if ((i+j)%2 == 0) {
      S1(i,j,(i-j)/2) ;
    }
  }
}
