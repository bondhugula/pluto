/* Generated from ../../../git/cloog/test/./non_optimal/nul_complex1.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.00s. */
if (n >= 0) {
  for (c1=0;c1<=5*n;c1++) {
    for (c2=max(ceild(2*c1,3),c1-n);c2<=min(floord(2*c1+2*n,3),c1);c2++) {
      if (c2%2 == 0) {
        S1((-2*c1+3*c2)/2,c1-c2) ;
      }
    }
  }
}
