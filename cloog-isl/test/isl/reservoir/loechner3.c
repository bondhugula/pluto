/* Generated from ../../../git/cloog/test/./reservoir/loechner3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.01s. */
if (M >= 1) {
  for (c2=1;c2<=M;c2++) {
    for (c4=2;c4<=c2+M;c4++) {
      for (c6=max(1,-c2+c4);c6<=min(M,c4-1);c6++) {
        S1(c2,c6,c4-c6) ;
      }
    }
  }
}
