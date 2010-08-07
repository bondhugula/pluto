/* Generated from ../../../git/cloog/test/./reservoir/fusion2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
if ((M >= 1) && (N >= 1)) {
  for (c4=1;c4<=M;c4++) {
    S1(1,c4) ;
  }
  for (c2=2;c2<=N;c2++) {
    for (c4=1;c4<=M;c4++) {
      S2(c2-1,c4) ;
    }
    for (c4=1;c4<=M;c4++) {
      S1(c2,c4) ;
    }
  }
  for (c4=1;c4<=M;c4++) {
    S2(N,c4) ;
  }
}
