/* Generated from ../../../git/cloog/test/./reservoir/jacobi3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
if ((M >= 1) && (N >= 3)) {
  for (c4=2;c4<=N-1;c4++) {
    for (c6=2;c6<=N-1;c6++) {
      S1(1,c4,c6) ;
    }
  }
  for (c2=3;c2<=2*M;c2++) {
    for (c4=2;c4<=N-1;c4++) {
      for (c6=2;c6<=N-1;c6++) {
        if (c2%2 == 0) {
          S1(c2/2,c4,c6) ;
        }
      }
    }
    for (c4=2;c4<=N-1;c4++) {
      for (c6=2;c6<=N-1;c6++) {
        if ((c2+1)%2 == 0) {
          S2((c2-1)/2,c4,c6) ;
        }
      }
    }
  }
  for (c4=2;c4<=N-1;c4++) {
    for (c6=2;c6<=N-1;c6++) {
      S2(M,c4,c6) ;
    }
  }
}
