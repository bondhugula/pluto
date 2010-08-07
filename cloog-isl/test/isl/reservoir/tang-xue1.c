/* Generated from ../../../git/cloog/test/reservoir/tang-xue1.cloog by CLooG 0.14.0-165-g01eb246 gmp bits in 0.01s. */
for (c2=0;c2<=9;c2++) {
  for (c4=max(-1,c2-9);c4<=min(4,c2+3);c4++) {
    for (c6=max(max(1,c2),c2-c4);c6<=min(min(9,c2+1),c2-c4+4);c6++) {
      for (c8=max(1,-c2+c4+c6);c8<=min(4,-c2+c4+c6+1);c8++) {
        if (c2%2 == 0) {
          if (c4%2 == 0) {
            S1(c2/2,(-c2+c4)/2,-c2+c6,-c4+c8);
          }
        }
      }
    }
  }
}
