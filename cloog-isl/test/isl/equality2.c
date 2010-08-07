/* Generated from ../../../git/cloog/test/equality2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.07s. */
for (i0=1;i0<=10000;i0++) {
  for (i1=1000;i1<=1016;i1++) {
    for (i2=1;i2<=min(-2*i1+2033,2*i1-1999);i2++) {
      if (2*i1 == i2+1999) {
        S2(i0,i1,i2,1,i0,2*i1-1000,1,2,i0,i1-499,2*i1-1999,i0,2*i1-1999,i1-999,i1-999) ;
      }
      if (i2 == 1) {
        if (i1%2 == 0) {
          S1(i0,i1,i2,2,i0,(i1+2)/2,i1-999,i0,i1-999,(i1-998)/2,(i1-998)/2) ;
        }
      }
    }
  }
}
