/* Generated from ../../../git/cloog/test/gesced2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.05s. */
for (c1=1;c1<=4;c1++) {
  for (c2=5;c2<=M-10;c2++) {
    S1(c1,c2) ;
  }
}
for (c1=5;c1<=min(9,M-10);c1++) {
  for (c2=-c1+1;c2<=4;c2++) {
    S2(c1+c2,c1) ;
  }
  for (c2=5;c2<=M-10;c2++) {
    S1(c1,c2) ;
    S2(c1+c2,c1) ;
  }
  for (c2=M-9;c2<=-c1+M;c2++) {
    S2(c1+c2,c1) ;
  }
}
if (M >= 20) {
  for (c2=-9;c2<=4;c2++) {
    S2(c2+10,10) ;
  }
  for (c2=5;c2<=M-10;c2++) {
    S1(10,c2) ;
    S2(c2+10,10) ;
  }
}
for (c1=11;c1<=M-10;c1++) {
  for (c2=-c1+1;c2<=4;c2++) {
    S2(c1+c2,c1) ;
  }
  for (c2=5;c2<=-c1+M;c2++) {
    S1(c1,c2) ;
    S2(c1+c2,c1) ;
  }
  for (c2=-c1+M+1;c2<=M-10;c2++) {
    S1(c1,c2) ;
  }
}
for (c1=M-9;c1<=M;c1++) {
  for (c2=5;c2<=M-10;c2++) {
    S1(c1,c2) ;
  }
}
