/* Generated from ../../../git/cloog/test/byu98-1-2-3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
for (i=2;i<=3;i++) {
  for (j=-i+6;j<=6;j++) {
    S1(i,j) ;
  }
}
for (j=3;j<=4;j++) {
  S1(4,j) ;
}
S1(4,5) ;
S2(4,5) ;
S1(4,6) ;
S1(5,4) ;
S2(5,4) ;
for (j=5;j<=6;j++) {
  S1(5,j) ;
}
for (i=6;i<=7;i++) {
  S2(i,-i+9) ;
  for (j=i-1;j<=6;j++) {
    S1(i,j) ;
  }
}
S2(8,1) ;
