/* Generated from ../../../git/cloog/test/classen.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 1.82s. */
if (m >= 1) {
  if (m >= 2) {
    S1(0,1,1,1) ;
    S2(0,1,1,1,1,1,2,1) ;
    S3(0,1,1,2,1,1,1,2) ;
    S4(0,1,2,2,1,1,2,2) ;
    S8(0,1) ;
  }
  if (m == 1) {
    S1(0,1,1,1) ;
    S8(0,1) ;
  }
  if (m >= 3) {
    S5(0,1,1,1,1,1,2,1) ;
    S1(1,1,2,1) ;
    S2(1,1,2,1,2,1,3,1) ;
    S3(1,1,2,2,2,1,2,2) ;
    S4(1,1,3,2,2,1,3,2) ;
    S6(0,1,1,2,1,1,1,2) ;
    S7(0,1,2,2,1,1,2,2) ;
    S1(1,2,1,2) ;
    S2(1,2,2,2,1,2,2,2) ;
    S3(1,2,2,3,1,2,1,3) ;
    S4(1,2,3,3,1,2,2,3) ;
    for (coordP1=1;coordP1<=2;coordP1++) {
      S8(1,coordP1) ;
    }
  }
  for (glT1=2;glT1<=m-2;glT1++) {
    S5(glT1-1,1,glT1,1,glT1,1,glT1+1,1) ;
    S1(glT1,1,glT1+1,1) ;
    S2(glT1,1,glT1+1,1,glT1+1,1,glT1+2,1) ;
    S3(glT1,1,glT1+1,2,glT1+1,1,glT1+1,2) ;
    S4(glT1,1,glT1+2,2,glT1+1,1,glT1+2,2) ;
    for (rp1=2;rp1<=glT1;rp1++) {
      S5(glT1-1,rp1,glT1,rp1,glT1-rp1+1,rp1,glT1-rp1+2,rp1) ;
      S6(glT1-1,rp1-1,glT1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+2,rp1) ;
      S7(glT1-1,rp1-1,glT1+1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+3,rp1) ;
      S1(glT1,rp1,glT1-rp1+2,rp1) ;
      S2(glT1,rp1,glT1+1,rp1,glT1-rp1+2,rp1,glT1-rp1+3,rp1) ;
      S3(glT1,rp1,glT1+1,rp1+1,glT1-rp1+2,rp1,glT1-rp1+2,rp1+1) ;
      S4(glT1,rp1,glT1+2,rp1+1,glT1-rp1+2,rp1,glT1-rp1+3,rp1+1) ;
    }
    S6(glT1-1,glT1,glT1,glT1+1,1,glT1,1,glT1+1) ;
    S7(glT1-1,glT1,glT1+1,glT1+1,1,glT1,2,glT1+1) ;
    S1(glT1,glT1+1,1,glT1+1) ;
    S2(glT1,glT1+1,glT1+1,glT1+1,1,glT1+1,2,glT1+1) ;
    S3(glT1,glT1+1,glT1+1,glT1+2,1,glT1+1,1,glT1+2) ;
    S4(glT1,glT1+1,glT1+2,glT1+2,1,glT1+1,2,glT1+2) ;
    for (coordP1=1;coordP1<=glT1+1;coordP1++) {
      S8(glT1,coordP1) ;
    }
  }
  if (m >= 3) {
    S5(m-2,1,m-1,1,m-1,1,m,1) ;
    S1(m-1,1,m,1) ;
    S3(m-1,1,m,2,m,1,m,2) ;
    for (rp1=2;rp1<=m-1;rp1++) {
      S5(m-2,rp1,m-1,rp1,-rp1+m,rp1,-rp1+m+1,rp1) ;
      S6(m-2,rp1-1,m-1,rp1,-rp1+m+1,rp1-1,-rp1+m+1,rp1) ;
      S7(m-2,rp1-1,m,rp1,-rp1+m+1,rp1-1,-rp1+m+2,rp1) ;
      S1(m-1,rp1,-rp1+m+1,rp1) ;
      S2(m-1,rp1,m,rp1,-rp1+m+1,rp1,-rp1+m+2,rp1) ;
      S3(m-1,rp1,m,rp1+1,-rp1+m+1,rp1,-rp1+m+1,rp1+1) ;
      S4(m-1,rp1,m+1,rp1+1,-rp1+m+1,rp1,-rp1+m+2,rp1+1) ;
    }
    S6(m-2,m-1,m-1,m,1,m-1,1,m) ;
    S7(m-2,m-1,m,m,1,m-1,2,m) ;
    S1(m-1,m,1,m) ;
    S2(m-1,m,m,m,1,m,2,m) ;
    for (coordP1=1;coordP1<=m;coordP1++) {
      S8(m-1,coordP1) ;
    }
  }
  for (glT1=m;glT1<=2*m-4;glT1++) {
    S5(glT1-1,glT1-m+2,glT1,glT1-m+2,m-1,glT1-m+2,m,glT1-m+2) ;
    S6(glT1-1,glT1-m+1,glT1,glT1-m+2,m,glT1-m+1,m,glT1-m+2) ;
    S1(glT1,glT1-m+2,m,glT1-m+2) ;
    S3(glT1,glT1-m+2,glT1+1,glT1-m+3,m,glT1-m+2,m,glT1-m+3) ;
    for (rp1=glT1-m+3;rp1<=m-1;rp1++) {
      S5(glT1-1,rp1,glT1,rp1,glT1-rp1+1,rp1,glT1-rp1+2,rp1) ;
      S6(glT1-1,rp1-1,glT1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+2,rp1) ;
      S7(glT1-1,rp1-1,glT1+1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+3,rp1) ;
      S1(glT1,rp1,glT1-rp1+2,rp1) ;
      S2(glT1,rp1,glT1+1,rp1,glT1-rp1+2,rp1,glT1-rp1+3,rp1) ;
      S3(glT1,rp1,glT1+1,rp1+1,glT1-rp1+2,rp1,glT1-rp1+2,rp1+1) ;
      S4(glT1,rp1,glT1+2,rp1+1,glT1-rp1+2,rp1,glT1-rp1+3,rp1+1) ;
    }
    S5(glT1-1,m,glT1,m,glT1-m+1,m,glT1-m+2,m) ;
    S6(glT1-1,m-1,glT1,m,glT1-m+2,m-1,glT1-m+2,m) ;
    S7(glT1-1,m-1,glT1+1,m,glT1-m+2,m-1,glT1-m+3,m) ;
    S1(glT1,m,glT1-m+2,m) ;
    S2(glT1,m,glT1+1,m,glT1-m+2,m,glT1-m+3,m) ;
    for (coordP1=glT1-m+2;coordP1<=m;coordP1++) {
      S8(glT1,coordP1) ;
    }
  }
  if (m >= 3) {
    S5(2*m-4,m-1,2*m-3,m-1,m-1,m-1,m,m-1) ;
    S6(2*m-4,m-2,2*m-3,m-1,m,m-2,m,m-1) ;
    S1(2*m-3,m-1,m,m-1) ;
    S3(2*m-3,m-1,2*m-2,m,m,m-1,m,m) ;
    S5(2*m-4,m,2*m-3,m,m-2,m,m-1,m) ;
    S6(2*m-4,m-1,2*m-3,m,m-1,m-1,m-1,m) ;
    S7(2*m-4,m-1,2*m-2,m,m-1,m-1,m,m) ;
    S1(2*m-3,m,m-1,m) ;
    S2(2*m-3,m,2*m-2,m,m-1,m,m,m) ;
    for (coordP1=m-1;coordP1<=m;coordP1++) {
      S8(2*m-3,coordP1) ;
    }
  }
  if (m == 2) {
    S5(0,1,1,1,1,1,2,1) ;
    S1(1,1,2,1) ;
    S3(1,1,2,2,2,1,2,2) ;
    S6(0,1,1,2,1,1,1,2) ;
    S7(0,1,2,2,1,1,2,2) ;
    S1(1,2,1,2) ;
    S2(1,2,2,2,1,2,2,2) ;
    for (coordP1=1;coordP1<=2;coordP1++) {
      S8(1,coordP1) ;
    }
  }
  if (m >= 2) {
    S5(2*m-3,m,2*m-2,m,m-1,m,m,m) ;
    S6(2*m-3,m-1,2*m-2,m,m,m-1,m,m) ;
    S1(2*m-2,m,m,m) ;
    S8(2*m-2,m) ;
  }
}
