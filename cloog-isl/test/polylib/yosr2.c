/* Generated from ../../../git/cloog/test/yosr2.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.02s. */
for (i=1;i<=M;i++) {
  S2(i) ;
}
for (proc=2;proc<=M-1;proc++) {
  for (i=1;i<=proc-1;i++) {
    S4(i,proc) ;
  }
  for (j=1;j<=proc-1;j++) {
    S1(proc,j) ;
  }
  for (j=proc+1;j<=M;j++) {
    for (k=1;k<=proc-1;k++) {
      S3(proc,j,k) ;
    }
  }
}
for (i=1;i<=M-1;i++) {
  S4(i,M) ;
}
for (j=1;j<=M-1;j++) {
  S1(M,j) ;
}
