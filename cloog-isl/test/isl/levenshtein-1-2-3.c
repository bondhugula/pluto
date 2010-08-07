/* Generated from ../../../git/cloog/test/levenshtein-1-2-3.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.12s. */
S1(0,0) ;
S2(1,0) ;
S3(1,1) ;
for (i=2;i<=N;i++) {
  S2(i,0) ;
  for (j=1;j<=i-1;j++) {
    S6(i,j) ;
  }
  S3(i,i) ;
}
S7(N+1,0) ;
for (j=1;j<=N;j++) {
  S6(N+1,j) ;
  S8(N+1,j) ;
}
for (i=N+2;i<=2*M-N-2;i++) {
  j = floord(i-N-1,2) ;
  S7(i,j) ;
  if ((i+N)%2 == 0) {
    S5(i,(i-N)/2) ;
    S8(i,(i-N)/2) ;
  }
  for (j=ceild(i-N+1,2);j<=floord(i+N-1,2);j++) {
    S6(i,j) ;
    S8(i,j) ;
  }
  if ((i+N)%2 == 0) {
    S4(i,(i+N)/2) ;
    S8(i,(i+N)/2) ;
  }
}
for (i=2*M-N-1;i<=2*M-2;i++) {
  for (j=i-M+1;j<=M-1;j++) {
    S6(i,j) ;
  }
}
