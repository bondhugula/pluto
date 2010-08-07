/* Generated from ../../../git/cloog/test/./reservoir/QR.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.27s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i) { hash(3); hash(i); }
#define S4(i,j) { hash(4); hash(i); hash(j); }
#define S5(i) { hash(5); hash(i); }
#define S6(i,j) { hash(6); hash(i); hash(j); }
#define S7(i,j,k) { hash(7); hash(i); hash(j); hash(k); }
#define S8(i,j) { hash(8); hash(i); hash(j); }
#define S9(i,j,k) { hash(9); hash(i); hash(j); hash(k); }
#define S10(i) { hash(10); hash(i); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if ((M <= -1) && (N >= 1)) {
    S1(0) ;
  }
  if ((M >= 0) && (N >= 1)) {
    S1(0) ;
  }
  if ((M >= 1) && (N >= 2)) {
    for (c4=0;c4<=M-1;c4++) {
      S2(0,c4) ;
    }
    S3(0) ;
    for (c4=0;c4<=M-1;c4++) {
      S4(0,c4) ;
    }
    S10(0) ;
    S1(1) ;
    S5(0) ;
  }
  if ((M <= 0) && (N >= 2)) {
    S3(0) ;
    S10(0) ;
    S1(1) ;
    S5(0) ;
  }
  if ((M >= 1) && (N == 1)) {
    for (c4=0;c4<=M-1;c4++) {
      S2(0,c4) ;
    }
    S3(0) ;
    for (c4=0;c4<=M-1;c4++) {
      S4(0,c4) ;
    }
    S10(0) ;
    S5(0) ;
  }
  if ((M <= 0) && (N == 1)) {
    S3(0) ;
    S10(0) ;
    S5(0) ;
  }
  for (c2=2;c2<=min(N-1,M);c2++) {
    for (c4=c2-1;c4<=N-1;c4++) {
      i = c2-2 ;
      S6(c2-2,c4) ;
      for (c6=c2-2;c6<=M-1;c6++) {
        i = c2-2 ;
        S7(c2-2,c4,c6) ;
      }
      i = c2-2 ;
      S8(c2-2,c4) ;
      for (c6=c2-2;c6<=M-1;c6++) {
        i = c2-2 ;
        S9(c2-2,c4,c6) ;
      }
    }
    for (c4=c2-1;c4<=M-1;c4++) {
      i = c2-1 ;
      S2(c2-1,c4) ;
    }
    i = c2-1 ;
    S3(c2-1) ;
    for (c4=c2-1;c4<=M-1;c4++) {
      i = c2-1 ;
      S4(c2-1,c4) ;
    }
    i = c2-1 ;
    S10(c2-1) ;
    S1(c2) ;
    i = c2-1 ;
    S5(c2-1) ;
  }
  if ((M >= 1) && (M <= N-2)) {
    c2 = M+1 ;
    for (c4=M;c4<=N-1;c4++) {
      i = M-1 ;
      S6(M-1,c4) ;
      c6 = M-1 ;
      i = M-1 ;
      k = M-1 ;
      S7(M-1,c4,M-1) ;
      i = M-1 ;
      S8(M-1,c4) ;
      c6 = M-1 ;
      i = M-1 ;
      k = M-1 ;
      S9(M-1,c4,M-1) ;
    }
    S3(M) ;
    S10(M) ;
    i = M+1 ;
    S1(M+1) ;
    S5(M) ;
  }
  if ((M >= N) && (N >= 2)) {
    c4 = N-1 ;
    i = N-2 ;
    j = N-1 ;
    S6(N-2,N-1) ;
    for (c6=N-2;c6<=M-1;c6++) {
      i = N-2 ;
      j = N-1 ;
      S7(N-2,N-1,c6) ;
    }
    i = N-2 ;
    j = N-1 ;
    S8(N-2,N-1) ;
    for (c6=N-2;c6<=M-1;c6++) {
      i = N-2 ;
      j = N-1 ;
      S9(N-2,N-1,c6) ;
    }
    for (c4=N-1;c4<=M-1;c4++) {
      i = N-1 ;
      S2(N-1,c4) ;
    }
    i = N-1 ;
    S3(N-1) ;
    for (c4=N-1;c4<=M-1;c4++) {
      i = N-1 ;
      S4(N-1,c4) ;
    }
    i = N-1 ;
    S10(N-1) ;
    i = N-1 ;
    S5(N-1) ;
  }
  if ((M == N-1) && (N >= 2)) {
    c4 = N-1 ;
    i = N-2 ;
    j = N-1 ;
    S6(N-2,N-1) ;
    c6 = N-2 ;
    i = N-2 ;
    j = N-1 ;
    k = N-2 ;
    S7(N-2,N-1,N-2) ;
    i = N-2 ;
    j = N-1 ;
    S8(N-2,N-1) ;
    c6 = N-2 ;
    i = N-2 ;
    j = N-1 ;
    k = N-2 ;
    S9(N-2,N-1,N-2) ;
    i = N-1 ;
    S3(N-1) ;
    i = N-1 ;
    S10(N-1) ;
    i = N-1 ;
    S5(N-1) ;
  }
  for (c2=max(M+2,2);c2<=N-1;c2++) {
    for (c4=c2-1;c4<=N-1;c4++) {
      i = c2-2 ;
      S6(c2-2,c4) ;
      i = c2-2 ;
      S8(c2-2,c4) ;
    }
    i = c2-1 ;
    S3(c2-1) ;
    i = c2-1 ;
    S10(c2-1) ;
    S1(c2) ;
    i = c2-1 ;
    S5(c2-1) ;
  }
  if ((M <= N-2) && (N >= 2)) {
    c4 = N-1 ;
    i = N-2 ;
    j = N-1 ;
    S6(N-2,N-1) ;
    i = N-2 ;
    j = N-1 ;
    S8(N-2,N-1) ;
    i = N-1 ;
    S3(N-1) ;
    i = N-1 ;
    S10(N-1) ;
    i = N-1 ;
    S5(N-1) ;
  }
}
