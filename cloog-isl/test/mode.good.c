/* Generated from ../../../git/cloog/test/mode.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int M, int N)
{
  /* Original iterators. */
  int i, j;
  for (i=0;i<=min(M,N-1);i++) {
    for (j=0;j<=i;j++) {
      S1(i,j) ;
      S2(i,j) ;
    }
    for (j=i+1;j<=N;j++) {
      S2(i,j) ;
    }
  }
  if ((M >= N) && (N >= 0)) {
    for (j=0;j<=N;j++) {
      S1(N,j) ;
      S2(N,j) ;
    }
  }
  if (N >= 0) {
    for (i=N+1;i<=M;i++) {
      for (j=0;j<=N;j++) {
        S1(i,j) ;
        S2(i,j) ;
      }
      for (j=N+1;j<=i;j++) {
        S1(i,j) ;
      }
    }
  }
  if (N <= -1) {
    for (i=0;i<=M;i++) {
      for (j=0;j<=i;j++) {
        S1(i,j) ;
      }
    }
  }
}
