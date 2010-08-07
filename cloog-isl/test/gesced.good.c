/* Generated from ../../../git/cloog/test/gesced.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c1;
  /* Original iterators. */
  int i, j;
  for (c1=1;c1<=N;c1++) {
    S1(c1) ;
  }
  for (c1=N+1;c1<=2*N;c1++) {
    for (i=1;i<=N;i++) {
      j = c1-N ;
      S2(i,c1-N) ;
    }
  }
  for (c1=2*N+1;c1<=M+N;c1++) {
    for (i=1;i<=N;i++) {
      j = c1-2*N ;
      S3(i,c1-2*N) ;
      j = c1-N ;
      S2(i,c1-N) ;
    }
  }
  for (c1=M+N+1;c1<=M+2*N;c1++) {
    for (i=1;i<=N;i++) {
      j = c1-2*N ;
      S3(i,c1-2*N) ;
    }
  }
}
