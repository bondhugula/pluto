/* Generated from ../../../git/cloog/test/guide.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i) { hash(2); hash(i); }

void test(int M, int N)
{
  /* Original iterators. */
  int i;
  for (i=1;i<=N;i++) {
    if (i >= M) {
      S1(i) ;
    }
    if (i <= min(2*M,M-1)) {
      S1(i) ;
    }
  }
  for (i=N+1;i<=2*N;i++) {
    S2(i) ;
  }
}
