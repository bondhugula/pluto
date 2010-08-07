/* Generated from ../../../git/cloog/test/iftest2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int M, int N)
{
  /* Original iterators. */
  int i, j;
  if (M >= 1) {
    for (i=1;i<=N;i++) {
      for (j=1;j<=M;j++) {
        if (i <= 2*M) {
          S1(i,j) ;
        }
        if (i >= max(M,2*M+1)) {
          S1(i,j) ;
        }
      }
    }
  }
}
