/* Generated from ../../../git/cloog/test/min-2-1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }

void test(int M, int N)
{
  /* Original iterators. */
  int i, j, k;
  if (M >= 0) {
    for (i=1;i<=N;i++) {
      for (j=0;j<=min(min(i,M),-i+N);j++) {
        for (k=0;k<=min(min(M,i),-i+N);k++) {
          S1(i,j,k) ;
        }
      }
    }
  }
}
