/* Generated from ../../../git/cloog/test/min-3-1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int M)
{
  /* Original iterators. */
  int i, j;
  for (i=0;i<=min(10,M);i++) {
    for (j=0;j<=min(10,M);j++) {
      S1(i,j) ;
    }
  }
}
