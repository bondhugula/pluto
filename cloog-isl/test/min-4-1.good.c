/* Generated from ../../../git/cloog/test/min-4-1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }

void test(int M, int N, int O)
{
  /* Original iterators. */
  int i;
  for (i=max(-M,-N);i<=min(O,N);i++) {
    S1(i) ;
  }
}
