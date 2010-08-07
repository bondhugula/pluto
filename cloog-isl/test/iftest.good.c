/* Generated from ../../../git/cloog/test/iftest.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }

void test(int m, int n)
{
  /* Original iterators. */
  int i;
  for (i=1;i<=n;i++) {
    if (i <= 2*m) {
      S1(i) ;
    }
    if (i >= max(m,2*m+1)) {
      S1(i) ;
    }
  }
}
