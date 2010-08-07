/* Generated from ../../../git/cloog/test/wavefront.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int n, int m)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j;
  if ((n >= 1) && (m >= 1)) {
    for (c1=2;c1<=n+m;c1++) {
      for (c2=max(1,c1-m);c2<=min(n,c1-1);c2++) {
        j = c1-c2 ;
        S1(c2,c1-c2) ;
      }
    }
  }
}
