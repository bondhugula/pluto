/* Generated from ../../../git/cloog/test/rectangle.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int n)
{
  /* Scattering iterators. */
  int c1;
  /* Original iterators. */
  int i, j;
  for (c1=0;c1<=2*n;c1++) {
    for (i=max(c1-n,0);i<=min(c1,n);i++) {
      j = c1-i ;
      S1(i,c1-i) ;
    }
  }
}
