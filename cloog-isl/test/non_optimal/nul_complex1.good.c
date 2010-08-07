/* Generated from ../../../git/cloog/test/./non_optimal/nul_complex1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
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
  int c1, c2;
  /* Original iterators. */
  int i, j;
  for (c1=0;c1<=5*n;c1++) {
    for (c2=max(c1-n,ceild(2*c1,3));c2<=min(c1,floord(2*c1+2*n,3));c2++) {
      if (c2%2 == 0) {
        i = (-2*c1+3*c2)/2 ;
        j = c1-c2 ;
        S1((-2*c1+3*c2)/2,c1-c2) ;
      }
    }
  }
}
