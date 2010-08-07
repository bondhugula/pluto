/* Generated from ../../../git/cloog/test/./reservoir/fusion1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i) { hash(2); hash(i); }
#define S3(i) { hash(3); hash(i); }

void test(int M)
{
  /* Scattering iterators. */
  int c2;
  /* Original iterators. */
  int i;
  for (c2=0;c2<=M;c2++) {
    S1(c2) ;
  }
  for (c2=1;c2<=M;c2++) {
    S2(c2) ;
  }
  for (c2=0;c2<=M;c2++) {
    S3(c2) ;
  }
}
