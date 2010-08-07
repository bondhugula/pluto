/* Generated from ../../../git/cloog/test/merge.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i) { hash(2); hash(i); }
#define S3(i) { hash(3); hash(i); }

void test()
{
  /* Scattering iterators. */
  int c1;
  /* Original iterators. */
  int i;
  for (c1=0;c1<=10;c1++) {
    if (c1 == 0) {
      S1(0) ;
    }
    if (c1 >= 2) {
      S2(c1) ;
    }
    S3(c1) ;
  }
}
