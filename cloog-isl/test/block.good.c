/* Generated from /home/skimo/git/cloog/test/block.cloog by CLooG 0.14.0-170-g72daac3 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1() { hash(1); }
#define S2() { hash(2); }
#define S3(i) { hash(3); hash(i); }

void test()
{
  /* Scattering iterators. */
  int c1;
  /* Original iterators. */
  int i;
  S1();
  S3(0);
  S2();
  S3(1);
}
