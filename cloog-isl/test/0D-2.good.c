/* Generated from ../../../git/cloog/test/0D-2.cloog by CLooG 0.14.0-117-g5444fca gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1() { hash(1); }

void test(int M)
{
  if (M >= 0) {
    S1() ;
  }
}
