/* Generated from ../../../git/cloog/test/no_lindep.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i;
  c1 = M+1 ;
  i = N+2 ;
  S1(N+2) ;
}
