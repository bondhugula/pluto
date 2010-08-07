/* Generated from ../../../git/cloog/test/gesced3.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i) { hash(2); hash(i); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c1;
  /* Original iterators. */
  int i;
  for (c1=M+1;c1<=2*M;c1++) {
    i = c1-M ;
    S1(c1-M) ;
  }
  for (c1=2*M+1;c1<=M+N;c1++) {
    i = c1-2*M ;
    S2(c1-2*M) ;
    i = c1-M ;
    S1(c1-M) ;
  }
  for (c1=M+N+1;c1<=2*M+N;c1++) {
    i = c1-2*M ;
    S2(c1-2*M) ;
  }
}
