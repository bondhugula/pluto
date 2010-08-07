/* Generated from ../../../git/cloog/test/./reservoir/fusion2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c2, c4;
  /* Original iterators. */
  int i, j;
  if ((M >= 1) && (N >= 1)) {
    for (c4=1;c4<=M;c4++) {
      S1(1,c4) ;
    }
  }
  if (M >= 1) {
    for (c2=2;c2<=N;c2++) {
      for (c4=1;c4<=M;c4++) {
        i = c2-1 ;
        S2(c2-1,c4) ;
      }
      for (c4=1;c4<=M;c4++) {
        S1(c2,c4) ;
      }
    }
  }
  if ((M >= 1) && (N >= 1)) {
    c2 = N+1 ;
    for (c4=1;c4<=M;c4++) {
      S2(N,c4) ;
    }
  }
}
