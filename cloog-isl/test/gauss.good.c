/* Generated from ../../../git/cloog/test/gauss.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test(int M)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j, k;
  if (M >= 2) {
    for (c2=2;c2<=M;c2++) {
      for (j=2;j<=M;j++) {
        S2(1,j,c2) ;
      }
    }
  }
  for (c1=2;c1<=M-1;c1++) {
    for (c2=c1+1;c2<=M;c2++) {
      for (j=1;j<=c1-1;j++) {
        S1(c1,j,c2) ;
      }
      for (j=c1+1;j<=M;j++) {
        S2(c1,j,c2) ;
      }
    }
  }
}
