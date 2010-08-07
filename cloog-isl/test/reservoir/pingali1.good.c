/* Generated from ../../../git/cloog/test/./reservoir/pingali1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if (N >= 2) {
    for (c2=1;c2<=M;c2++) {
      for (c4=1;c4<=2;c4++) {
        if ((c4+1)%2 == 0) {
          j = (c4+1)/2 ;
          S2(c2,(c4+1)/2) ;
        }
      }
      for (c4=3;c4<=2*N-1;c4++) {
        for (c6=max(1,c4-N);c6<=floord(c4-1,2);c6++) {
          j = c4-c6 ;
          S1(c2,c4-c6,c6) ;
        }
        if ((c4+1)%2 == 0) {
          j = (c4+1)/2 ;
          S2(c2,(c4+1)/2) ;
        }
      }
    }
  }
  if (N == 1) {
    for (c2=1;c2<=M;c2++) {
      S2(c2,1) ;
    }
  }
}
