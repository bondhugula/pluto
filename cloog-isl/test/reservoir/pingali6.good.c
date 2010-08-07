/* Generated from ../../../git/cloog/test/./reservoir/pingali6.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.03s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if (N >= 3) {
    for (c4=2;c4<=N-1;c4++) {
      for (c6=2;c6<=N-1;c6++) {
        S1(1,c4,c6) ;
      }
    }
  }
  if (N >= 3) {
    for (c2=3;c2<=2*M;c2++) {
      for (c4=2;c4<=N-1;c4++) {
        for (c6=2;c6<=N-1;c6++) {
          if (c2%2 == 0) {
            S1(c2/2,c4,c6) ;
          }
        }
      }
      for (c4=2;c4<=N-1;c4++) {
        for (c6=2;c6<=N-1;c6++) {
          if ((c2+1)%2 == 0) {
            i = (c2-1)/2 ;
            S2((c2-1)/2,c4,c6) ;
          }
        }
      }
    }
  }
  if (N >= 3) {
    c2 = 2*M+1 ;
    for (c4=2;c4<=N-1;c4++) {
      for (c6=2;c6<=N-1;c6++) {
        S2(M,c4,c6) ;
      }
    }
  }
}
