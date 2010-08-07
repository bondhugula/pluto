/* Generated from ../../../git/cloog/test/./reservoir/mg-interp2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.07s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }
#define S4(i,j,k) { hash(4); hash(i); hash(j); hash(k); }

void test(int M, int N, int O, int P, int Q, int R, int S, int T, int U)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if ((M >= P+1) && (N >= Q+1)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=Q;c4<=N-1;c4++) {
        for (c6=P;c6<=M-1;c6++) {
          S1(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= Q+1)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=Q;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S2(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= P+1) && (N >= 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=P;c6<=M-1;c6++) {
          S3(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S4(c2,c4,c6) ;
        }
      }
    }
  }
}
