/* Generated from ../../../git/cloog/test/./reservoir/cholesky2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.05s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  for (c2=2;c2<=min(3,3*M-4);c2++) {
    if ((c2+1)%3 == 0) {
      i = (c2+1)/3 ;
      S1((c2+1)/3) ;
    }
    for (c4=ceild(c2+4,3);c4<=min(c2,M);c4++) {
      if ((c2+c4)%2 == 0) {
        i = (c2-c4+2)/2 ;
        S2((c2-c4+2)/2,c4) ;
      }
    }
  }
  for (c2=4;c2<=3*M-4;c2++) {
    if ((c2+1)%3 == 0) {
      i = (c2+1)/3 ;
      S1((c2+1)/3) ;
    }
    for (c4=ceild(c2+2,3);c4<=min(c2-2,M);c4++) {
      for (c6=ceild(c2-c4+2,2);c6<=min(c2-c4,c4);c6++) {
        i = c2-c4-c6+1 ;
        S3(c2-c4-c6+1,c4,c6) ;
      }
    }
    for (c4=ceild(c2+4,3);c4<=min(M,c2);c4++) {
      if ((c2+c4)%2 == 0) {
        i = (c2-c4+2)/2 ;
        S2((c2-c4+2)/2,c4) ;
      }
    }
  }
  for (c2=max(2,3*M-3);c2<=min(3,3*M-2);c2++) {
    if ((c2+1)%3 == 0) {
      i = (c2+1)/3 ;
      S1((c2+1)/3) ;
    }
  }
  for (c2=max(3*M-3,4);c2<=3*M-2;c2++) {
    if ((c2+1)%3 == 0) {
      i = (c2+1)/3 ;
      S1((c2+1)/3) ;
    }
    for (c4=ceild(c2+2,3);c4<=min(M,c2-2);c4++) {
      for (c6=ceild(c2-c4+2,2);c6<=min(c2-c4,c4);c6++) {
        i = c2-c4-c6+1 ;
        S3(c2-c4-c6+1,c4,c6) ;
      }
    }
  }
  if (M >= 1) {
    c2 = 3*M-1 ;
    S1(M) ;
  }
}
