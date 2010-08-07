/* Generated from ../../../git/cloog/test/./reservoir/pingali5.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.03s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  for (c2=3;c2<=2*M-3;c2++) {
    for (c4=ceild(c2+3,2);c4<=M;c4++) {
      for (i=ceild(c2+1,2);i<=min(c4-1,c2-1);i++) {
        j = c2-i ;
        S1(i,c2-i,c4) ;
      }
    }
    for (c4=max(1,c2-M);c4<=floord(c2-1,2);c4++) {
      i = c2-c4 ;
      S2(c2-c4,c4) ;
    }
    for (c4=ceild(c2+3,2);c4<=M;c4++) {
      for (i=ceild(c2+1,2);i<=min(c4-1,c2-1);i++) {
        j = c2-i ;
        S3(i,c2-i,c4) ;
      }
    }
  }
  for (c2=max(2*M-2,3);c2<=2*M-1;c2++) {
    for (c4=max(1,c2-M);c4<=floord(c2-1,2);c4++) {
      i = c2-c4 ;
      S2(c2-c4,c4) ;
    }
  }
}
