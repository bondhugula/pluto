/* Generated from ../../../git/cloog/test/lux.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k,l) { hash(1); hash(i); hash(j); hash(k); hash(l); }
#define S2(i,j,k,l,m) { hash(2); hash(i); hash(j); hash(k); hash(l); hash(m); }

void test(int M)
{
  /* Original iterators. */
  int i, j, k, l, m;
  if (M >= 2) {
    for (l=2;l<=M;l++) {
      S1(1,1,M,l) ;
    }
  }
  for (i=2;i<=M-1;i++) {
    for (j=1;j<=i-1;j++) {
      for (k=j+1;k<=M;k++) {
        S2(i,j,k,k,i) ;
      }
    }
    for (l=i+1;l<=M;l++) {
      S1(i,i,M,l) ;
    }
  }
  if (M >= 2) {
    for (j=1;j<=M-1;j++) {
      for (k=j+1;k<=M;k++) {
        S2(M,j,k,k,M) ;
      }
    }
  }
}
