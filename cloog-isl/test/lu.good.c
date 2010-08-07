/* Generated from ../../../git/cloog/test/lu.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test(int n)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j, k;
  if (n >= 2) {
    for (j=2;j<=n;j++) {
      S1(1,j) ;
    }
  }
  for (c1=2;c1<=n-1;c1++) {
    for (c2=2;c2<=n-1;c2++) {
      for (i=1;i<=min(c2-1,c1-1);i++) {
        S2(i,c2,c1) ;
      }
    }
    for (i=1;i<=c1-1;i++) {
      S2(i,n,c1) ;
    }
    for (j=c1+1;j<=n;j++) {
      S1(c1,j) ;
    }
  }
  if (n >= 2) {
    for (c2=2;c2<=n;c2++) {
      for (i=1;i<=c2-1;i++) {
        S2(i,c2,n) ;
      }
    }
  }
}
