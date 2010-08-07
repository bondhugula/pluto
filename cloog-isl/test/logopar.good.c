/* Generated from ../../../git/cloog/test/logopar.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int m, int n)
{
  /* Original iterators. */
  int i, j;
  for (j=0;j<=m;j++) {
    S1(1,j) ;
  }
  if (m >= n+1) {
    for (i=2;i<=n;i++) {
      for (j=0;j<=i-2;j++) {
        S2(i,j) ;
      }
      for (j=i-1;j<=n;j++) {
        S1(i,j) ;
        S2(i,j) ;
      }
      for (j=n+1;j<=m;j++) {
        S1(i,j) ;
      }
    }
  }
  if (m == n) {
    for (i=2;i<=n;i++) {
      for (j=0;j<=i-2;j++) {
        S2(i,j) ;
      }
      for (j=i-1;j<=n;j++) {
        S1(i,j) ;
        S2(i,j) ;
      }
    }
  }
  for (i=n+1;i<=m+1;i++) {
    for (j=i-1;j<=m;j++) {
      S1(i,j) ;
    }
  }
}
