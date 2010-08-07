/* Generated from ../../../git/cloog/test/lu2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k,l) { hash(1); hash(i); hash(j); hash(k); hash(l); }
#define S2(i,j,k,l,m) { hash(2); hash(i); hash(j); hash(k); hash(l); hash(m); }

void test(int n)
{
  /* Original iterators. */
  int i, j, k, l, m;
  if (n >= 2) {
    for (l=2;l<=n;l++) {
      S1(1,n,1,l) ;
    }
  }
  for (i=2;i<=n-1;i++) {
    for (j=2;j<=n-1;j++) {
      for (k=1;k<=min(j-1,i-1);k++) {
        S2(i,j,k,j,i) ;
      }
    }
    for (k=1;k<=i-1;k++) {
      S2(i,n,k,n,i) ;
    }
    for (l=i+1;l<=n;l++) {
      S1(i,n,i,l) ;
    }
  }
  if (n >= 2) {
    for (j=2;j<=n;j++) {
      for (k=1;k<=j-1;k++) {
        S2(n,j,k,j,n) ;
      }
    }
  }
}
