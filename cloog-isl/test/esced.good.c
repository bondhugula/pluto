/* Generated from ../../../git/cloog/test/esced.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int n, int m)
{
  /* Original iterators. */
  int i, j;
  if (n >= 1) {
    for (i=1;i<=m;i++) {
      S1(i) ;
      for (j=1;j<=n;j++) {
        S2(i,j) ;
      }
    }
  }
  if (n <= 0) {
    for (i=1;i<=m;i++) {
      S1(i) ;
    }
  }
}
