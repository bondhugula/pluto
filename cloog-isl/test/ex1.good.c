/* Generated from ../../../git/cloog/test/ex1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int n)
{
  /* Original iterators. */
  int i, j;
  for (i=0;i<=14;i++) {
    for (j=0;j<=n-15;j++) {
      S1(i,j) ;
    }
  }
  for (i=15;i<=n;i++) {
    for (j=0;j<=9;j++) {
      S1(i,j) ;
    }
    for (j=10;j<=n-15;j++) {
      S1(i,j) ;
      S2(i,j) ;
    }
    for (j=n-14;j<=n;j++) {
      S2(i,j) ;
    }
  }
}
