/* Generated from ../../../git/cloog/test/cholesky.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.05s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i) { hash(3); hash(i); }
#define S4(i,j) { hash(4); hash(i); hash(j); }
#define S5(i,j,k) { hash(5); hash(i); hash(j); hash(k); }
#define S6(i,j) { hash(6); hash(i); hash(j); }

void test(int n)
{
  /* Scattering iterators. */
  int c1, c3, c5;
  /* Original iterators. */
  int i, j, k;
  if (n >= 2) {
    S1(1) ;
    S3(1) ;
    for (c3=2;c3<=n;c3++) {
      S4(1,c3) ;
      S6(1,c3) ;
    }
  }
  if (n == 1) {
    S1(1) ;
    S3(1) ;
  }
  for (c1=2;c1<=n-1;c1++) {
    S1(c1) ;
    for (c3=1;c3<=c1-1;c3++) {
      S2(c1,c3) ;
    }
    S3(c1) ;
    for (c3=c1+1;c3<=n;c3++) {
      S4(c1,c3) ;
      for (c5=1;c5<=c1-1;c5++) {
        S5(c1,c3,c5) ;
      }
      S6(c1,c3) ;
    }
  }
  if (n >= 2) {
    S1(n) ;
    for (c3=1;c3<=n-1;c3++) {
      S2(n,c3) ;
    }
    S3(n) ;
  }
}
