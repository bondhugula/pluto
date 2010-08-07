/* Generated from ../../../git/cloog/test/forwardsub-3-1-2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i,j) { hash(4); hash(i); hash(j); }

void test(int M)
{
  /* Original iterators. */
  int i, j;
  S3(2,1) ;
  S1(3,1) ;
  S1(4,1) ;
  S4(4,2) ;
  for (i=5;i<=M+1;i++) {
    S1(i,1) ;
    for (j=2;j<=floord(i-1,2);j++) {
      S2(i,j) ;
    }
    if (i%2 == 0) {
      S4(i,i/2) ;
    }
  }
  for (i=M+2;i<=2*M-1;i++) {
    for (j=i-M;j<=floord(i-1,2);j++) {
      S2(i,j) ;
    }
    if (i%2 == 0) {
      S4(i,i/2) ;
    }
  }
  i = 2*M ;
  S4(2*M,M) ;
}
