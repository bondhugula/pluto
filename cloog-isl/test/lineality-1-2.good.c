/* Generated from ../../../git/cloog/test/lineality-1-2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int M)
{
  /* Original iterators. */
  int i, j;
  S1(1,1) ;
  S2(1,1) ;
  for (j=2;j<=M;j++) {
    S1(1,j) ;
  }
  for (i=2;i<=M-1;i++) {
    for (j=1;j<=i-1;j++) {
      S1(i,j) ;
    }
    S1(i,i) ;
    S2(i,i) ;
    for (j=i+1;j<=M;j++) {
      S1(i,j) ;
    }
  }
  for (j=1;j<=M-1;j++) {
    S1(M,j) ;
  }
  S1(M,M) ;
  S2(M,M) ;
}
