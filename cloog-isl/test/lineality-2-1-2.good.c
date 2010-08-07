/* Generated from ../../../git/cloog/test/lineality-2-1-2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
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
  for (i=1;i<=M-3;i++) {
    for (j=1;j<=i+1;j++) {
      S1(i,j) ;
    }
    j = i+2 ;
    S1(i,i+2) ;
    S2(i,i+2) ;
    for (j=i+3;j<=M;j++) {
      S1(i,j) ;
    }
  }
  if (M >= 3) {
    i = M-2 ;
    for (j=1;j<=M-1;j++) {
      S1(M-2,j) ;
    }
    S1(M-2,M) ;
    S2(M-2,M) ;
  }
  for (i=M-1;i<=M;i++) {
    for (j=1;j<=M;j++) {
      S1(i,j) ;
    }
    j = i+2 ;
    S2(i,i+2) ;
  }
}
