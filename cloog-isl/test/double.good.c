/* Generated from ../../../git/cloog/test/double.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i) { hash(4); hash(i); }

void test(int M, int N)
{
  /* Original iterators. */
  int i, j;
  for (i=0;i<=M;i++) {
    S1(i) ;
    for (j=0;j<=N;j++) {
      S2(i,j) ;
      S3(i,j) ;
    }
    S4(i) ;
  }
}
