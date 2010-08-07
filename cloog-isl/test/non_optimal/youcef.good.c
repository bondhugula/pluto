/* Generated from ../../../git/cloog/test/./non_optimal/youcef.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }

void test()
{
  /* Original iterators. */
  int i, j;
  for (i=0;i<=3;i++) {
    S1(i,i) ;
    S2(i,i) ;
    for (j=i+1;j<=4;j++) {
      S2(i,j) ;
    }
    S2(i,5) ;
    S3(i,5) ;
  }
  S1(4,4) ;
  S2(4,4) ;
  S2(4,5) ;
  S3(4,5) ;
  S1(5,5) ;
  S2(5,5) ;
  S3(5,5) ;
}
