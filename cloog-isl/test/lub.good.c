/* Generated from ../../../git/cloog/test/lub.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }
#define S4(i,j) { hash(4); hash(i); hash(j); }

void test(int M)
{
  /* Original iterators. */
  int i, j, k;
  for (i=1;i<=M-1;i++) {
    for (j=i+1;j<=M;j++) {
      S1(i,j) ;
      for (k=i+1;k<=M;k++) {
        S2(i,j,k) ;
        S3(i,j,k) ;
      }
      S4(i,j) ;
    }
  }
}
