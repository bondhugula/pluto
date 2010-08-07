/* Generated from ../../../git/cloog/test/forwardsub-2-1-2-3.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i,j) { hash(4); hash(i); hash(j); }

void test(int M)
{
  /* Original iterators. */
  int i, j, k;
  S3(1,0) ;
  for (k=2;k<=M;k++) {
    S1(1,1,k) ;
  }
  for (i=2;i<=M-1;i++) {
    S4(i,0) ;
    for (k=i+1;k<=M;k++) {
      S2(i,1,k) ;
    }
  }
  S4(M,0) ;
}
