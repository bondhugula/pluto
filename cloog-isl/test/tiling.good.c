/* Generated from ../../../git/cloog/test/tiling.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(ii,i) { hash(1); hash(ii); hash(i); }

void test(int n)
{
  /* Original iterators. */
  int ii, i;
  for (ii=0;ii<=floord(n,10);ii++) {
    for (i=max(10*ii,0);i<=min(10*ii+9,n);i++) {
      S1(ii,i) ;
    }
  }
}
