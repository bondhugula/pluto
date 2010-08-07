/* Generated from ../../../git/cloog/test/./reservoir/two.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }

void test()
{
  /* Original iterators. */
  int i, j, k;
  for (i=0;i<=1;i++) {
    if ((i+1)%2 == 0) {
      j = (-i+3)/2 ;
      k = (i+9)/2 ;
      S1(i,(-i+3)/2,(i+9)/2) ;
    }
  }
}
