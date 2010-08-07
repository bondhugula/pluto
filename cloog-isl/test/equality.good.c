/* Generated from ../../../git/cloog/test/equality.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i0,i1) { hash(1); hash(i0); hash(i1); }
#define S2(i0,i1) { hash(2); hash(i0); hash(i1); }

void test()
{
  /* Original iterators. */
  int i0, i1;
  for (i0=0;i0<=5;i0++) {
    for (i1=ceild(4*i0,5);i1<=floord(6*i0+20,5);i1++) {
      if (2*i0 == i1) {
        S1(i0,i1) ;
      }
      if (i1 == 4) {
        S2(i0,i1) ;
      }
    }
  }
}
