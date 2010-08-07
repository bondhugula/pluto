/* Generated from stride.cloog by CLooG 0.14.0-200-g26bdb56 gmp bits in 0.00s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test()
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j;
  for (c1=3;c1<=100;c1++) {
    if (c1 == 25) {
      S1(25);
    }
    if (c1%3 == 0) {
      S2(c1,c1/3);
    }
  }
}
