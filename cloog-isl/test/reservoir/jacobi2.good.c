/* Generated from ../../../git/cloog/test/./reservoir/jacobi2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4;
  /* Original iterators. */
  int i, j;
  for (c2=0;c2<=M-1;c2++) {
    for (c4=0;c4<=M-1;c4++) {
      S1(c2,c4) ;
    }
  }
}
