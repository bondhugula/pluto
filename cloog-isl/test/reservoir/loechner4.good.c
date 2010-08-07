/* Generated from ../../../git/cloog/test/./reservoir/loechner4.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k,l) { hash(1); hash(i); hash(j); hash(k); hash(l); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4, c6, c8;
  /* Original iterators. */
  int i, j, k, l;
  for (c2=2;c2<=2*M;c2++) {
    for (c4=1;c4<=M;c4++) {
      for (c6=1;c6<=M;c6++) {
        for (c8=max(1,c2-M);c8<=min(c2-1,M);c8++) {
          l = c2-c8 ;
          S1(c6,c4,c8,c2-c8) ;
        }
      }
    }
  }
}
