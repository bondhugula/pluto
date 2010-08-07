/* Generated from ../../../git/cloog/test/./reservoir/tang-xue1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k,l) { hash(1); hash(i); hash(j); hash(k); hash(l); }

void test()
{
  /* Scattering iterators. */
  int c2, c4, c6, c8;
  /* Original iterators. */
  int i, j, k, l;
  for (c2=0;c2<=9;c2++) {
    for (c4=max(-1,c2-9);c4<=min(4,c2+3);c4++) {
      for (c6=max(max(c2,1),c2-c4);c6<=min(min(c2+1,9),c2-c4+4);c6++) {
        for (c8=max(1,-c2+c4+c6);c8<=min(4,-c2+c4+c6+1);c8++) {
          if (c2%2 == 0) {
            if ((c2+c4)%2 == 0) {
              j = (-c2+c4)/2 ;
              k = -c2+c6 ;
              l = -c4+c8 ;
              S1(c2/2,(-c2+c4)/2,-c2+c6,-c4+c8) ;
            }
          }
        }
      }
    }
  }
}
