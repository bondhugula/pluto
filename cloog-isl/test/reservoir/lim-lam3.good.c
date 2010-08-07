/* Generated from ../../../git/cloog/test/./reservoir/lim-lam3.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.04s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i) { hash(4); hash(i); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  for (c2=5;c2<=min(5*M,8);c2++) {
    if (c2%5 == 0) {
      S4(c2/5) ;
    }
  }
  for (c2=9;c2<=min(13,5*M-1);c2++) {
    for (c4=max(1,ceild(c2-M-3,4));c4<=floord(c2-4,5);c4++) {
      i = c2-4*c4-3 ;
      S2(c2-4*c4-3,c4) ;
    }
    if (c2%5 == 0) {
      S4(c2/5) ;
    }
    for (c4=max(1,ceild(c2-3*M-1,2));c4<=floord(c2-4,5);c4++) {
      if ((c2+c4+2)%3 == 0) {
        i = (c2-2*c4-1)/3 ;
        S3((c2-2*c4-1)/3,c4) ;
      }
    }
  }
  for (c2=14;c2<=5*M-1;c2++) {
    for (c4=max(2,ceild(c2-M-3,4));c4<=min(M-1,floord(c2-8,3));c4++) {
      for (c6=max(1,ceild(c2-2*c4-M-5,2));c6<=min(c4-1,floord(c2-3*c4-6,2));c6++) {
        i = c2-2*c4-2*c6-5 ;
        S1(c2-2*c4-2*c6-5,c4,c6) ;
      }
    }
    for (c4=max(ceild(c2-M-3,4),1);c4<=floord(c2-4,5);c4++) {
      i = c2-4*c4-3 ;
      S2(c2-4*c4-3,c4) ;
    }
    if (c2%5 == 0) {
      S4(c2/5) ;
    }
    for (c4=max(ceild(c2-3*M-1,2),1);c4<=floord(c2-4,5);c4++) {
      if ((c2+c4+2)%3 == 0) {
        i = (c2-2*c4-1)/3 ;
        S3((c2-2*c4-1)/3,c4) ;
      }
    }
  }
  if (M >= 2) {
    c2 = 5*M ;
    S4(M) ;
  }
}
