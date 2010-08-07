/* Generated from ../../../git/cloog/test/./reservoir/lim-lam4.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test(int M)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if (M >= 2) {
    S1(1,0,0) ;
  }
  for (c2=2;c2<=2*M-2;c2++) {
    for (c4=max(-M+1,-c2+1);c4<=-1;c4++) {
      for (i=max(1,c2-M+1);i<=min(c2+c4,M-1);i++) {
        j = c2+c4-i ;
        S1(i,c2+c4-i,-c4) ;
      }
      for (c6=max(-c4,c2-M+1);c6<=min(c2-1,M-1);c6++) {
        i = c2-c6 ;
        j = c4+c6 ;
        S2(c2-c6,c4+c6,c6) ;
      }
    }
    for (i=max(1,c2-M+1);i<=min(M-1,c2);i++) {
      j = c2-i ;
      S1(i,c2-i,0) ;
    }
  }
}
