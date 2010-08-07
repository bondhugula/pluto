/* Generated from ../../../git/cloog/test/./reservoir/lim-lam1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.02s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test()
{
  /* Scattering iterators. */
  int c2, c4;
  /* Original iterators. */
  int i, j;
  S1(1,100) ;
  for (c2=-98;c2<=0;c2++) {
    c4 = -2*c2+2 ;
    j = -c2+1 ;
    S1(1,-c2+1) ;
    for (c4=-2*c2+3;c4<=199;c4++) {
      if (c4%2 == 0) {
        i = (2*c2+c4)/2 ;
        S1((2*c2+c4)/2,c4/2) ;
      }
      if ((c4+1)%2 == 0) {
        i = (2*c2+c4-1)/2 ;
        j = (c4+1)/2 ;
        S2((2*c2+c4-1)/2,(c4+1)/2) ;
      }
    }
    i = c2+100 ;
    S1(c2+100,100) ;
  }
  for (c2=1;c2<=99;c2++) {
    S2(c2,1) ;
    for (c4=2;c4<=-2*c2+200;c4++) {
      if (c4%2 == 0) {
        i = (2*c2+c4)/2 ;
        S1((2*c2+c4)/2,c4/2) ;
      }
      if ((c4+1)%2 == 0) {
        i = (2*c2+c4-1)/2 ;
        j = (c4+1)/2 ;
        S2((2*c2+c4-1)/2,(c4+1)/2) ;
      }
    }
    c4 = -2*c2+201 ;
    j = -c2+101 ;
    S2(100,-c2+101) ;
  }
  S2(100,1) ;
}
