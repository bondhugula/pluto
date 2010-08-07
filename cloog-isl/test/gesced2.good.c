/* Generated from ../../../git/cloog/test/gesced2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.04s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }

void test(int M)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j;
  for (c1=1;c1<=4;c1++) {
    for (c2=5;c2<=M-10;c2++) {
      S1(c1,c2) ;
    }
  }
  for (c1=5;c1<=min(M-10,9);c1++) {
    for (c2=-c1+1;c2<=4;c2++) {
      i = c1+c2 ;
      S2(c1+c2,c1) ;
    }
    for (c2=5;c2<=M-10;c2++) {
      S1(c1,c2) ;
      i = c1+c2 ;
      S2(c1+c2,c1) ;
    }
    for (c2=M-9;c2<=-c1+M;c2++) {
      i = c1+c2 ;
      S2(c1+c2,c1) ;
    }
  }
  if (M >= 20) {
    for (c2=-9;c2<=4;c2++) {
      i = c2+10 ;
      S2(c2+10,10) ;
    }
    for (c2=5;c2<=M-10;c2++) {
      S1(10,c2) ;
      i = c2+10 ;
      S2(c2+10,10) ;
    }
  }
  for (c1=11;c1<=M-10;c1++) {
    for (c2=-c1+1;c2<=4;c2++) {
      i = c1+c2 ;
      S2(c1+c2,c1) ;
    }
    for (c2=5;c2<=-c1+M;c2++) {
      S1(c1,c2) ;
      i = c1+c2 ;
      S2(c1+c2,c1) ;
    }
    for (c2=-c1+M+1;c2<=M-10;c2++) {
      S1(c1,c2) ;
    }
  }
  for (c1=M-9;c1<=M;c1++) {
    for (c2=5;c2<=M-10;c2++) {
      S1(c1,c2) ;
    }
  }
}
