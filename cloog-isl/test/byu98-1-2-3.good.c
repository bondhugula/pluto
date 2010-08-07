/* Generated from ../../../git/cloog/test/byu98-1-2-3.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
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
  /* Original iterators. */
  int i, j;
  for (i=2;i<=3;i++) {
    for (j=-i+6;j<=6;j++) {
      S1(i,j) ;
    }
  }
  for (j=3;j<=4;j++) {
    S1(4,j) ;
  }
  S1(4,5) ;
  S2(4,5) ;
  S1(4,6) ;
  S1(5,4) ;
  S2(5,4) ;
  for (j=5;j<=6;j++) {
    S1(5,j) ;
  }
  for (i=6;i<=7;i++) {
    j = -i+9 ;
    S2(i,-i+9) ;
    for (j=i-1;j<=6;j++) {
      S1(i,j) ;
    }
  }
  S2(8,1) ;
}
