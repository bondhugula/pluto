/* Generated from ../../../git/cloog/test/orc.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.06s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i) { hash(4); hash(i); }
#define S5(i) { hash(5); hash(i); }
#define S6(i,j) { hash(6); hash(i); hash(j); }
#define S7(i) { hash(7); hash(i); }

void test()
{
  /* Scattering iterators. */
  int p1, p2;
  /* Original iterators. */
  int i, j;
  S1(0) ;
  S2(0,0) ;
  for (p2=1;p2<=22;p2++) {
    if ((p2+1)%2 == 0) {
      j = (p2-1)/2 ;
      S3(0,(p2-1)/2) ;
    }
    if (p2%2 == 0) {
      S2(0,p2/2) ;
    }
  }
  S3(0,11) ;
  for (p1=2;p1<=6;p1++) {
    if ((p1+1)%3 == 0) {
      i = (p1-2)/3 ;
      S4((p1-2)/3) ;
    }
    if ((p1+2)%3 == 0) {
      i = (p1-1)/3 ;
      S2((p1-1)/3,0) ;
    }
    if (p1%3 == 0) {
      S1(p1/3) ;
    }
    for (p2=1;p2<=floord(-2*p1+68,3);p2++) {
      if ((p1+2)%3 == 0) {
        i = (p1-1)/3 ;
        if ((p2+1)%2 == 0) {
          j = (p2-1)/2 ;
          S3((p1-1)/3,(p2-1)/2) ;
        }
        if (p2%2 == 0) {
          S2((p1-1)/3,p2/2) ;
        }
      }
    }
    p2 = floord(-2*p1+71,3) ;
    if ((p1+2)%3 == 0) {
      i = (p1-1)/3 ;
      if ((p2+1)%2 == 0) {
        j = (p2-1)/2 ;
        S3((p1-1)/3,(p2-1)/2) ;
      }
    }
  }
  S2(2,0) ;
  for (p2=1;p2<=18;p2++) {
    if ((p2+1)%2 == 0) {
      j = (p2-1)/2 ;
      S3(2,(p2-1)/2) ;
    }
    if (p2%2 == 0) {
      S2(2,p2/2) ;
    }
  }
  S3(2,9) ;
  S4(2) ;
  S5(0) ;
  S6(0,0) ;
  for (p2=1;p2<=9;p2++) {
    S6(0,p2) ;
  }
  for (p1=2;p1<=42;p1++) {
    if ((p1+1)%3 == 0) {
      i = (p1-2)/3 ;
      S7((p1-2)/3) ;
    }
    if ((p1+2)%3 == 0) {
      i = (p1-1)/3 ;
      S6((p1-1)/3,0) ;
    }
    if (p1%3 == 0) {
      S5(p1/3) ;
    }
    for (p2=1;p2<=9;p2++) {
      if ((p1+2)%3 == 0) {
        i = (p1-1)/3 ;
        S6((p1-1)/3,p2) ;
      }
    }
  }
  S6(14,0) ;
  for (p2=1;p2<=9;p2++) {
    S6(14,p2) ;
  }
  S7(14) ;
}
