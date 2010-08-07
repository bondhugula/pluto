/* Generated from ../../../git/cloog/test/vasilache.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.15s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1() { hash(1); }
#define S2() { hash(2); }
#define S3() { hash(3); }
#define S4(i,j) { hash(4); hash(i); hash(j); }
#define S5(i,j) { hash(5); hash(i); hash(j); }
#define S6(i,j,k,l) { hash(6); hash(i); hash(j); hash(k); hash(l); }
#define S7(i,j,k,l) { hash(7); hash(i); hash(j); hash(k); hash(l); }
#define S8() { hash(8); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int p1, p3, p5, p7;
  /* Original iterators. */
  int i, j, k, l;
  S1() ;
  S2() ;
  for (p1=0;p1<=N-1;p1++) {
    for (p3=0;p3<=N-1;p3++) {
      S4(p1,p3) ;
      S5(p1,p3) ;
    }
  }
  for (p1=0;p1<=N-1;p1++) {
    for (p3=0;p3<=N-1;p3++) {
      for (p5=0;p5<=floord(N-1,32);p5++) {
        if (p5 >= 0) {
          p7 = 32*p5 ;
          l = 32*p5 ;
          S7(p1,p3,p5,32*p5) ;
        }
        if (p5 <= -1) {
          S7(p1,p3,p5,0) ;
        }
        for (p7=max(32*p5+1,1);p7<=min(N-1,32*p5+31);p7++) {
          l = p7-1 ;
          S6(p1,p3,p5,p7-1) ;
          S7(p1,p3,p5,p7) ;
        }
        if (p5 >= ceild(N-32,32)) {
          l = N-1 ;
          S6(p1,p3,p5,N-1) ;
        }
        if (p5 <= floord(N-33,32)) {
          p7 = 32*p5+32 ;
          l = 32*p5+31 ;
          S6(p1,p3,p5,32*p5+31) ;
        }
      }
    }
  }
  S8() ;
}
