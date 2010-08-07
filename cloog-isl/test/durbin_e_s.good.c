/* Generated from ../../../git/cloog/test/durbin_e_s.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.05s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }
#define S4(i,j,k) { hash(4); hash(i); hash(j); hash(k); }
#define S5(i,j,k) { hash(5); hash(i); hash(j); hash(k); }
#define S6(i,j,k) { hash(6); hash(i); hash(j); hash(k); }
#define S7(i,j,k) { hash(7); hash(i); hash(j); hash(k); }
#define S8(i,j,k) { hash(8); hash(i); hash(j); hash(k); }

void test()
{
  /* Original iterators. */
  int i, j, k;
  S4(1,0,0) ;
  S7(1,0,0) ;
  S8(1,0,3) ;
  S2(2,-7,0) ;
  S3(2,-7,1) ;
  S6(2,-7,2) ;
  S8(2,0,3) ;
  S5(2,1,3) ;
  S2(3,-7,0) ;
  S3(3,-7,1) ;
  S3(3,-6,1) ;
  S6(3,-6,2) ;
  S8(3,0,3) ;
  for (j=1;j<=2;j++) {
    S5(3,j,3) ;
  }
  for (i=4;i<=8;i++) {
    S2(i,-7,0) ;
    S3(i,-7,1) ;
    for (j=-6;j<=i-10;j++) {
      S3(i,j,1) ;
    }
    j = i-9 ;
    S3(i,i-9,1) ;
    S6(i,i-9,2) ;
    S8(i,0,3) ;
    for (j=1;j<=i-1;j++) {
      S5(i,j,3) ;
    }
  }
  S2(9,-7,0) ;
  S3(9,-7,1) ;
  for (j=-6;j<=-1;j++) {
    S3(9,j,1) ;
  }
  S3(9,0,1) ;
  S6(9,0,2) ;
  S8(9,0,3) ;
  for (j=1;j<=8;j++) {
    S5(9,j,3) ;
  }
  S2(10,-7,0) ;
  S3(10,-7,1) ;
  for (j=-6;j<=0;j++) {
    S3(10,j,1) ;
  }
  S3(10,1,1) ;
  S6(10,1,2) ;
  S5(10,1,3) ;
  S1(10,1,4) ;
  for (j=2;j<=9;j++) {
    S5(10,j,3) ;
    S1(10,j,4) ;
  }
  S1(10,10,4) ;
}
