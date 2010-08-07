/* Generated from ../../../git/cloog/test/cholesky2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.11s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i) { hash(3); hash(i); }
#define S4(i,j) { hash(4); hash(i); hash(j); }
#define S5(i,j,k) { hash(5); hash(i); hash(j); hash(k); }
#define S6(i,j) { hash(6); hash(i); hash(j); }

void test(int M)
{
  /* Scattering iterators. */
  int c1, c2, c3;
  /* Original iterators. */
  int i, j, k;
  if (M >= 2) {
    for (c2=1;c2<=M-1;c2++) {
      S1(c2) ;
      for (c3=c2+1;c3<=M;c3++) {
        S4(c2,c3) ;
      }
    }
    S1(M) ;
  }
  if (M == 1) {
    S1(1) ;
  }
  if (M >= 2) {
    S3(1) ;
  }
  if (M >= 3) {
    S6(1,2) ;
    for (c2=3;c2<=M;c2++) {
      S6(1,c2) ;
      for (i=2;i<=c2-1;i++) {
        S5(i,c2,1) ;
      }
    }
  }
  if (M == 2) {
    S6(1,2) ;
  }
  for (c1=3;c1<=3*M-7;c1++) {
    if ((c1+2)%3 == 0) {
      i = (c1+2)/3 ;
      S3((c1+2)/3) ;
    }
    if (c1%3 == 0) {
      c2 = (c1+3)/3 ;
      i = (c1+3)/3 ;
      S2((c1+3)/3,c1/3) ;
    }
    c2 = floord(c1+6,3) ;
    if ((c1+1)%3 == 0) {
      i = (c1+1)/3 ;
      S6((c1+1)/3,c2) ;
    }
    if (c1%3 == 0) {
      S2(c2,c1/3) ;
    }
    for (c2=ceild(c1+7,3);c2<=M;c2++) {
      if ((c1+1)%3 == 0) {
        i = (c1+1)/3 ;
        S6((c1+1)/3,c2) ;
      }
      if (c1%3 == 0) {
        S2(c2,c1/3) ;
      }
      if ((c1+1)%3 == 0) {
        c3 = (c1+1)/3 ;
        for (i=ceild(c1+4,3);i<=c2-1;i++) {
          k = (c1+1)/3 ;
          S5(i,c2,(c1+1)/3) ;
        }
      }
    }
  }
  for (c1=max(3*M-6,3);c1<=3*M-4;c1++) {
    if ((c1+2)%3 == 0) {
      i = (c1+2)/3 ;
      S3((c1+2)/3) ;
    }
    if (c1%3 == 0) {
      c2 = (c1+3)/3 ;
      i = (c1+3)/3 ;
      S2((c1+3)/3,c1/3) ;
    }
    for (c2=ceild(c1+4,3);c2<=M;c2++) {
      if ((c1+1)%3 == 0) {
        i = (c1+1)/3 ;
        S6((c1+1)/3,c2) ;
      }
      if (c1%3 == 0) {
        S2(c2,c1/3) ;
      }
    }
  }
  if (M >= 2) {
    c1 = 3*M-3 ;
    j = M-1 ;
    S2(M,M-1) ;
  }
  if (M >= 1) {
    c1 = 3*M-2 ;
    S3(M) ;
  }
}
