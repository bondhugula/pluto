/* Generated from ../../../git/cloog/test/yosr2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i) { hash(2); hash(i); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }
#define S4(i,j) { hash(4); hash(i); hash(j); }

void test(int M)
{
  /* Scattering iterators. */
  int proc;
  /* Original iterators. */
  int i, j, k;
  for (i=1;i<=M;i++) {
    S2(i) ;
  }
  for (proc=2;proc<=M-1;proc++) {
    for (i=1;i<=proc-1;i++) {
      S4(i,proc) ;
    }
    for (j=1;j<=proc-1;j++) {
      S1(proc,j) ;
    }
    for (j=proc+1;j<=M;j++) {
      for (k=1;k<=proc-1;k++) {
        S3(proc,j,k) ;
      }
    }
  }
  for (i=1;i<=M-1;i++) {
    S4(i,M) ;
  }
  for (j=1;j<=M-1;j++) {
    S1(M,j) ;
  }
}
