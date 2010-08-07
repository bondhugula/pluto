/* Generated from ../../../git/cloog/test/yosr.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test(int n)
{
  /* Scattering iterators. */
  int proc;
  /* Original iterators. */
  int i, j, k;
  if (n >= 2) {
    for (j=2;j<=n;j++) {
      S1(1,j) ;
    }
  }
  for (proc=2;proc<=n-1;proc++) {
    for (i=1;i<=proc-1;i++) {
      for (j=i+1;j<=n;j++) {
        S2(i,j,proc) ;
      }
    }
    for (j=proc+1;j<=n;j++) {
      S1(proc,j) ;
    }
  }
  if (n >= 2) {
    for (i=1;i<=n-1;i++) {
      for (j=i+1;j<=n;j++) {
        S2(i,j,n) ;
      }
    }
  }
}
