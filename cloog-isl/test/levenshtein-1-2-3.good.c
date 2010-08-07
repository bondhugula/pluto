/* Generated from ../../../git/cloog/test/levenshtein-1-2-3.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.03s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }
#define S4(i,j) { hash(4); hash(i); hash(j); }
#define S5(i,j) { hash(5); hash(i); hash(j); }
#define S6(i,j) { hash(6); hash(i); hash(j); }
#define S7(i,j) { hash(7); hash(i); hash(j); }
#define S8(i,j) { hash(8); hash(i); hash(j); }

void test(int M, int N)
{
  /* Original iterators. */
  int i, j;
  S1(0,0) ;
  S2(1,0) ;
  S3(1,1) ;
  for (i=2;i<=N;i++) {
    S2(i,0) ;
    for (j=1;j<=i-1;j++) {
      S6(i,j) ;
    }
    S3(i,i) ;
  }
  i = N+1 ;
  S7(N+1,0) ;
  for (j=1;j<=N;j++) {
    S6(N+1,j) ;
    S8(N+1,j) ;
  }
  for (i=N+2;i<=2*M-N-2;i++) {
    j = floord(i-N-1,2) ;
    S7(i,j) ;
    if ((i+N)%2 == 0) {
      j = (i-N)/2 ;
      S5(i,(i-N)/2) ;
      S8(i,(i-N)/2) ;
    }
    for (j=ceild(i-N+1,2);j<=floord(i+N-1,2);j++) {
      S6(i,j) ;
      S8(i,j) ;
    }
    if ((i+N)%2 == 0) {
      j = (i+N)/2 ;
      S4(i,(i+N)/2) ;
      S8(i,(i+N)/2) ;
    }
  }
  for (i=2*M-N-1;i<=2*M-2;i++) {
    for (j=i-M+1;j<=M-1;j++) {
      S6(i,j) ;
    }
  }
}
