/* Generated from ../../../git/cloog/test/4-param.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i) { hash(1); hash(i); }
#define S2(i) { hash(2); hash(i); }

void test(int m, int n, int p, int q)
{
  /* Original iterators. */
  int i;
  for (i=m;i<=min(min(n,p-1),q);i++) {
    S1(i) ;
  }
  for (i=p;i<=min(min(q,m-1),n);i++) {
    S2(i) ;
  }
  for (i=max(m,p);i<=min(n,q);i++) {
    S1(i) ;
    S2(i) ;
  }
  for (i=max(m,q+1);i<=n;i++) {
    S1(i) ;
  }
  for (i=max(p,n+1);i<=q;i++) {
    S2(i) ;
  }
}
