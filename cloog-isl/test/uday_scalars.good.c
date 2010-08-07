/* Generated from ../../../git/cloog/test/uday_scalars.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(j,l,m) { hash(1); hash(j); hash(l); hash(m); }
#define S2(j,l,m) { hash(2); hash(j); hash(l); hash(m); }

void test(int n)
{
  /* Scattering iterators. */
  int p3;
  /* Original iterators. */
  int j, l, m;
  for (p3=0;p3<=n;p3++) {
    S1(p3,0,0) ;
  }
  for (p3=0;p3<=n;p3++) {
    S2(0,p3,0) ;
  }
}
