/* Generated from /home/skimo/git/cloog/test/constbound.cloog by CLooG 0.14.0-170-g72daac3 64 bits in 0.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }

void test()
{
  /* Scattering iterators. */
  int t0, t2, t3;
  /* Original iterators. */
  int i, j, k;
  for (t0=0;t0<=199;t0++) {
    for (t2=max(0,50*t0);t2<=50*t0+24;t2++) {
      for (t3=0;t3<=t2;t3++) {
        S1(t0,t2,t3);
      }
    }
    for (t2=50*t0+25;t2<=min(9999,50*t0+49);t2++) {
      for (t3=0;t3<=t2;t3++) {
        S2(t0,t2,t3);
      }
    }
  }
}
