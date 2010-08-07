/* Generated from ../../../git/cloog/test/equality2.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.05s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) { hash(1); hash(i0); hash(i1); hash(i2); hash(i3); hash(i4); hash(i5); hash(i6); hash(i7); hash(i8); hash(i9); hash(i10); }
#define S2(i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) { hash(2); hash(i0); hash(i1); hash(i2); hash(i3); hash(i4); hash(i5); hash(i6); hash(i7); hash(i8); hash(i9); hash(i10); hash(i11); hash(i12); hash(i13); hash(i14); }

void test()
{
  /* Original iterators. */
  int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14;
  for (i0=1;i0<=10000;i0++) {
    for (i1=1000;i1<=1016;i1++) {
      for (i2=1;i2<=min(-2*i1+2033,2*i1-1999);i2++) {
        if (2*i1 == i2+1999) {
          S2(i0,i1,i2,1,i0,2*i1-1000,1,2,i0,i1-499,2*i1-1999,i0,2*i1-1999,i1-999,i1-999) ;
        }
        if (i2 == 1) {
          if (i1%2 == 0) {
            S1(i0,i1,i2,2,i0,(i1+2)/2,i1-999,i0,i1-999,(i1-998)/2,(i1-998)/2) ;
          }
        }
      }
    }
  }
}
