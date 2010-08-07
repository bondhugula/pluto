	#define S1(zT0,zT1,zT2,i,j,k)	c[j][k]+=a[i][j]*b[i][k]+b[i][j]*a[i][k];

		int t1, t2, t3, t4, t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.02s. */
for (t1=0;t1<=floord(N-1,32);t1++) {
  for (t2=t1;t2<=floord(N-1,32);t2++) {
    for (t3=0;t3<=floord(N-1,32);t3++) {
      for (t4=32*t1;t4<=min(N-1,32*t1+31);t4++) {
        for (t5=32*t3;t5<=min(N-1,32*t3+31);t5++) {
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
            S1(t1,t2,t3,t5,t4,t6);
          }
}
        }
      }
    }
  }
}
/* End of CLooG code */
