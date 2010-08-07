	#define S1(zT0,zT1,zT2,i,j,k)	c[i][k]+=a[j][k]*b[i][j];
	#define S2(zT0,zT1,zT2,i,j,k)	c[i][j]+=a[j][j]*b[i][j];
	#define S3(zT0,zT1,zT2,i,j)	c[i][j]+=a[j][j]*b[i][j];

		int t1, t2, t3, t4, t5, t6, t7;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.05s. */
for (t2=0;t2<=floord(N-1,32);t2++) {
  for (t3=0;t3<=floord(N-1,32);t3++) {
    for (t4=0;t4<=min(floord(N-3,32),t3);t4++) {
      for (t5=32*t2;t5<=min(N-1,32*t2+31);t5++) {
        for (t6=32*t4;t6<=min(min(N-3,32*t3+29),32*t4+31);t6++) {
{
	lbv=max(32*t3,t6+2); 	ubv=min(N-1,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
            S2(t2,t3,t4,t5,t7,t6);
          }
}
        }
      }
    }
  }
}
for (t2=0;t2<=floord(N-1,32);t2++) {
  for (t3=0;t3<=floord(N-1,32);t3++) {
    for (t5=32*t2;t5<=min(N-1,32*t2+31);t5++) {
{
	lbv=32*t3; 	ubv=min(N-1,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S3(t2,t3,0,t5,t7);
      }
}
    }
  }
}
for (t2=0;t2<=floord(N-1,32);t2++) {
  for (t3=0;t3<=floord(N-3,32);t3++) {
    for (t4=t3;t4<=floord(N-1,32);t4++) {
      for (t5=32*t2;t5<=min(N-1,32*t2+31);t5++) {
        for (t6=max(32*t4,32*t3+2);t6<=min(N-1,32*t4+31);t6++) {
{
	lbv=32*t3; 	ubv=min(32*t3+31,t6-2);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
            S1(t2,t3,t4,t5,t6,t7);
          }
}
        }
      }
    }
  }
}
/* End of CLooG code */
