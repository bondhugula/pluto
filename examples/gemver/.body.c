	#define S1(zT0,zT1,zT2,zT3,i,j)	B[i][j]=A[i][j]+u1[i]*v1[j]+u2[i]*v2[j];
	#define S2(zT0,zT1,zT2,zT3,i,j)	x[i]=x[i]+beta*B[j][i]*y[j];
	#define S3(i)	x[i]=x[i]+z[i];
	#define S4(i,j)	w[i]=w[i]+alpha*B[i][j]*x[j];

		int t1, t2, t3, t4, t5, t6, t7, t8, t9;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.06s. */
for (t2=0;t2<=floord(N-1,8000);t2++) {
  for (t3=0;t3<=floord(N-1,256);t3++) {
    for (t4=20*t2;t4<=min(floord(N-1,400),20*t2+19);t4++) {
      for (t5=16*t3;t5<=min(floord(N-1,16),16*t3+15);t5++) {
        for (t6=16*t5;t6<=min(N-1,16*t5+15);t6++) {
{
	lbv=400*t4; 	ubv=min(N-1,400*t4+399);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
            S1(t2,t3,t4,t5,t6,t7);
            S2(t2,t3,t4,t5,t7,t6);
          }
}
        }
      }
    }
  }
}
for (t2=0;t2<=N-1;t2++) {
  S3(t2);
}
for (t2=0;t2<=N-1;t2++) {
  for (t3=0;t3<=N-1;t3++) {
    S4(t2,t3);
  }
}
/* End of CLooG code */
