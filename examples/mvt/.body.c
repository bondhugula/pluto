
	#define S1(zT0,zT1,zT2,zT3,i,j)	x1[i]=x1[i]+a[i][j]*y_1[j];
	#define S2(zT0,zT1,zT2,zT3,i,j)	x2[i]=x2[i]+a[j][i]*y_2[j];

		int t1, t1t, newlb_t1, newub_t1, t2, t3, t4, t5, t6, t6t, newlb_t6, newub_t6, t7;

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.02s. */
if (N >= 1) {
	lb1=0;
	ub1=floord(N-1,256);
#pragma omp parallel for shared(t1,lb1,ub1) private(t2,t3,t4,t5,t6,t7)
	for (t2=lb1; t2<=ub1; t2++) {
    for (t3=0;t3<=floord(N-1,256);t3++) {
      for (t4=8*t2;t4<=min(floord(N-1,32),8*t2+7);t4++) {
        for (t5=8*t3;t5<=min(floord(N-1,32),8*t3+7);t5++) {
          for (t6=32*t5;t6<=min(N-1,32*t5+31);t6++) {
{
	lbv=32*t4; 	ubv=min(N-1,32*t4+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
              S1(t2,t3,t4,t5,t7,t6);
              S2(t2,t3,t4,t5,t7,t6);
            }
}
          }
        }
      }
    }
  }
}
/* End of CLooG code */
